/*
  PLEASE NOTE that this service uses a utility method in utils/firestore-paginator.
  The utility class uses Firestore's pagination API but to my knowledge it only
  works on single queries at a time. This means you can't have multiple queries
  in a row and smush them together to aggregate some data.

  This introduces potential issues later on when we scale the software up.
  Currently I am storing the "type" data (dataset, code, text) directly in the
  documents within the collection, but this is static and changing them will
  require updating the entire database.

  A solution would be to place the extension mappings (".csv" => dataset) in some
  other container like a Firestore document, but this would make it hard to work
  with the Firestore pagination API.

  Therefore, later on if we scale this up we will need to rethink how the pagination/
  file retrieval system will work. For now, it is all based on single queries and
  static data within the Firestore documents.

  Although, once we decide on the categories they likely won't change for a while,
  so this is really not as big an issue as it seems. Worst case, you run a script to
  update the entire Firestore DB if you want to update the mappings. But still not
  very scalable.
*/

import { Injectable } from '@angular/core';
import { Auth } from '@angular/fire/auth';
import {
  Firestore,
  OrderByDirection,
  QueryConstraint,
  orderBy,
  where,
} from '@angular/fire/firestore';
import { Functions, httpsCallable } from '@angular/fire/functions';
import { Storage } from '@angular/fire/storage';
import {
  Observable,
  map,
  of,
  shareReplay,
  switchMap,
} from 'rxjs';

import { FirestorePaginator } from '../utils/firestore-paginator';
import { lex_next_string } from '../utils/lex-next-string';

import { UserService } from '../services/user.service';

import { User } from '../models/user';
import { UserFile, UserFileSortKey, UserFileType } from '../models/user-file';


@Injectable({
  providedIn: 'root',
})
export class FileStorageService {
  private files$: Observable<UserFile[]>;
  private paginator?: FirestorePaginator<UserFile>;

  // path does not apply to firestore path. only used for queries to documents
  // inside the uploads/{userID}/files collection to "simulate" file system
  private path = '/'; // '/' means root dir
  private searchFilter = '';
  private typeFilter: UserFileType | null = null;
  private sortField: UserFileSortKey = 'name';
  private sortDirection: OrderByDirection = 'asc';

  constructor(
    private auth: Auth,
    private firestore: Firestore,
    private functions: Functions,
    private storage: Storage,
    private userService: UserService
  ) {
    const defaultPageSize = 10;

    this.files$ = this.userService.getCurrentUser().pipe(
      switchMap((user: User | null) => {
        if (user) {
          this.paginator = new FirestorePaginator<UserFile>(
            firestore,
            `users/${user.id}/files`,
            defaultPageSize,
            this.buildQuery(),
            false,
            true
          );

          return this.paginator.getItems().pipe(
            map((files: UserFile[]) =>
              files.map((file: any) => {
                return {
                  ...file,
                  uploadedOn: file.uploadedOn?.toDate(),
                } as UserFile;
              })
            )
          );
        } else {
          return of([]);
        }
      })
    );
  }

  /*
    Builds the query and sends it to this.paginator in order to update the
    observable.

    Due to Firestore's composite indexing, combinations of queries are important
    and each new one requires creating a new composite index. For now, the
    possible query field combinations are:

    1) path (no search query) + sort by (name/type/uploaded) asc./desc.
    2) type filter + sort by (name/type/uploaded) asc./desc.
    3) name filter + sort by (name/type/uploaded) asc./desc.
    4) type filter + name filter + sort by (name/type/uploaded) asc./desc.
  */
  private buildQuery(): QueryConstraint[] {
    let queryConstraints: QueryConstraint[] = [];

    // if there is no search filter or type filter, we should only display
    // the folders in the current directory. otherwise, no constraint on path.
    if (this.searchFilter === '' && this.typeFilter === null) {
      queryConstraints.push(where('path', '==', this.path));
    }

    // query for files with the desired type
    if (this.typeFilter !== null) {
      queryConstraints.push(where('type', '==', this.typeFilter));
    }

    // query for files that start with search filter
    if (this.searchFilter !== '') {
      queryConstraints.push(
        where('name', '>=', this.searchFilter),
        where('name', '<', lex_next_string(this.searchFilter))
      );
    }

    // sort accordingly.
    queryConstraints.push(orderBy(this.sortField, this.sortDirection));
    return queryConstraints;
  }

  /*
    Builds the QueryConstraint[] based on the current settings and
    sends it to the paginator. This will send query requests to Firestore
    and trigger the Observable to change.
  */
  private updateQuery(): void {
    this.paginator?.setQueries(this.buildQuery());
  }

  getFiles(): Observable<UserFile[] | null> {
    return this.files$;
  }

  setPageSize(pageSize: number) {
    this.paginator!.setPageSize(pageSize);
  }

  setPath(path: string) {
    this.path = path;
    this.updateQuery();
  }

  setSearchFilter(search: string) {
    this.searchFilter = search;
    this.updateQuery();
  }

  setTypeFilter(type: UserFileType | null) {
    this.typeFilter = type;
    this.updateQuery();
  }

  setSortField(field: UserFileSortKey, direction: OrderByDirection) {
    this.sortField = field;
    this.sortDirection = direction;
    this.updateQuery();
  }

  downloadFile(file: UserFile) {}

  deletePath(path: string) {}
}
