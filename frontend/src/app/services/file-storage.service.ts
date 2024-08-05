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
import {
  Firestore,
  OrderByDirection,
  QueryConstraint,
  orderBy,
  where,
} from '@angular/fire/firestore';
import { Functions, httpsCallable } from '@angular/fire/functions';
import {
  Storage,
  getBlob,
  getDownloadURL,
  listAll,
  ref,
} from '@angular/fire/storage';
import { Observable, map, of, switchMap } from 'rxjs';

// import { saveAs } from 'file-saver';
import JSZip from 'jszip';

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
  private pageSize = 10;

  constructor(
    private firestore: Firestore,
    private functions: Functions,
    private storage: Storage,
    private userService: UserService
  ) {
    this.files$ = this.userService.getCurrentUser().pipe(
      switchMap((user: User | null) => {
        if (user) {
          this.paginator = new FirestorePaginator<UserFile>(
            this.firestore,
            `users/${user.id}/files`,
            this.pageSize,
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

  // Sets pagination size to pageSize + 1 to check if there's a next page
  setPageSize(pageSize: number) {
    // in case the paginator doesn't exist yet (UserService is null)
    // just set the instance variable so it can read from it later.
    this.pageSize = pageSize;
    this.paginator?.setPageSize(pageSize);
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

  deletePath(path: string) {
    const cloudFunctionCallable = httpsCallable(
      this.functions,
      'request_user_delete_path'
    );
    cloudFunctionCallable({ path });
  }

  createFolder(name: string, path: string) {
    if (name == '') {
      return;
    }

    const cloudFunctionCallable = httpsCallable(
      this.functions,
      'request_user_create_folder'
    );
    cloudFunctionCallable({ name, path });
  }

  // TODO: these functions don't work. downloadFile() seems
  // to work occasionally (sometimes downloads, other times opens a new
  // window with file contents). since the saveAs package says it's meant
  // for locally generated files there's probably a better solution I should research.

  // for downloadFolder, I need to configure the CORS stuff on the bucket
  // to be able to use getBlob() but i'm not sure what direction we should go
  // on this at the moment. should we download directly from cloud storage?
  // or move this out to backend?

  async downloadFile(file: UserFile) {
    // const storageRef = ref(this.storage, file.storageLink);
    // const downloadUrl = await getDownloadURL(storageRef);
    // saveAs(downloadUrl, file.name);
  }

  async downloadFolder(folder: UserFile) {
    // const zip = new JSZip();
    // await this.addFilesFromDirectoryToZip(folder.storageLink, zip);
    // return await zip.generateAsync({ type: 'blob' });
  }

  private async addFilesFromDirectoryToZip(
    cloudStoragePath: string,
    zip: JSZip
  ) {
    // const storageRef = ref(this.storage, cloudStoragePath);
    // const directoryContents = await listAll(storageRef);
    // for (const file of directoryContents.items) {
    //   const fileRef = ref(this.storage, file.fullPath);
    //   const fileBlob = await getBlob(fileRef);
    //   zip.file(file.fullPath, fileBlob);
    // }
    // for (const folder of directoryContents.prefixes) {
    //   await this.addFilesFromDirectoryToZip(folder.fullPath, zip);
    // }
  }

  // Paginator doesn't seem to work. Previous and next methods don't show the proper items.
  // Since the existing paginator is limited by only going next and before anyway,
  // a new pagination system might be in order.
  nextPage() {
    // this.paginator?.next();
  }

  previousPage() {
    // this.paginator?.previous();
  }
}
