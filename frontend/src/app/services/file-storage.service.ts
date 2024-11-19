/*
  TODO: pagination, downloading files/folders

  PLEASE NOTE that this service currently loads all uploaded files from Firestore
  into memory. Some sort of pagination system is needed, although using the
  Firestore pagination API is debatable since it only allows to go to pages
  immediately back and forward, no page jumping. For now, sorting filters
  and such are done on the frontend.

  An idea is to load in staggered amounts of data at once which is sorted/filtered
  on the frontend. Once the user reaches the end of the loaded data, load another
  set so it's kind of like staggered pagination. Or, we could just paginate API for
  everything but this wouldn't let you go forward 2 pages at once, for example.
  Gmail is like this.

  Separate service for retrieving user files from Firestore. Includes retrieving
  files, deleting files, downloading files (TODO) and making folders.
*/

import { Injectable } from '@angular/core';
import {
  Firestore,
  OrderByDirection,
  collection,
  collectionData,
} from '@angular/fire/firestore';
import { Functions, httpsCallable } from '@angular/fire/functions';
import {
  Storage,
  getBlob,
  getDownloadURL,
  listAll,
  ref,
} from '@angular/fire/storage';
import {
  BehaviorSubject,
  Observable,
  combineLatest,
  map,
  of,
  switchMap,
} from 'rxjs';

import { saveAs } from 'file-saver';
import JSZip from 'jszip';
import { posix } from 'path-browserify';

import { UserService } from '../services/user.service';

import { ExtensionType, ExtensionTypeMap } from '../models/extension-type';
import { User } from '../models/user';
import { UserFile } from '../models/user-file';

@Injectable({
  providedIn: 'root',
})
export class FileStorageService {
  private files$: Observable<UserFile[]>;

  private searchFilter$ = new BehaviorSubject<string>('');
  private pathFilter$ = new BehaviorSubject<string[]>(['/']);
  private typeFilter$ = new BehaviorSubject<ExtensionType | null>(null);
  private sortFilter$ = new BehaviorSubject<keyof UserFile>('uploadedOn');
  private sortDirectionFilter$ = new BehaviorSubject<OrderByDirection>('desc');

  private pageSize$ = new BehaviorSubject<number>(10);
  private pageNumber$ = new BehaviorSubject<number>(1);
  private totalPages$ = new BehaviorSubject<number>(1);

  constructor(
    private firestore: Firestore,
    private functions: Functions,
    private storage: Storage,
    private userService: UserService
  ) {
    this.files$ = this.userService.getCurrentUser().pipe(
      switchMap((user: User | null) => {
        if (user) {
          const collectionRef = collection(
            this.firestore,
            `users/${user.id}/files`
          );
          return (collectionData(collectionRef) as Observable<UserFile[]>).pipe(
            map((files: UserFile[]) =>
              files.map((file: any) => {
                return {
                  ...file,
                  type:
                    ExtensionTypeMap[
                      file.extension as keyof typeof ExtensionTypeMap
                    ] || 'unknown',
                  uploadedOn: file.uploadedOn.toDate(),
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

  /**
   * Retrieves the user files as an Observable array. The search, path, type,
   * page, and sort filters are automatically applied.
   *
   * @returns An Observable array of the user's files with filters applied.
   */
  getFiles(): Observable<UserFile[]> {
    return combineLatest([
      this.files$,
      this.searchFilter$,
      this.pathFilter$,
      this.typeFilter$,
      this.sortFilter$,
      this.sortDirectionFilter$,
      this.pageSize$,
      this.pageNumber$,
    ]).pipe(
      map(
        ([
          files,
          search,
          path,
          type,
          sort,
          sortDirection,
          pageSize,
          pageNumber,
        ]) => {
          let filterFunc: (file: UserFile) => boolean;

          // file name should start with the search query, non-case-sensitive
          if (search.length > 0) {
            filterFunc = (f) =>
              f.name.toLowerCase().startsWith(search.toLowerCase());
          }
          // if no search filter, filter based on path
          else {
            filterFunc = (f) => f.path === posix.join(...path);
          }
          // finally, apply type filter
          if (type !== null) {
            filterFunc = (f) => filterFunc(f) && f.type === type;
          }

          let sortFunc: (a: UserFile, b: UserFile) => number;

          if (sortDirection === 'asc') {
            sortFunc = (a, b) => {
              const valueA: any = a[sort];
              const valueB: any = b[sort];
              return typeof valueA === 'string'
                ? valueA.localeCompare(valueB)
                : valueA - valueB;
            };
          }
          // sort by desc instead
          else {
            sortFunc = (a, b) => {
              const valueA: any = a[sort];
              const valueB: any = b[sort];
              return typeof valueA === 'string'
                ? valueB.localeCompare(valueA)
                : valueB - valueA;
            };
          }

          const filter = files.filter(filterFunc).sort(sortFunc);

          this.totalPages$.next(
            Math.max(1, Math.ceil(filter.length / pageSize))
          );

          const pageIndexStart = pageSize * (pageNumber - 1);
          const pageIndexLast = pageIndexStart + pageSize;
          return filter.slice(pageIndexStart, pageIndexLast);
        }
      )
    );
  }

  /**
   * Retrieves the current search filter as an Observable.
   *
   * @returns An Observable of the search filter.
   */
  getSearchFilter(): Observable<string> {
    return this.searchFilter$.asObservable();
  }

  /**
   * Sets the search filter.
   *
   * @param search The search filter to apply.
   */
  setSearchFilter(search: string): void {
    this.searchFilter$.next(search);
    this.pageNumber$.next(1);
  }

  /**
   * Retrieves the current path filter as on Observable.
   * Initializes to`['/']`
   *
   * @returns An Observable of the path filter.
   */
  getPathFilter(): Observable<string[]> {
    return this.pathFilter$.asObservable();
  }

  /**
   * Sets the path filter. Set to`['/']`at minimum.
   *
   * @param path The path filter to apply.
   */
  setPathFilter(path: string[]): void {
    this.pathFilter$.next(path);
    this.pageNumber$.next(1);
  }

  /**
   * Retrieves the current type filter as on Observable.
   *
   * @returns An Observable of the type filter.
   */
  getTypeFilter(): Observable<ExtensionType | null> {
    return this.typeFilter$.asObservable();
  }

  /**
   * Sets the type filter. Setting to the same value
   * will revert the type filter back to null.
   *
   * @param type The type filter to apply.
   */
  setTypeFilter(type: ExtensionType): void {
    type === this.typeFilter$.value
      ? this.typeFilter$.next(null)
      : this.typeFilter$.next(type);

    this.pageNumber$.next(1);
  }

  /**
   * Retrieves the current sort filter as on Observable.
   *
   * @returns An Observable of the sort filter.
   */
  getSortFilter(): Observable<keyof UserFile> {
    return this.sortFilter$.asObservable();
  }

  /**
   * Retrieves the current sort direction filter as on Observable.
   *
   * @returns An Observable of the sort direction filter.
   */
  getSortDirectionFilter(): Observable<OrderByDirection> {
    return this.sortDirectionFilter$.asObservable();
  }

  /**
   * Sets the sort filter. Setting to the same value will reverse
   * the sort filter's direction.
   *
   * @param sort The sort filter to apply.
   */
  setSortFilter(sort: keyof UserFile): void {
    sort == this.sortFilter$.value
      ? this.sortDirectionFilter$.value === 'asc'
        ? this.sortDirectionFilter$.next('desc')
        : this.sortDirectionFilter$.next('asc')
      : (this.sortFilter$.next(sort), this.sortDirectionFilter$.next('asc'));

    this.pageNumber$.next(1);
  }

  /**
   * Retrieves the current size of each page as an Observable.
   *
   * @returns An Observable of the current page size.
   */
  getPageSize(): Observable<number> {
    return this.pageSize$.asObservable();
  }

  /**
   * Sets the size of each page. Set to 1 at minimum.
   *
   * @param size The page size to apply.
   */
  setPageSize(size: number): void {
    if (size <= 0) {
      return;
    }

    this.pageSize$.next(size);
  }

  /**
   * Retrieves the current page number as an Observable.
   *
   * @returns An Observable of the current page number.
   */
  getPageNumber(): Observable<number> {
    return this.pageNumber$.asObservable();
  }

  /**
   * Sets the page number. Set between`[1, totalPages]`
   *
   * @param page The page number to apply.
   */
  setPageNumber(page: number): void {
    if (page <= 0 || page > this.totalPages$.value) {
      return;
    }

    this.pageNumber$.next(page);
  }

  /**
   * Retrieves the total pages available as an Observable.
   *
   * @returns An Observable of the total pages.
   */
  getTotalPages(): Observable<number> {
    return this.totalPages$.asObservable();
  }

  /**
   * Deletes all user files under a specified path. Can be used
   * to delete entire folders. Set`['/']` at minimum.
   *
   * @param path The path to delete.
   */
  async deletePath(path: string[]) {
    const cloudFunctionCallable = httpsCallable(
      this.functions,
      'request_user_delete_path'
    );
    await cloudFunctionCallable({ path: posix.join(...path) });
  }

  /**
   * Retrieves the download URL of a user file from Cloud Storage.
   *
   * @param file The file to retrieve.
   */
  async getDownloadUrl(file: UserFile): Promise<string> {
    const storageRef = ref(this.storage, file.storageLink);
    return await getDownloadURL(storageRef);
  }

  /**
   * Retrieves the blob of a user file from Cloud Storage.
   *
   * @param file The file to retrieve.
   */
  async getFileBlob(file: UserFile): Promise<Blob> {
    const storageRef = ref(this.storage, file.storageLink);
    return await getBlob(storageRef);
  }

  // TODO: these functions don't work. downloadFile() seems
  // to work occasionally (sometimes downloads, other times opens a new
  // window with file contents). since the saveAs package says it's meant
  // for locally generated files there's probably a better solution I should research.

  // for downloadFolder, I need to configure the CORS stuff on the bucket
  // to be able to use getBlob() but i'm not sure what direction we should go
  // on this at the moment. should we download directly from cloud storage?
  // or move this out to backend?

  private async downloadFile(file: UserFile) {
    // const storageRef = ref(this.storage, file.storageLink);
    // const downloadUrl = await getDownloadURL(storageRef);
    // saveAs(downloadUrl, file.name);
  }

  private async downloadFolder(folder: UserFile) {
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
}
