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
import { Firestore, collection, collectionData } from '@angular/fire/firestore';
import { Functions, httpsCallable } from '@angular/fire/functions';
import {
  Storage,
  getBlob,
  getDownloadURL,
  listAll,
  ref,
} from '@angular/fire/storage';
import { Observable, map, of, switchMap } from 'rxjs';

import { saveAs } from 'file-saver';
import JSZip from 'jszip';

import { UserService } from '../services/user.service';

import { ExtensionTypeMap } from '../models/extension-type';
import { User } from '../models/user';
import { UserFile } from '../models/user-file';

@Injectable({
  providedIn: 'root',
})
export class FileStorageService {
  private files$: Observable<UserFile[]>;

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
                  type: ExtensionTypeMap[file.extension as keyof typeof ExtensionTypeMap] || 'unknown',
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

  getFiles(): Observable<UserFile[] | null> {
    return this.files$;
  }

  async deletePath(path: string) {
    const cloudFunctionCallable = httpsCallable(
      this.functions,
      'request_user_delete_path'
    );
    await cloudFunctionCallable({ path });
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
}
