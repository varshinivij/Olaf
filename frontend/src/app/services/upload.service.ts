// TODO:
// 1) upload entire folders at once
// 2) check for dupe file names
//    and rename them like "file1.txt -> file1(1).txt or similar."

// Separate service for uploading files/folders specifically for user file
// storage system, exposes an "upload queue" w/ real-time upload progress.
// Files are uploaded to the user's Cloud Storage bucket at uploads/{userID}/

// (Currently can only upload files as once, not entire folders.
// How to incorporate folders into the upload queue to show progress bar
// that makes sense?)

import { Injectable } from '@angular/core';
import { Auth } from '@angular/fire/auth';
import { Functions, httpsCallable } from '@angular/fire/functions';
import {
  Storage,
  getDownloadURL,
  ref,
  uploadBytesResumable,
} from '@angular/fire/storage';
import { BehaviorSubject, Observable } from 'rxjs';

import { posix } from 'path-browserify';

import { UserUploadTask } from '../models/user-upload-task';

@Injectable({
  providedIn: 'root',
})
export class UploadService {
  private uploadSubject = new BehaviorSubject<UserUploadTask[]>([]);
  private uploadProgress$: Observable<UserUploadTask[]> =
    this.uploadSubject.asObservable();

  constructor(
    private auth: Auth,
    private functions: Functions,
    private storage: Storage
  ) {}

  /**
   * Retrieves the current upload progress queue as an observable.
   *
   * @returns The upload progress queue observable.
   */
  getUploadProgress(): Observable<UserUploadTask[]> {
    return this.uploadProgress$;
  }

  /**
   * Removes an upload progress object from the queue.
   *
   * @param upload - A reference to the UserFileUpload object inside the array.
   */
  removeUpload(upload: UserUploadTask): void {
    const currentUploads = this.uploadSubject.value;
    const removed = currentUploads.filter((u) => u.id != upload.id);
    this.uploadSubject.next(removed);
  }

  /**
   * Takes a file with a specified user upload path and
   * adds corresponding UserUploadTask object to the upload progress queue.
   *
   * @param file - File object to be uploaded.
   * @param path - path of where to upload within the user's dashboard file system
   * @param onCompleted - callback function to run when file successfully uplaods
   * @param onError - callback function to run if file errors during upload
   */
  uploadFile(
    file: File,
    path: string,
    onCompleted?: (uploadRef: UserUploadTask) => void,
    onError?: (uploadRef: UserUploadTask) => void
  ): void {
    // add new file UploadTask to the observable
    const newUpload: UserUploadTask = {
      id: crypto.randomUUID(),
      file,
      name: file.name,
      type: 'file',
      progress: 0,
      uploadPath: path, // the upload path relative to user's directory
      downloadURL: null,
      status: 'pending' as const,
      onCompleted,
      onError,
    };

    this.uploadSubject.next([...this.uploadSubject.value, newUpload]);
    // queue the newly added file to upload
    this.uploadFileToStorage(newUpload);
  }

  /**
   * Takes a new folder name with a specified user upload path and
   * adds corresponding UserUploadTask object to the upload progress queue.
   *
   * @param name - name of the folder to be created.
   * @param path - path of where to upload within the user's dashboard file system
   * @param onCompleted - callback function to run when folder successfully uplaods
   * @param onError - callback function to run if folder errors during creation
   */
  createNewFolder(
    name: string,
    path: string,
    onCompleted?: (uploadRef: UserUploadTask) => void,
    onError?: (uploadRef: UserUploadTask) => void
  ): void {
    if (name == '') {
      return;
    }

    const folderUpload: UserUploadTask = {
      id: crypto.randomUUID(),
      file: null,
      name: name,
      type: 'folder',
      progress: 0,
      uploadPath: path, // the upload path relative to user's directory
      downloadURL: null,
      status: 'pending' as const,
      onCompleted,
      onError,
    };

    this.uploadSubject.next([...this.uploadSubject.value, folderUpload]);

    const cloudFunctionCallable = httpsCallable(
      this.functions,
      'request_user_create_folder'
    );

    cloudFunctionCallable({ name, path })
      .then(() => {
        this.updateUpload(folderUpload, { status: 'completed' });
      })
      .catch((error) => {
        this.updateUpload(folderUpload, { status: 'error' });
      });
  }

  /**
   * Takes a UserUploadTask of type='file' and uploads to Cloud Storage.
   * Calls updateUpload() during file upload progress accordingly.
   *
   * @param upload - the UserUploadTask to update.
   */
  private uploadFileToStorage(upload: UserUploadTask): void {
    if (!upload.file) {
      throw new Error(
        'UserUploadTask did not contain a file. You may have passed in a folder upload.'
      );
    }

    const cloudStoragePath = posix.join(
      'uploads',
      // using the Auth library directly isn't amazing since it can be null
      // for the first few secs when loading, but making a subscription inside
      // a service seems to be discouraged. for now, i will do this.
      this.auth.currentUser!.uid,
      upload.uploadPath,
      upload.file.name // add name check here
    );

    const storageRef = ref(this.storage, cloudStoragePath);
    const uploadTask = uploadBytesResumable(storageRef, upload.file);

    uploadTask.on(
      'state_changed',
      (snapshot) => {
        const progress =
          (snapshot.bytesTransferred / snapshot.totalBytes) * 100;
        this.updateUpload(upload, { progress, status: 'uploading' });
      },
      (error) => {
        this.updateUpload(upload, { status: 'error' });
      },
      async () => {
        try {
          const downloadURL = await getDownloadURL(uploadTask.snapshot.ref);
          this.updateUpload(upload, { downloadURL, status: 'completed' });
        } catch (error) {
          console.error('Failed to get download URL: ', error);
          this.updateUpload(upload, { status: 'error' });
        }
      }
    );
  }

  /**
   * Takes a UserUploadTask object with a set of updates and updates them in
   * the Observable accordingly. Also calls the upload's onCompleted and
   * onError handlers if applicable.
   *
   * @param upload - the UserUploadTask to update.
   * @param updates - object containing attribute updates. Do not update 'id' field. If 'status' is updated to 'completed' or 'error' the upload's onCompleted/onError handlers will fire.
   */
  private updateUpload(
    upload: UserUploadTask,
    updates: Omit<Partial<UserUploadTask>, 'id'>  // do not edit ID.
  ): void {
    const currentUploads = this.uploadSubject.value;
    let newUpload: UserUploadTask | undefined = undefined;

    const updated = currentUploads.map((u) =>
      u.id === upload.id
        ? (() => {
            newUpload = { ...upload, ...updates };
            return newUpload;
          })()
        : u
    );

    this.uploadSubject.next(updated);

    // once a particular UserFileUpload is updated, check if
    // we should run its onCompleted/onError and run it
    // (my way of making event handlers like Cloud Storage)
    if (newUpload && updates.status === 'completed') {
      upload.onCompleted && upload.onCompleted(newUpload);
    } else if (newUpload && updates.status === 'error') {
      upload.onError && upload.onError(newUpload);
    }
  }
}
