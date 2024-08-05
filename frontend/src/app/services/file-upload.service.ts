import { Injectable } from '@angular/core';
import { Auth } from '@angular/fire/auth';
import {
  Storage,
  getDownloadURL,
  ref,
  uploadBytesResumable,
} from '@angular/fire/storage';
import { BehaviorSubject, Observable } from 'rxjs';

import { posix } from 'path-browserify';

import { UserFileUpload } from '../models/user-file-upload';

@Injectable({
  providedIn: 'root',
})
export class FileUploadService {
  private uploadSubject = new BehaviorSubject<UserFileUpload[]>([]);
  private uploadProgress$: Observable<UserFileUpload[]> =
    this.uploadSubject.asObservable();

  constructor(private auth: Auth, private storage: Storage) {}

  /**
   * Retrieves the current upload progress queue as an observable.
   */
  getUploadProgress(): Observable<UserFileUpload[]> {
    return this.uploadProgress$;
  }

  /**
   * Removes an upload progress object from the queue.
   *
   * @param upload - A reference to the UserFileUpload object inside the array.
   */
  removeUpload(upload: UserFileUpload): void {
    const currentUploads = this.uploadSubject.value;
    const removed = currentUploads.filter((u) => u != upload);
    this.uploadSubject.next(removed);
  }

  /**
   * Takes an array of files with a specified user upload path and
   * adds corresponding UserFileUpload objects to the upload progress queue.
   *
   * @param files - Array of file objects to be uploaded.
   * @param path - path of where to upload within the user's dashboard file system
   */
  uploadFiles(files: File[], path: string): void {
    let currentUploads = this.uploadSubject.value;

    // add new files to the observable
    const newUploads = files.map(
      (file) =>
        ({
          file,
          progress: 0,
          uploadPath: path, // the upload path relative to user's directory
          downloadURL: null,
          status: 'pending' as const,
        } as UserFileUpload)
    );

    this.uploadSubject.next([...currentUploads, ...newUploads]);

    // queue the newly added files to upload
    newUploads.forEach((upload) => {
      this.uploadFile(upload);
    });
  }

  private uploadFile(upload: UserFileUpload): void {
    // TODO: check for dupe file names. also add upload folders.

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
        console.error('Upload failed: ', error);
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

  private updateUpload(
    upload: UserFileUpload,
    updates: Partial<UserFileUpload>
  ): void {
    const currentUploads = this.uploadSubject.value;
    const updated = currentUploads.map((u) =>
      u == upload ? { ...upload, ...updates } : upload
    );

    this.uploadSubject.next(updated);
  }
}
