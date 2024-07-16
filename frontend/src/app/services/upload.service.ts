import { Injectable } from '@angular/core';
import { Auth } from '@angular/fire/auth';
import {
  Storage,
  getDownloadURL,
  ref,
  uploadBytesResumable,
} from '@angular/fire/storage';
import { BehaviorSubject, Observable } from 'rxjs';

import { UploadProgress } from '../models/upload-progress';

@Injectable({
  providedIn: 'root',
})
export class UploadService {
  private uploadSubject = new BehaviorSubject<UploadProgress[]>([]);
  private uploadProgress$: Observable<UploadProgress[]> =
    this.uploadSubject.asObservable();

  constructor(private auth: Auth, private storage: Storage) {}

  /**
   * Retrieves the current upload queue progress as an Observable.
   *
   * @returns An Observable of UploadProgress[].
   */
  getUploadProgress(): Observable<UploadProgress[]> {
    return this.uploadProgress$;
  }

  /**
   * Adds new files to the upload queue without removing existing ones.
   *
   * @param files An array of File objects to be added to the upload queue.
   */
  addFiles(files: File[]): void {
    const currentUploads = this.uploadSubject.value;
    const newUploads = files.map((file) => ({
      file,
      progress: 0,
      downloadURL: null,
      status: 'pending' as const,
    }));
    this.uploadSubject.next([...currentUploads, ...newUploads]);
  }

  /**
   * Replaces the current upload queue with a new set of files.
   *
   * @param files An array of File objects to set as the new upload queue.
   */
  setFiles(files: File[]): void {
    const newUploads = files.map((file) => ({
      file,
      progress: 0,
      downloadURL: null,
      status: 'pending' as const,
    }));
    this.uploadSubject.next(newUploads);
  }

  /**
   * Removes a file from the upload queue at the specified index.
   *
   * @param index The index of the file to be removed from the upload queue.
   */
  removeFile(index: number): void {
    const currentUploads = this.uploadSubject.value;
    currentUploads.splice(index, 1);
    this.uploadSubject.next(currentUploads);
  }

  /**
   * Initiates the upload process for all pending files in the queue.
   * Uploads to currently logged in user's Firebase Storage folder at
   * uploads/{userID}/{fileName}.
   */
  uploadFiles(): void {
    const uploads = this.uploadSubject.value;
    uploads.forEach((upload, index) => {
      if (upload.status === 'pending') {
        this.uploadFile(upload.file, index);
      }
    });
  }

  private updateUpload(index: number, updates: Partial<UploadProgress>): void {
    const currentUploads = this.uploadSubject.value;
    currentUploads[index] = { ...currentUploads[index], ...updates };
    this.uploadSubject.next(currentUploads);
  }

  private uploadFile(file: File, index: number): void {
    const filePath = `uploads/${this.auth.currentUser!.uid}/${file.name}`;
    const storageRef = ref(this.storage, filePath);
    const uploadTask = uploadBytesResumable(storageRef, file);

    uploadTask.on(
      'state_changed',
      (snapshot) => {
        const progress =
          (snapshot.bytesTransferred / snapshot.totalBytes) * 100;
        this.updateUpload(index, { progress, status: 'uploading' });
      },
      (error) => {
        console.error('Upload failed: ', error);
        this.updateUpload(index, { status: 'error' });
      },
      async () => {
        try {
          const downloadURL = await getDownloadURL(uploadTask.snapshot.ref);
          this.updateUpload(index, { downloadURL, status: 'completed' });
        } catch (error) {
          console.error('Failed to get download URL: ', error);
          this.updateUpload(index, { status: 'error' });
        }
      }
    );
  }
}
