import { Injectable } from '@angular/core';
import { Auth } from '@angular/fire/auth';
import {
  Firestore,
  deleteDoc,
  doc,
  docData,
  getDoc,
  setDoc,
} from '@angular/fire/firestore';

@Injectable({
  providedIn: 'root',
})
export class UploadService {
  // currently uses angularfire auth.currentUser to read user state. may
  // need to relook later, not sure if currentUser is always synced.
  // i also heard people sometimes save current user in localStorage.
  constructor(private auth: Auth) {}

  async uploadFile(): Promise<void> {
    console.log(this.auth.currentUser);
  }
}
