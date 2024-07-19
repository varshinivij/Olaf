/*
  TODO:
    1) Write updateUserInfo() function
    2) Any user services that will be necessary in the future (edit account, etc.)
    3) (if necessary) potentially also expose a signal that holds currentUser as well
    4) Consider placing a .tap(?) on Observable and storing state in local storage?
       (maybe replace with readUserInfo? not sure.)
       https://medium.com/@aayyash/authentication-in-angular-why-it-is-so-hard-to-wrap-your-head-around-it-23ea38a366de
*/

import { Injectable } from '@angular/core';
import {
  Auth,
  GoogleAuthProvider,
  User as FirebaseUser,
  createUserWithEmailAndPassword,
  signInWithEmailAndPassword,
  signInWithPopup,
  updateEmail,
  updatePassword,
  updateProfile,
  user,
} from '@angular/fire/auth';
import {
  Firestore,
  Timestamp,
  deleteDoc,
  doc,
  docData,
  getDoc,
  setDoc,
} from '@angular/fire/firestore';
import { Observable, of, switchMap } from 'rxjs';

import { User } from '../models/user';
import {
  FirebaseStorage,
  getDownloadURL,
  getStorage,
  ref,
  uploadBytesResumable,
} from '@angular/fire/storage';

@Injectable({
  providedIn: 'root',
})
export class UserService {
  private user$: Observable<User | null>;
  private storage: FirebaseStorage = getStorage();

  constructor(private auth: Auth, private firestore: Firestore) {
    this.user$ = user(this.auth).pipe(
      switchMap((user: FirebaseUser | null) => {
        if (user) {
          const userDocRef = doc(this.firestore, 'users', user.uid);
          return docData(userDocRef) as Observable<User>;
        } else {
          return of(null);
        }
      })
    );
  }

  /**
   * Retrieves the current user as an observable.
   *
   * This observable updates on authentication state changes as well as Firestore user document changes.
   *
   * @returns An Observable of the current user or null if not authenticated.
   */
  getCurrentUser(): Observable<User | null> {
    return this.user$;
  }

  /**
   * Authenticates a user using Google Sign-In.
   *
   * If the user does not exist in Firestore, their information is automatically created.
   */
  async loginWithGoogle(): Promise<void> {
    const provider = new GoogleAuthProvider();
    const credential = await signInWithPopup(this.auth, provider);
    // login with google w/o account will register in Auth anyway (special case)
    if (!(await this.checkUserExists(credential.user))) {
      await this.createUserInfo(credential.user);
    }
  }

  /**
   * Authenticates a user with email and password.
   *
   * @param email The user's email address.
   * @param password The user's password.
   */
  async loginWithEmail(email: string, password: string): Promise<void> {
    const credential = await signInWithEmailAndPassword(
      this.auth,
      email,
      password
    );
    // just in case an account exists in Auth but somehow not in the database
    if (!(await this.checkUserExists(credential.user))) {
      await this.createUserInfo(credential.user);
    }
  }

  /**
   * Signs up a new user with an email and password.
   *
   * After successful authentication, it creates user information in Firestore.
   *
   * @param email The email address of the new user.
   * @param password The password for the new user.
   */
  async signupWithEmail(email: string, password: string): Promise<void> {
    const credential = await createUserWithEmailAndPassword(
      this.auth,
      email,
      password
    );
    await this.createUserInfo(credential.user);
  }

  /**
   * Logs out the current user.
   *
   * This will end the user's session and clear any authentication tokens.
   */
  async logout(): Promise<void> {
    await this.auth.signOut();
  }

  /**
   * Deletes the current user's account.
   *
   * This will remove the user's Firebase Auth record and delete their information from Firestore.
   * Note: This action is irreversible.
   */
  async deleteAccount(): Promise<void> {
    await this.deleteUserInfo(this.auth.currentUser!);
  }

  /*
    Firestore 'users' container private utility methods
  */

  // checks if a given FirebaseUser exists in the Firestore users container
  private async checkUserExists(user: FirebaseUser): Promise<boolean> {
    const userDocRef = doc(this.firestore, 'users', user.uid);
    const docSnap = await getDoc(userDocRef);
    return docSnap.exists();
  }

  // creates new Firestore user data based on given FirebaseUser
  private async createUserInfo(user: FirebaseUser): Promise<void> {
    const data = {
      id: user.uid,
      email: user.email || 'noemail@example.com',
      name: user.displayName || 'No Name',
      role: 'someRole',
      status: 'someStatus',
      createdAt: new Date(),
      updatedAt: new Date(),
      organization: 'someOrganization',
      profilePictureUrl: null,
    };

    const userDocument = doc(this.firestore, 'users', user.uid);
    await setDoc(userDocument, data);
  }

  // reads Firestore user data based on given FirebaseUser as UserProfile
  private async readUserInfo(user: FirebaseUser): Promise<User> {
    const docSnap = await getDoc(doc(this.firestore, 'users', user.uid));
    const data = docSnap.data()!;
    return data as User;
  }

  // updates Firestore user data based on given FirebaseUser and update data
  public async updateUserInfo(
    user: FirebaseUser,
    data: Partial<User>
  ): Promise<void> {
    const userDocument = doc(this.firestore, 'users', user.uid);
    // Update Firebase Auth profile fields
    if (data.email) await updateEmail(user, data.email);
    if (data.name) await updateProfile(user, { displayName: data.name });
    // Note: Updating password should be done through a separate method for security reasons

    await setDoc(
      userDocument,
      { ...data, updatedAt: new Date() },
      { merge: true }
    );
  }

  // deletes Firestore and Firebase Auth user data based on given FirebaseUser
  private async deleteUserInfo(user: FirebaseUser): Promise<void> {
    await user.delete();
    await deleteDoc(doc(this.firestore, 'users', user.uid));
  }

  async uploadProfilePicture(
    file: File,
    onProgress: (progress: number) => void
  ): Promise<string> {
    const storageRef = ref(this.storage, `user/profile/${file.name}`);
    const uploadTask = uploadBytesResumable(storageRef, file);

    return new Promise((resolve, reject) => {
      uploadTask.on(
        'state_changed',
        (snapshot) => {
          const progress =
            (snapshot.bytesTransferred / snapshot.totalBytes) * 100;
          onProgress(progress);
        },
        (error) => {
          reject(error);
        },
        async () => {
          const downloadUrl = await getDownloadURL(uploadTask.snapshot.ref);
          resolve(downloadUrl);
        }
      );
    });
  }
}
