/*
  TODO:
    1) Any user services that will be necessary in the future (edit password, etc.)
    2) (if necessary) potentially also expose a signal that holds currentUser as well
    3) Consider storing state in local storage to persist data after a refresh.
       https://medium.com/@aayyash/authentication-in-angular-why-it-is-so-hard-to-wrap-your-head-around-it-23ea38a366de
    4) Consider working on uploadProfilePicture to shrink image size before upload
    5) Promise resolveAll instead of consecutive awaits when possible
    6) Using auth.currentUser can be inconsistent since it can take a little while to
       update after component init/refreshing. Maybe just stick to observables.
*/

import { Injectable, Injector } from '@angular/core';
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
  deleteDoc,
  doc,
  docData,
  getDoc,
  setDoc,
} from '@angular/fire/firestore';
import {
  Storage,
  deleteObject,
  getDownloadURL,
  ref,
  uploadBytesResumable,
} from '@angular/fire/storage';
import { Observable, of, map, switchMap, shareReplay } from 'rxjs';

import { FileStorageService } from './file-storage.service';
import { SessionsService } from './sessions.service';

import { User } from '../models/user';

@Injectable({
  providedIn: 'root',
})
export class UserService {
  private user$: Observable<User | null>;

  constructor(
    private auth: Auth,
    private firestore: Firestore,
    private storage: Storage,
    private injector: Injector
  ) {
    this.user$ = user(this.auth).pipe(
      switchMap((user: FirebaseUser | null) => {
        if (user) {
          const userDocRef = doc(this.firestore, 'users', user.uid);
          return docData(userDocRef).pipe(
            map((user: any) => {
              return {
                ...user,
                createdAt: user.createdAt?.toDate(),
                updatedAt: user.updatedAt?.toDate(),
              } as User;
            }),
            shareReplay(1)
          ) as Observable<User>;
        } else {
          return of(null);
        }
      })
    );
  }

  /**
   * Retrieves the current user as an observable.
   *
   * This observable updates on authentication state changes as well as Firestore user document changes. It (should) always emits an initial value upon subscription.
   *
   * @returns An Observable of the current user or null if not authenticated.
   */
  getCurrentUser(): Observable<User | null> {
    return this.user$;
  }

  /**
   * Returns a promise to retrieve an ID token for the user. Used for
   * Firebase user authentication in pure HTTP endpoints.
   *
   * @returns A Promise returning an ID token for the user.
   */
  getIdToken(): Promise<string> {
    return this.auth.currentUser!.getIdToken();
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
   * Updates the account information for the current user.
   *
   * @param data A subset of a User object containing user properties to
   * update in Firestore. Updating user.id will have no effect.
   */
  async updateAccount(data: Partial<User>): Promise<void> {
    await this.updateUserInfo(this.auth.currentUser!, data);
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

  /**
   * Uploads a file to Cloud Storage under profilePictures/{userID}.
   * Updates user profile picture in Firestore record on completion.
   * Promise returns the downloadURL to the uploaded picture.
   *
   * TODO: somehow minify the image size so we don't end up with GBs of pics
   *
   * @param file The file to be uploaded.
   * @param onProgress Callback to be run on progress snapshot.
   */
  async uploadProfilePicture(
    file: File,
    onProgress: (progress: number) => void
  ): Promise<string> {
    const storageRef = ref(
      this.storage,
      `profilePictures/${this.auth.currentUser!.uid}`
    );

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
          this.updateAccount({
            profilePictureUrl: downloadUrl,
          });
          resolve(downloadUrl);
        }
      );
    });
  }

  /**
   * Converts an error thrown by Firebase Auth into a user-friendly
   * error message. If the error was not thrown by Firebase Auth
   * (it should have the 'code' property) it returns a default
   * error message. A list of codes can be found here:
   * https://firebase.google.com/docs/auth/admin/errors
   *
   * @param error The error thrown by Firebase Auth.
   * @returns A user-friendly error message.
   */
  static convertAuthErrorToMessage(error: any) {
    switch (error.code) {
      case 'auth/user-disabled': {
        return 'Account has been disabled';
      }
      case 'auth/invalid-credential': {
        return 'Incorrect email or password';
      }
      case 'auth/email-already-in-use': {
        return 'Account already exists';
      }
      default: {
        return 'Something went wrong, try again later';
      }
    }
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
      name: user.displayName,
      role: 'someRole',
      status: 'someStatus',
      createdAt: new Date(),
      updatedAt: new Date(),
      organization: null,
      profilePictureUrl: user.photoURL,
    };

    const userDocument = doc(this.firestore, 'users', user.uid);
    await setDoc(userDocument, data);
  }

  // updates Firestore user data based on given FirebaseUser and update data
  private async updateUserInfo(
    user: FirebaseUser,
    data: Partial<User>
  ): Promise<void> {
    const userDocument = doc(this.firestore, 'users', user.uid);

    // should not be able to update user id
    delete data.id;
    // Update Firebase Auth profile fields
    // TODO: Updating password should be done through a separate method for security reasons
    if (data.email) await updateEmail(user, data.email);
    if (data.name) await updateProfile(user, { displayName: data.name });

    await setDoc(
      userDocument,
      { ...data, updatedAt: new Date() },
      { merge: true }
    );
  }

  private async deleteProfilePicture(user: FirebaseUser) {
    const storageRef = ref(this.storage, `profilePictures/${user.uid}`);
    await deleteObject(storageRef);
  }

  // deletes Firestore and Firebase Auth user data based on given FirebaseUser
  private async deleteUserInfo(user: FirebaseUser): Promise<void> {
    // delete Firebase auth record
    await user.delete();
    // delete Firestore record
    await deleteDoc(doc(this.firestore, 'users', user.uid));
    // delete profile picture
    await this.deleteProfilePicture(user);
    // use the FileStorageService to delete uploaded files
    const fileStorageService = this.injector.get(FileStorageService);
    await fileStorageService.deletePath('/');
    // use the SessionService to delete all sessions
    const sessionService = this.injector.get(SessionsService);
    await sessionService.deleteAllSessions();
  }
}
