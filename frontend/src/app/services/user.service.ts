/*
  TODO:
    1) Write updateUserInfo() function (if necessary)?
    2) Any user services that will be necessary in the future (edit/delete account, etc.)?
    3) Route guard so you can't go to Home Page without logging in (+future pages)?
*/

import {
  Injectable,
  WritableSignal,
  inject,
  signal,
} from '@angular/core';
import {
  Auth,
  GoogleAuthProvider,
  User as FirebaseUser,
  createUserWithEmailAndPassword,
  signInWithPopup,
  signInWithEmailAndPassword,
  updateEmail,
  updatePassword,
  updateProfile,
  user,
} from '@angular/fire/auth';
import {
  Firestore,
  deleteDoc,
  doc,
  getDoc,
  setDoc,
} from '@angular/fire/firestore';
import { Observable, Subscription } from 'rxjs';
import { User as UserProfile } from '../models/user';

@Injectable({
  providedIn: 'root',
})
export class UserService {
  private auth: Auth = inject(Auth);
  private firestore: Firestore = inject(Firestore);

  // observable which streams events on login, logout, and token refresh.
  // emits raw User object from @angular/fire. likely unused outside of class.
  user$: Observable<FirebaseUser | null> = user(this.auth);

  // signal emitting the current user as a models/user.ts/User object.
  // access through UserService().currentUser().
  // from my understanding: use this signal within a template and it
  // auto-reacts to changes. use effect() to react within component logic.
  currentUser: WritableSignal<UserProfile | null> = signal<UserProfile | null>(
    null
  );

  constructor() {}

  async loginWithGoogle(): Promise<void> {
    const provider = new GoogleAuthProvider();
    const credential = await signInWithPopup(this.auth, provider);
    // login with google without an account will register in Auth anyway (special case)
    if (!(await this.checkUserExists(credential.user))) {
      await this.createUserInfo(credential.user);
    }
    this.currentUser.set(await this.readUserInfo(credential.user));
  }

  async loginWithEmail(email: string, password: string) {
    const credential = await signInWithEmailAndPassword(this.auth, email, password);
    // just in case an account exists in Auth but somehow not in the database
    if (!(await this.checkUserExists(credential.user))) {
      await this.createUserInfo(credential.user);
    }
    this.currentUser.set(await this.readUserInfo(credential.user));
  }

  async signupWithEmail(email: string, password: string) {
    const credential = await createUserWithEmailAndPassword(this.auth, email, password);
    await this.createUserInfo(credential.user);
    this.currentUser.set(await this.readUserInfo(credential.user));
  }

  async logout(): Promise<void> {
    await this.auth.signOut();
    this.currentUser.set(null);
  }

  async deleteAccount(): Promise<void> {
    await this.deleteUserInfo(this.auth.currentUser!);
    this.currentUser.set(null);
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
      email: user.email || 'noemail@example.com',
      name: user.displayName || 'No Name',
      role: 'someRole',
      status: 'someStatus',
      createdAt: new Date(),
      updatedAt: new Date(),
    }; // userID is stored as documentID so no need to put in document

    const userDocument = doc(this.firestore, 'users', user.uid);
    await setDoc(userDocument, data);
  }

  // reads Firestore user data based on given FirebaseUser as UserProfile
  private async readUserInfo(user: FirebaseUser): Promise<UserProfile> {
    const docSnap = await getDoc(doc(this.firestore, 'users', user.uid));
    const data = docSnap.data()!;

    return {
      id: docSnap.id,
      email: data['email'],
      name: data['name'],
      role: data['role'],
      status: data['status'],
      createdAt: data['createdAt'],
      updatedAt: data['updatedAt'],
    } as UserProfile;
  }

  // TODO: updates Firestore user data based on given FirebaseUser and update data
  private async updateUserInfo(
    user: FirebaseUser,
    data: Partial<Omit<UserProfile, 'id'>>
  ): Promise<void> {
    // should probably update Firestore data for user as well as relevant Auth
    // fields (updateEmail(), updatePassword(), updateProfile())?
  }

  // deletes Firestore and Firebase Auth user data based on given FirebaseUser
  private async deleteUserInfo(user: FirebaseUser): Promise<void> {
    await user.delete();
    await deleteDoc(doc(this.firestore, 'users', user.uid));
  }
}
