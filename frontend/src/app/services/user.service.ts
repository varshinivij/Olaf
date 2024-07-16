import { Injectable, WritableSignal, inject, signal } from '@angular/core';
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
  Timestamp,
  deleteDoc,
  doc,
  getDoc,
  setDoc,
} from '@angular/fire/firestore';
import { Observable } from 'rxjs';
import { User as UserProfile } from '../models/user';
import {
  FirebaseStorage,
  getDownloadURL,
  getStorage,
  ref,
  uploadBytes,
  uploadBytesResumable,
} from '@angular/fire/storage';

@Injectable({
  providedIn: 'root',
})
export class UserService {
  private auth: Auth = inject(Auth);
  private firestore: Firestore = inject(Firestore);
  private storage: FirebaseStorage = getStorage();

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
    const credential = await signInWithEmailAndPassword(
      this.auth,
      email,
      password
    );
    // just in case an account exists in Auth but somehow not in the database
    if (!(await this.checkUserExists(credential.user))) {
      await this.createUserInfo(credential.user);
    }
    this.currentUser.set(await this.readUserInfo(credential.user));
  }

  async signupWithEmail(email: string, password: string) {
    const credential = await createUserWithEmailAndPassword(
      this.auth,
      email,
      password
    );
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
      organization: 'someOrganization',
      profilePictureUrl: null,
    }; // userID is stored as documentID so no need to put in document

    const userDocument = doc(this.firestore, 'users', user.uid);
    await setDoc(userDocument, data);
  }

  // reads Firestore user data based on given FirebaseUser as UserProfile
  public async readUserInfo(user: FirebaseUser): Promise<UserProfile> {
    const docSnap = await getDoc(doc(this.firestore, 'users', user.uid));
    const data = docSnap.data()!;

    return {
      id: docSnap.id,
      email: data['email'],
      name: data['name'],
      role: data['role'],
      status: data['status'],
      createdAt: (data['createdAt'] as Timestamp).toDate(),
      updatedAt: (data['updatedAt'] as Timestamp).toDate(),
      organization: data['organization'],
      profilePictureUrl: data['profilePictureUrl'],
    } as UserProfile;
  }

  // updates Firestore user data based on given FirebaseUser and update data
  public async updateUserInfo(
    user: FirebaseUser,
    data: Partial<Omit<UserProfile, 'id'>>
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
    userId: string,
    file: File,
    onProgress: (progress: number) => void
  ): Promise<string> {
    const storageRef = ref(this.storage, `user/profile/${userId}/${file.name}`);
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
