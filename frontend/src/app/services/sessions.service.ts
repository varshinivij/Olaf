import { Injectable } from '@angular/core';
import { Session, createNewSession } from '../models/session';
import { ChatMessage } from '../models/chat-message';
import { UserService } from './user.service';
import {
  Firestore,
  addDoc,
  getDocs,
  setDoc,
  collection,
  doc,
  query,
  where
} from '@angular/fire/firestore';
import { firstValueFrom } from 'rxjs';

@Injectable({
  providedIn: 'root'
})
export class SessionsService {

  constructor(private firestore: Firestore, private userService: UserService) {
    this.loadUserSessions();
  }

  userSessions: Session[] = [];
  activeSession: Session = createNewSession();

  changeUserSession(session: Session) {
    this.activeSession = session;
  }

  changeUserSessionById(id: string) {
    this.activeSession = this.userSessions.find((session) => session.id === id) || createNewSession();
  }

  async saveActiveSession() {
    console.log("saving");
    await this.ensureUserIsSet();
    console.log("user set");
    if (!this.activeSession.id) {
      console.log(this.activeSession);
      //this is a double write below - we should refactor this to only write once
      const newSessionRef = await this.addNewSession(this.activeSession);
      this.activeSession.id = newSessionRef.id;
      this.userSessions.push(this.activeSession);
    } else {
      await this.updateSession(this.activeSession);
    }
  }

  private async ensureUserIsSet() {
    if (!this.activeSession.userId) {
      const user = await firstValueFrom(this.userService.getCurrentUser());
      if (user) {
        this.activeSession.userId = user.id;
      }
    }
  }

  private addNewSession(session: Session) {
    console.log("adding");
    return addDoc(collection(this.firestore, 'sessions'), session);
  }

  private updateSession(session: Session) {
    return setDoc(doc(this.firestore, 'sessions', session.id), session);
  }

  private async queryUserSessions(userId: string) {
    const sessionsQuery = query(
      collection(this.firestore, 'sessions'),
      where('userId', '==', userId)
    );
    const sessionsSnapshot = await getDocs(sessionsQuery);
    return sessionsSnapshot.docs.map(doc => doc.data() as Session);
  }

  /** Load all sessions for the current user */
  async loadUserSessions() {
    console.log("loading");
    const currentUser = await firstValueFrom(this.userService.getCurrentUser())
    if (currentUser) {
      console.log("user found");
      this.userSessions = await this.queryUserSessions(currentUser.id);
      console.log(this.userSessions);
    }
  }

  newSession() {
    this.activeSession = createNewSession();
  }

  addMessageToActiveSession(message: ChatMessage) {
    this.activeSession.history.push(message);
  }
}