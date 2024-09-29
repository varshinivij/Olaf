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
  where,
  deleteDoc
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
    console.log(this.activeSession.id)
    if (!this.activeSession.id || this.activeSession.id === '') {
      //this is a double write below - we should refactor this to only write once
      const newSessionRef = await this.addNewSession(this.activeSession);
      this.activeSession.id = newSessionRef.id;
      console.log(this.activeSession.id);
      await this.userSessions.push(this.activeSession);
    }
    await this.updateSession(this.activeSession);
  }

  renameSession(session: Session, name: string) {
    session.name = name;
    this.updateSession(session);
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

  async addMessageToActiveSession(message: ChatMessage) {
    this.activeSession.history.push(message);
    await this.saveActiveSession();
  }

  setSessionSandBoxId(sandboxId: string) {
    this.activeSession.sandboxId = sandboxId;
    this.saveActiveSession();
  }

  async deleteSession(sessionId: string) {
    await deleteDoc(doc(this.firestore, 'sessions', sessionId));

    this.userSessions = this.userSessions.filter(session => session.id !== sessionId);

    if (this.activeSession.id === sessionId) {
      this.newSession();
    }

    console.log(`Session ${sessionId} deleted successfully.`);
  }
}
