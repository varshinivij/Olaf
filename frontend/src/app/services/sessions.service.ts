import { Injectable } from '@angular/core';
import {
  Firestore,
  addDoc,
  collection,
  deleteDoc,
  doc,
  getDocs,
  query,
  setDoc,
  where,
} from '@angular/fire/firestore';
import { firstValueFrom } from 'rxjs';

import { UserService } from './user.service';

import { ChatMessage } from '../models/chat-message';
import { Project } from '../models/project';
import { Session } from '../models/session';

@Injectable({
  providedIn: 'root',
})
export class SessionsService {
  /**
   * Returns a list of all the user's sessions regardless of project.
   */
  userSessions: Session[] = [];

  constructor(
    private firestore: Firestore,
    private userService: UserService,
  ) {
    this.loadAllSessions();
  }

  blankSession(project: Project): Session {
    return {
      id: '',
      userId: '',
      projectId: project.id,
      name: '<untitled session>',
      context: '',
      history: [
        {
          type: 'text',
          role: 'assistant',
          content: 'Hello, how can I help you today?',
        },
      ],
      sandboxId: null,
    };
  }

  /**
   * Renames an existing session and writes to Firestore.
   */
  async renameSession(session: Session, name: string) {
    if (!name.trim()) {
      return;
    }
    session.name = name;
    await this.saveSession(session);
  }

  /**
   * Saves changes to the current active session. Writes a new session
   * to Firestore if the session is blank.
   */
  async saveSession(session: Session) {
    await this.ensureUserIsSet(session);
    if (session.id === '') {
      //this is a double write below - we should refactor this to only write once
      const newSessionRef = await this.createSession(session);
      session.id = newSessionRef.id;
      await this.userSessions.push(session);
    }
    await this.updateSession(session);
  }

  /**
   * Loads all sessions for the current user into userSessions.
   */
  async loadAllSessions() {
    const currentUser = await firstValueFrom(this.userService.getCurrentUser());
    if (currentUser) {
      this.userSessions = await this.queryUserSessions(currentUser.id);
      console.log(this.userSessions);
    }
  }

  /**
   * Pushes a message to the current session and writes to Firestore.
   */
  async addMessageToSession(session: Session, message: ChatMessage) {
    session.history.push(message);
    await this.saveSession(session);
  }

  /**
   * Sets the sandbox id of the current session and writes to Firestore.
   */
  async setSessionSandBoxId(session: Session, sandboxId: string) {
    session.sandboxId = sandboxId;
    await this.saveSession(session);
  }

  /**
   * Deletes a specified session from Firestore and updates userSessions/activeSession.
   */
  async deleteSession(session: Session) {
    await deleteDoc(doc(this.firestore, 'sessions', session.id));

    this.userSessions = this.userSessions.filter(
      (s) => s.id !== session.id,
    );

    console.log(`Session ${session.id} deleted successfully.`);
  }

  /**
   * Deletes all user sessions from Firestore and updates userSessions/activeSession.
   */
  async deleteAllSessions() {
    let promises = [];
    for (let session of this.userSessions) {
      promises.push(this.deleteSession(session));
    }
    await Promise.all(promises);
  }

  /**
   * Sets the userId of the session to the current user. Errors if user is
   * not logged in.
   */
  private async ensureUserIsSet(session: Session) {
    if (!session.userId) {
      session.userId = (await firstValueFrom(this.userService.getCurrentUser()))!.id;
    }
  }

  /*
    Firestore 'sessions' container private utility methods
  */

  private async queryUserSessions(userId: string) {
    const sessionsQuery = query(
      collection(this.firestore, 'sessions'),
      where('userId', '==', userId),
    );
    const sessionsSnapshot = await getDocs(sessionsQuery);
    return sessionsSnapshot.docs.map((doc) => doc.data() as Session);
  }

  private createSession(session: Session) {
    return addDoc(collection(this.firestore, 'sessions'), session);
  }

  private updateSession(session: Session) {
    return setDoc(doc(this.firestore, 'sessions', session.id), session);
  }
}
