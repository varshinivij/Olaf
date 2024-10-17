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
import { Session } from '../models/session';
import { ChatMessage } from '../models/chat-message';
import { UserService } from './user.service';
import { HttpClient } from '@angular/common/http';
import { firstValueFrom } from 'rxjs';

// Make sure there are no duplicate imports
// import { UserService } from './user.service';  // Already imported
// import { ChatMessage } from '../models/chat-message'; // Already imported
import { Project } from '../models/project';
// import { Session } from '../models/session'; // Already imported

@Injectable({
  providedIn: 'root',
})
export class SessionsService {
  /**
   * Returns a list of all the user's sessions regardless of project.
   */
  userSessions: Session[] = [];
  activeSession: Session = this.createNewSession();
  private readonly API_URL = 'http://127.0.0.1:5001/twocube-web/us-central1'; // Update this with your GCP function URL

  constructor(
    private http: HttpClient,
    private userService: UserService,
    private firestore: Firestore, // Make sure Firestore is injected
  ) {
    this.loadAllSessions();
  }

  createNewSession(): Session {
    return {
      id: '',
      userId: '',
      projectId: '',
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

  changeUserSession(session: Session) {
    this.activeSession = session;
  }

  changeUserSessionById(id: string) {
    this.activeSession =
      this.userSessions.find((session) => session.id === id) ||
      this.createNewSession();
  }

  async saveActiveSession() {
    console.log('saving');
    await this.ensureUserIsSet(this.activeSession);
    console.log('user set');
    console.log(this.activeSession.id);
    if (!this.activeSession.id || this.activeSession.id === '') {
      const newSessionRef = await this.addNewSession(this.activeSession);
      if (newSessionRef) {
        this.activeSession.id = newSessionRef.id;
      } else {
        console.error('newSessionRef is undefined');
      }
      console.log(this.activeSession.id);
      this.userSessions.push(this.activeSession);
    }
    await this.saveSession(this.activeSession);
  }

  /**
   * Saves changes to the current active session. Writes a new session
   * to Firestore if the session is blank.
   */
  async saveSession(session: Session) {
    await this.ensureUserIsSet(session);
    if (session.id === '') {
      // This is a double write below - we should refactor this to only write once
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
      this.userSessions = (await this.queryUserSessions(currentUser.id)) || [];
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

  private addNewSession(session: Session) {
    console.log('adding');
    return this.http
      .post<{ id: string }>(`${this.API_URL}/add_session`, session)
      .toPromise();
  }

  private updateSession(session: Session) {
    return this.http
      .put(`${this.API_URL}/update_session`, session, {
        params: { id: session.id },
      })
      .toPromise();
  }

  /**
   * Deletes a specified session from Firestore and updates userSessions/activeSession.
   */
  async deleteSession(sessionId: string) {
    await this.http
      .delete(`${this.API_URL}/delete_session`, { params: { id: sessionId } })
      .toPromise();

    this.userSessions = this.userSessions.filter(
      (session) => session.id !== sessionId,
    );

    if (this.activeSession.id === sessionId) {
      this.newSession();
    }

    console.log(`Session ${sessionId} deleted successfully.`);
  }

  /**
   * Deletes all user sessions from Firestore and updates userSessions/activeSession.
   */
  async deleteAllSessions() {
    let promises = [];
    for (let session of this.userSessions) {
      promises.push(this.deleteSession(session.id));
    }
    await Promise.all(promises);
  }

  /**
   * Sets the userId of the session to the current user. Errors if user is
   * not logged in.
   */
  private async ensureUserIsSet(session: Session) {
    if (!session.userId) {
      session.userId = (await firstValueFrom(
        this.userService.getCurrentUser(),
      ))!.id;
    }
  }

  /*
    Firestore 'sessions' container private utility methods
  */

  private async queryUserSessions(userId: string) {
    const sessions = await this.http
      .get<Session[]>(`${this.API_URL}/get_sessions`, { params: { userId } })
      .toPromise();
    return sessions || [];
  }

  private createSession(session: Session) {
    return addDoc(collection(this.firestore, 'sessions'), session);
  }

  private newSession() {
    this.activeSession = this.createNewSession();
  }

  async addMessageToActiveSession(message: ChatMessage) {
    this.activeSession.history.push(message);
    await this.saveActiveSession();
  }
}
