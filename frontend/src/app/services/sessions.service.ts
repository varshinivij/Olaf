import { Injectable } from '@angular/core';
import { HttpClient } from '@angular/common/http';
import { firstValueFrom } from 'rxjs';

import { UserService } from './user.service';

import { ChatMessage } from '../models/chat-message';
import { Project } from '../models/project';
import { Session } from '../models/session';

@Injectable({
  providedIn: 'root',
})
export class SessionsService {
  userSessions: Session[] = [];
  private baseUrl = 'http://127.0.0.1:5001/twocube-web/us-central1/';

  constructor(
    private http: HttpClient,
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

  async renameSession(session: Session, name: string) {
    if (!name.trim()) {
      return;
    }
    session.name = name;
    await this.saveSession(session);
  }

  async saveSession(session: Session) {
    await this.ensureUserIsSet(session);
    if (session.id === '') {
      const response = await this.http
        .post<any>(`${this.baseUrl}create_session`, session)
        .toPromise();
      session.id = response.id;
      this.userSessions.push(session);
    } else {
      await this.http
        .put(`${this.baseUrl}update_session?id=${session.id}`, session)
        .toPromise();
    }
  }

  async loadAllSessions() {
    const currentUser = await firstValueFrom(this.userService.getCurrentUser());
    if (currentUser) {
      this.userSessions =
        (await this.http
          .get<
            Session[]
          >(`${this.baseUrl}get_sessions?userId=${currentUser.id}`)
          .toPromise()) ?? [];

      console.log(this.userSessions);
    }
  }

  async addMessageToSession(session: Session, message: ChatMessage) {
    session.history.push(message);
    // await this.saveSession(session);
  }

  async setSessionSandBoxId(session: Session, sandboxId: string) {
    session.sandboxId = sandboxId;
    await this.saveSession(session);
  }

  async deleteSession(session: Session) {
    await this.http
      .delete(`${this.baseUrl}delete_session?id=${session.id}`)
      .toPromise();
    this.userSessions = this.userSessions.filter((s) => s.id !== session.id);
    console.log(`Session ${session.id} deleted successfully.`);
  }

  async deleteAllSessions() {
    const currentUser = await firstValueFrom(this.userService.getCurrentUser());
    if (currentUser) {
      await this.http
        .delete(`${this.baseUrl}delete_all_sessions?userId=${currentUser.id}`)
        .toPromise();
      this.userSessions = [];
    }
  }

  private async ensureUserIsSet(session: Session) {
    if (!session.userId) {
      session.userId = (await firstValueFrom(
        this.userService.getCurrentUser(),
      ))!.id;
    }
  }
}
