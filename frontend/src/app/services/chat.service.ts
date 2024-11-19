import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import { Observable, firstValueFrom } from 'rxjs';
import { ChatMessage } from '../models/chat-message';

@Injectable({
  providedIn: 'root',
})
export class ChatService {
  private chatAPIEndpoint =
    'REMOVED'; // generalist chatder agent
  private nameMakerAPIEndpoint = 'REMOVED';

  constructor(private http: HttpClient) {}

  sendMessage(
    message: string,
    sessionId: string,
    userId: string,
    projectId: string
  ): Observable<any> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });

    // Set up query parameters
    const params = new HttpParams()
      .set('session_id', sessionId)
      .set('user_id', userId)
      .set('project_id', projectId);

    const body = {
      message,
    };

    return this.http.post(this.chatAPIEndpoint, body, {
      headers,
      params,
      observe: 'events',
      reportProgress: true,
      responseType: 'text',
    });
  }

  async generateChatNameFromHistory(history: ChatMessage[]): Promise<string> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });

    const response = await firstValueFrom(
      this.http.post<{ message: string }>(
        this.nameMakerAPIEndpoint,
        { history },
        {
          headers,
          responseType: 'json',
        }
      )
    );

    // Return the message from the response
    return response.message;
  }
}
