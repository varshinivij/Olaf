import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
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

  sendMessage(history: ChatMessage[]): Observable<any> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });

    return this.http.post(
      this.chatAPIEndpoint,
      { history },
      {
        headers,
        observe: 'events',
        reportProgress: true,
        responseType: 'text',
      }
    );
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
          reportProgress: true,
          responseType: 'json',
        }
      )
    );

    // Return the message from the response
    return response.message;
  }
}
