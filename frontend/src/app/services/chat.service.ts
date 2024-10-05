import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';

import { ChatMessage } from '../models/chat-message';

@Injectable({
  providedIn: 'root',
})
export class ChatService {
  private chatAPIEndpoint =
    'REMOVED'; // generalist chatder agent

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
      },
    );
  }
}
