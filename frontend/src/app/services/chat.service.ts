import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { ChatMessage } from '../models/chat-message';

@Injectable({
  providedIn: 'root'
})
export class ChatService {

  private chatAPIEndpoint = 'https://ask-agent-7drpntdska-uc.a.run.app'; // Replace with your cloud function URL

  constructor(private http: HttpClient) { }

  sendMessage(history: ChatMessage[]): Observable<ChatMessage[]> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<ChatMessage[]>(this.chatAPIEndpoint, { history }, { headers });
  }

}