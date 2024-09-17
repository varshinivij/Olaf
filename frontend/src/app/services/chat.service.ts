import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { ChatMessage } from '../models/chat-message';

@Injectable({
  providedIn: 'root',
})
export class ChatService {
  private chatAPIEndpoint =
    'http://127.0.0.1:5001/twocube-web/us-central1/generate_plan'; // generalist chat
  private plannerAPIEndpoint = 'https://generate-plan-7drpntdska-uc.a.run.app';
  private coderAPIEndpoint = 'https://generate-code-7drpntdska-uc.a.run.app'; // coder agent

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

  requestPlan(history: ChatMessage[]): Observable<ChatMessage[]> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });
    return this.http.post<ChatMessage[]>(
      this.plannerAPIEndpoint,
      { history },
      { headers }
    );
  }

  requestCode(history: ChatMessage[]): Observable<ChatMessage[]> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });
    return this.http.post<ChatMessage[]>(
      this.coderAPIEndpoint,
      { history },
      { headers }
    );
  }
}
