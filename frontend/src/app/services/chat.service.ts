import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { ChatMessage } from '../models/chat-message';

@Injectable({
  providedIn: 'root'
})
export class ChatService {

  private chatAPIEndpoint = 'https://ask-agent-7drpntdska-uc.a.run.app'; // generalist chat
  private plannerAPIEndpoint = 'https://generate-plan-7drpntdska-uc.a.run.app'
  private coderAPIEndpoint = 'https://generate-code-7drpntdska-uc.a.run.app'; // coder agent

  constructor(private http: HttpClient) { }

  sendMessage(history: ChatMessage[]): Observable<ChatMessage[]> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<ChatMessage[]>(this.chatAPIEndpoint, { history }, { headers });
  }

  requestPlan(history: ChatMessage[]): Observable<ChatMessage[]> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<ChatMessage[]>(this.plannerAPIEndpoint, { history }, { headers });
  }

  requestCode(history: ChatMessage[]): Observable<ChatMessage[]> {
    // remove all images from history
    history = history.filter((message) => message.type !== 'image');
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<ChatMessage[]>(this.coderAPIEndpoint, { history }, { headers });
  }
  

}