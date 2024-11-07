import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders, HttpParams } from '@angular/common/http';
import { Observable } from 'rxjs';

@Injectable({
  providedIn: 'root',
})
export class ChatService {
  private chatAPIEndpoint =
    'http://127.0.0.1:5001/twocube-web/us-central1/master_agent_interaction';

  constructor(private http: HttpClient) {}

  sendMessage(
    message: string,
    sessionId: string,
    userId: string,
    projectId: string,
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
}
