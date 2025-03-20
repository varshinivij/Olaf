import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';
import { throwError } from 'rxjs';
import { environment } from '../../environments/environment'
@Injectable({
  providedIn: 'root'
})
export class SandboxService {
  private sandboxId: string | null = null;

  constructor(private http: HttpClient) {}

  boxRequestApi:string = environment.boxRequestApi;
  boxExecuteApi:string = environment.boxExecuteApi;
  boxCloseApi: string = environment.boxCloseApi
  boxStatusApi: string = environment.boxStatusApi
  boxUploadApi: string = environment.boxUploadApi
  addFirebaseFileApi: string = environment.addFirebaseFileApi
  boxTerminalApi: string = environment.boxTerminalApi

  createSandbox(): Observable<any> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<any>(this.boxRequestApi, { headers });
  }

  setSandboxId(id: string) {
    this.sandboxId = id;
  }

  getSandboxId(): string | null {
    return this.sandboxId;
  }
  
  closeSandbox(): Observable<any> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<any>(this.boxCloseApi, { sandboxId: this.sandboxId }, { headers });
  }

  executeCode(code: string): Observable<any> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<any>(this.boxExecuteApi, { sandboxId: this.sandboxId, code });
  }

  isSandboxConnected(){
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<any>(this.boxStatusApi, { sandboxId: this.sandboxId }, { headers });
  }

  uploadFile(file: File|Blob): Observable<any> {
    const formData = new FormData();
    formData.append('file', file);
    formData.append('sandboxId', this.sandboxId || '');
    return this.http.post<any>(this.boxUploadApi, formData);
  }

  addFirebaseFilesToSandBox(filepaths: [string]): Observable<any> {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json'
    });
    return this.http.post<any>(this.addFirebaseFileApi,  { sandboxId: this.sandboxId, filePaths: filepaths}, { headers });
  }

  runTerminalCommand(command: string): Observable<any> {
    const sandboxId = this.getSandboxId();
    if (!sandboxId) {
      return throwError(() => new Error('No sandbox ID. Sandbox not ready.'));
    }
    return this.http.post<any>(
      this.boxTerminalApi,
      { sandboxId, command }
    );
  }

}
