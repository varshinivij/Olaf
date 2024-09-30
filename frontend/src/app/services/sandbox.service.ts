import { Injectable } from '@angular/core';
import { HttpClient, HttpHeaders } from '@angular/common/http';
import { Observable } from 'rxjs';

@Injectable({
  providedIn: 'root'
})
export class SandboxService {
  private sandboxId: string | null = null;

  constructor(private http: HttpClient) {}

  boxRequestApi:string = 'REMOVED'
  boxExecuteApi:string = 'REMOVED'
  boxCloseApi: string = 'REMOVED'
  boxStatusApi: string = 'REMOVED'
  boxUploadApi: string = 'REMOVED'
  addFirebaseFileApi: string = 'REMOVED'

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
    return this.http.post<any>(this.boxExecuteApi, { sandboxId: this.sandboxId, code }, { headers });
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

}
