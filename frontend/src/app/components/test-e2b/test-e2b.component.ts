import { Component } from '@angular/core';
import { CommonModule, JsonPipe, KeyValuePipe } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient, HttpHeaders } from '@angular/common/http';

import { saveAs } from 'file-saver';
import { basename } from 'path-browserify';

@Component({
  selector: 'app-test-e2b',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './test-e2b.component.html',
  styleUrl: './test-e2b.component.scss',
})
export class TestE2bComponent {
  sandboxId = '';
  results = [];
  selectedUploadFiles: File[] = [];
  downloadPath = '';

  executeUrl =
    'https://us-central1-twocube-web.cloudfunctions.net/execute_on_sandbox';
  boxUrl = 'https://us-central1-twocube-web.cloudfunctions.net/request_sandbox';
  uploadUrl =
    'https://us-central1-twocube-web.cloudfunctions.net/upload_to_sandbox';
  downloadUrl =
    'https://us-central1-twocube-web.cloudfunctions.net/download_from_sandbox';

  constructor(private http: HttpClient) {}

  getSandboxId() {
    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });
    this.http.post<any>(this.boxUrl, { headers }).subscribe((response: any) => {
      console.log(response);
      this.sandboxId = response.sandboxId;
    });
  }

  postCode() {
    const code = `import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(0, 10, 100)
y = np.sin(x)

plt.figure(figsize=(10, 6))
plt.plot(x, y, label='sin(x)')

plt.xlabel('x')
plt.ylabel('y')
plt.title('Simple Sine Wave')

plt.legend()

plt.grid(True)

plt.show()`;

    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });

    this.http
      .post<any>(
        this.executeUrl,
        { sandboxId: this.sandboxId, code },
        { headers }
      )
      .subscribe((response: any) => {
        console.log(response);
        console.log(this.results);
        this.results = response.results;
      });
  }

  onUploadFilesSelected(event: Event) {
    const target = event.target as HTMLInputElement;
    const files = target.files as FileList;
    this.selectedUploadFiles = Array.from(files);
  }

  uploadFiles() {
    // making the request as formData will send
    // as contentType="multipart/form-data"
    // however, explicity setting the header seems
    // to be bad. reference: https://github.com/pallets/flask/discussions/4901
    let formData = new FormData();
    formData.append('sandboxId', this.sandboxId);
    for (let file of this.selectedUploadFiles) {
      formData.append('file', file, file.name);
    }

    this.http.post(this.uploadUrl, formData).subscribe((response: any) => {
      console.log(response);
    });
  }

  downloadFile() {
    const path = this.downloadPath;

    const headers = new HttpHeaders({
      'Content-Type': 'application/json',
    });

    this.http
      .post(
        this.downloadUrl,
        { sandboxId: this.sandboxId, path },
        { headers, responseType: 'blob' }
      )
      .subscribe((response: any) => {
        saveAs(response, basename(path));
      });
  }
}
