import { Component } from '@angular/core';
import { CommonModule, JsonPipe, KeyValuePipe } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient, HttpHeaders } from '@angular/common/http';

@Component({
  selector: 'app-test-e2b',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './test-e2b.component.html',
  styleUrl: './test-e2b.component.scss',
})
export class TestE2bComponent {
  sandboxId = '';
  executeResponse: any;
  results = [];

  executeUrl =
    'https://us-central1-twocube-web.cloudfunctions.net/execute_on_sandbox';
  boxUrl = 'https://us-central1-twocube-web.cloudfunctions.net/request_sandbox';

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
        this.executeResponse = response;
        console.log(response);
        this.results = this.executeResponse.results;
        console.log(this.results);
      });
  }
}
