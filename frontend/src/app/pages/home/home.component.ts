import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent {


  constructor(private http: HttpClient) { }

  messages = [
    { type: 'text', role: 'assistant', content: 'Hello, how can I help you today?' },
    { type: 'code', role: 'user', content: 'const x = 10;\nconsole.log(x);' }
  ];

  newMessage: string = '';

  sendMessage() {
    if (this.newMessage.trim()) {
      this.messages.push({ type: 'text', role:'user', content: this.newMessage.trim() });
      this.newMessage = '';
    }
    // run get request on this api https://on-request-example-7drpntdska-uc.a.run.app
    this.http.get('https://on-request-example-7drpntdska-uc.a.run.app', {responseType: 'text'}).subscribe((res: any) => {
      this.messages.push({ type: 'text', role:'assistant', content: res});
    })
  }

}