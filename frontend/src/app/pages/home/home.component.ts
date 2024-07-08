import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent {
  messages = [
    { type: 'text', content: 'Hello, how can I help you today?' },
    { type: 'code', content: 'const x = 10;\nconsole.log(x);' }
  ];

  newMessage: string = '';

  sendMessage() {
    if (this.newMessage.trim()) {
      this.messages.push({ type: 'text', content: this.newMessage.trim() });
      this.newMessage = '';
    }
  }
}