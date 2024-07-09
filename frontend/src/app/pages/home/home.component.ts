import { Component } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { ChatService } from '../../services/chat.service';
import { ChatMessage } from '../../models/chat-message';
@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss']
})
export class HomeComponent {


  constructor(
    private http: HttpClient,
    private chatService: ChatService) { }

  messages:ChatMessage[] = [
    { type: 'text', role: 'assistant', content: 'Hello, how can I help you today?' },
  ];

  newMessage: string = '';

  sendMessage() {
    if (this.newMessage.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: this.newMessage
      };
      this.messages.push(userMessage);
      this.newMessage = '';

      this.chatService.sendMessage(this.messages).subscribe(
        (response: ChatMessage[]) => {
          this.messages = [...this.messages, ...response];
        },
        error => {
          console.error('Error:', error);
        }
      );
    }
  }
}