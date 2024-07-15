import { Component } from '@angular/core';
import { AsyncPipe, CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { Router } from '@angular/router';
import { Observable } from 'rxjs';

import { ChatService } from '../../services/chat.service';
import { SandboxService } from '../../services/sandbox.service';
import { UserService } from '../../services/user.service';

import { ChatMessage } from '../../models/chat-message';
import { User } from '../../models/user';

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './home.component.html',
  styleUrls: ['./home.component.scss'],
})
export class HomeComponent {
  loading: boolean = false;
  isConnected: boolean = false;
  user$: Observable<User | null>;

  constructor(
    private http: HttpClient,
    private router: Router,
    private chatService: ChatService,
    private sandboxService: SandboxService,
    private userService: UserService,
  ) {
    // example of how to use userService observable
    this.user$ = userService.getCurrentUser();
  }

  messages: ChatMessage[] = [
    {
      type: 'text',
      role: 'assistant',
      content: 'Hello, how can I help you today?',
    },
  ];

  newMessage: string = '';

  connectToSandBox() {
    this.loading = true;
    this.sandboxService.createSandbox().subscribe(
      (response: any) => {
        this.sandboxService.setSandboxId(response.sandboxId);
        this.isConnected = true;
        this.loading = false;
        console.log(this.sandboxService.getSandboxId());
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  sendMessage() {
    if (this.newMessage.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: this.newMessage,
      };
      this.messages.push(userMessage);
      this.newMessage = '';
      this.loading = true;

      this.chatService.sendMessage(this.messages).subscribe(
        (response: ChatMessage[]) => {
          this.processResponse(response);
          this.messages = [...this.messages, ...response];
          this.loading = false;
        },
        (error) => {
          console.error('Error:', error);
          this.loading = false;
        }
      );
    }
  }

  continue() {
    this.loading = true;
    this.chatService.sendMessage(this.messages).subscribe(
      (response: ChatMessage[]) => {
        this.messages = [...this.messages, ...response];
        this.loading = false;
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  executeCode(code: string) {
    this.sandboxService.executeCode(code).subscribe(
      (result: any) => {
        const codeResultMessage: ChatMessage = {
          type: 'text',
          role: 'assistant',
          content: result.output || 'Code executed goes here',
        };
        this.messages.push(codeResultMessage);
        this.loading = false;
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  processResponse(response: ChatMessage[]) {
    response.forEach((message) => {
      if (message.type === 'code') {
        this.executeCode(message.content);
      }
    });
  }

  async logout() {
    try {
      await this.userService.logout();
    } catch (error) {
      console.error('Error logging out: ', error);
    }
  }

  async deleteAccount() {
    try {
      await this.userService.deleteAccount();
    } catch (error) {
      console.error('Error deleting account: ', error);
    }
  }
}
