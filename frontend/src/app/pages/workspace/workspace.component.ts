import { Component, OnDestroy, OnInit } from '@angular/core';
import { AsyncPipe, CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { Router } from '@angular/router';
import { Observable, Subscription } from 'rxjs';
import Split from 'split.js';

import { ChatService } from '../../services/chat.service';
import { SandboxService } from '../../services/sandbox.service';
import { UploadService } from '../../services/upload.service';
import { UserService } from '../../services/user.service';
import { SessionsService } from '../../services/sessions.service';
import { ChatMessage } from '../../models/chat-message';


@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './workspace.component.html',
  styleUrls: ['./workspace.component.scss'],
})
export class WorkspaceComponent implements OnInit, OnDestroy {
  loading: boolean = false;
  isConnected: boolean = false;
  selectedUploadFiles: File[] = [];
  uploadSubscription: Subscription | undefined;
  
  selectedTab: string = 'planner'; // Default tab

  // message scheme ONLY FOR REFERENCE
  // messages: ChatMessage[] = [
  //   {
  //     type: 'text',
  //     role: 'assistant',
  //     content: 'Hello, how can I help you today?',
  //   },
  // ];

  plans: string[] = []
  codes: string[] = []

  newMessage: string = '';

  constructor(
    private http: HttpClient,
    private router: Router,
    private chatService: ChatService,
    private sandboxService: SandboxService,
    public sessionsService: SessionsService,
    public uploadService: UploadService,
    public userService: UserService
  ) {}

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content', '#den-sidebar'], {
      sizes: [25, 50, 25], // Initial sizes of the columns in percentage
      minSize: 200, // Minimum size of each column in pixels
      gutterSize: 10, // Size of the gutter (the draggable area between columns)
      cursor: 'col-resize', // Cursor to show when hovering over the gutter
    });
  }

  ngOnInit() {
    this.uploadSubscription = this.uploadService.getUploadProgress().subscribe(
      (uploads) => {
        console.log(uploads);
      }
    ); 
  }

  ngOnDestroy() {
    this.uploadSubscription?.unsubscribe();
  }

  selectTab(tab: string) {
    this.selectedTab = tab;
  }

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
      this.sessionsService.addMessageToActiveSession(userMessage);
      this.newMessage = '';
      this.loading = true;

      this.chatService.sendMessage(this.sessionsService.activeSession.history).subscribe(
        (responnse:any) => {
          let responseMessages: ChatMessage = {
            type: 'text',
            role: 'assistant',
            content: responnse["message"],
          };
          this.sessionsService.activeSession.history.push(responseMessages);
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
    this.chatService.sendMessage(this.sessionsService.activeSession.history).subscribe(
      (response: ChatMessage[]) => {
        this.sessionsService.activeSession.history = [...this.sessionsService.activeSession.history, ...response];
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
        this.sessionsService.activeSession.history.push(codeResultMessage);
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

  requestPlan() {
    this.loading = true;
    this.chatService.requestPlan(this.sessionsService.activeSession.history).subscribe(
      (response: any) => {
        //add the plan to message as a ChatMessage with type Plan
        let planMessage: ChatMessage = {
          type: 'plan',
          role: 'assistant',
          content: response["message"],
        }
        this.sessionsService.activeSession.history.push(planMessage);

      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  requestCode() {
    console.log(this.sessionsService.activeSession.history);
    this.loading = true;
    this.chatService.requestCode(this.sessionsService.activeSession.history).subscribe(
      (response: any) => {
        //add the code to message as a ChatMessage with type Code
        let codeMessage: ChatMessage = {
          type: 'code',
          role: 'assistant',
          content: response["code"],
        }
        this.sessionsService.activeSession.history.push(codeMessage);
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }
}