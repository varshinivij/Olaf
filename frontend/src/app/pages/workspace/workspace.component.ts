import { Component, OnDestroy, OnInit } from '@angular/core';
import { AsyncPipe, CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient, HttpHeaders, HttpEventType } from '@angular/common/http';
import { Router } from '@angular/router';
import { Observable, Subscription } from 'rxjs';
import Split from 'split.js';
import { jsonrepair } from 'jsonrepair';

import { ChatService } from '../../services/chat.service';
import { SandboxService } from '../../services/sandbox.service';
import { FileStorageService } from '../../services/file-storage.service';
import { UploadService } from '../../services/upload.service';
import { UserService } from '../../services/user.service';
import { SessionsService } from '../../services/sessions.service';
import { ChatMessage } from '../../models/chat-message';
import { UserFile } from '../../models/user-file';
import { response } from 'express';


interface CodeStream {
  isOpen: boolean; // Indicates if a code block has started
  buffer: string;  // Temporary buffer to hold incomplete code
}


const codeStream: CodeStream = {
  isOpen: false,
  buffer: ''
};

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './workspace.component.html',
  styleUrls: ['./workspace.component.scss'],
})
export class WorkspaceComponent implements OnInit, OnDestroy {
  loading: boolean = false;
  executingCode: boolean = false;
  isConnected: boolean = false;
  userFiles: UserFile[] = [];
  selectedUploadFiles: File[] = [];
  uploadSubscription: Subscription | undefined;

  selectedTab: string = 'planner'; // Default tab
  public latestPlanMessage: any;

  // message scheme ONLY FOR REFERENCE
  // messages: ChatMessage[] = [
  //   {
  //     type: 'text',
  //     role: 'assistant',
  //     content: 'Hello, how can I help you today?',
  //   },
  // ];

  plans: string[] = [];
  codes: string[] = [];

  newMessage: string = '';

  constructor(
    private http: HttpClient,
    private router: Router,
    private chatService: ChatService,
    private sandboxService: SandboxService,
    public sessionsService: SessionsService,
    public uploadService: UploadService,
    public userService: UserService,
    private fileStorageService: FileStorageService
  ) { }

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content', '#den-sidebar'], {
      sizes: [25, 50, 25], // Initial sizes of the columns in percentage
      minSize: 200, // Minimum size of each column in pixels
      gutterSize: 10, // Size of the gutter (the draggable area between columns)
      cursor: 'col-resize', // Cursor to show when hovering over the gutter
    });
  }

  ngOnInit() {
    this.uploadSubscription = this.uploadService
      .getUploadProgress()
      .subscribe((uploads) => {
        console.log(uploads);
      });
    this.getUserFiles();
    this.getLatestPlanMessage();
  }

  ngOnDestroy() {
    this.uploadSubscription?.unsubscribe();
  }

  selectTab(tab: string) {
    this.selectedTab = tab;
  }

  getUserFiles() {
    this.fileStorageService.getFiles().subscribe((files) => {
      console.log(files);
      this.userFiles = files || [];
    });
  }

  /**
   * Add a file to the e2b sandbox using the firebase storage link
   * @param fileUrl storage download link url of the file
   */
  addFileToSandbox(file: UserFile) {
    const uploadText = `Uploaded ${file.name}`;
    const uploadMessage: ChatMessage = {
      type: 'text',
      role: 'user',
      content: uploadText,
    };
    this.sessionsService.addMessageToActiveSession(uploadMessage);
    // this.loading = true;
    // const downloadUrl = this.fileStorageService.getDownloadUrl(file);
    // this.sandboxService.uploadFile(downloadUrl).subscribe(
    //   (response: any) => {
    //     console.log(response);
    //     this.loading = false;
    //   },
    //   (error) => {
    //     console.error('Error:', error);
    //     this.loading = false;
    //   }
    // );
  }

  connectToSandBox() {
    this.loading = true;
    if (this.sessionsService.activeSession.sandboxId) {
      this.checkSandboxConnection();
    } else {
      this.createSandbox();
    }
    // Destroy the sandbox when the tab is closed
    window.addEventListener('unload', () => this.sandboxService.closeSandbox());
  }

  private checkSandboxConnection() {
    this.sandboxService.isSandboxConnected().subscribe(
      (response: any) => {
        if (response.alive) {
          this.onSandboxConnected();
        } else {
          this.createSandbox();
        }
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  private onSandboxConnected() {
    this.isConnected = true;
    this.loading = false;
    if (this.sessionsService.activeSession.sandboxId) {
      this.sandboxService.setSandboxId(
        this.sessionsService.activeSession.sandboxId
      );
    }
  }

  private createSandbox() {
    this.sandboxService.createSandbox().subscribe(
      (response: any) => {
        this.sandboxService.setSandboxId(response.sandboxId);
        this.onSandboxCreated();
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      }
    );
  }

  private onSandboxCreated() {
    this.isConnected = true;
    this.loading = false;
    console.log(this.sandboxService.getSandboxId());
  }

  /**
   * Add a new message to the active session from message bar
   **/
  async addUserMessage(): Promise<void> {
    if (this.newMessage.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: this.newMessage,
      };
      await this.sessionsService.addMessageToActiveSession(userMessage);
      this.newMessage = '';
    }
  }

  /**
   * Send a message to the generalist chat service
   */
  async sendMessage() {
    if (!this.isConnected) {
      this.connectToSandBox();
    }
    if (this.newMessage.trim()) {
      await this.addUserMessage();
      console.log(this.sessionsService.activeSession.history);
      this.loading = true;
      let responseType:
        | 'text'
        | 'code'
        | 'plan'
        | 'error'
        | 'image'
        | 'result' = 'text'; // Default type
      let responseContent = ''; // To store the content from the rest of the chunks

      // Create an initial placeholder message in the chat
      const responseMessage: ChatMessage = {
        type: responseType,
        role: 'assistant',
        content: '', // Initially empty
      };
      // Add the placeholder message to the session and store a reference to it
      await this.sessionsService.addMessageToActiveSession(responseMessage);
      const messageIndex =
        this.sessionsService.activeSession.history.length - 1; // Get the index of the added message

      this.chatService
        .sendMessage(this.sessionsService.activeSession.history)
        .subscribe({
          next: (event: any) => {
            if (event.type === HttpEventType.DownloadProgress && event.loaded) {
              const chunk = event.partialText || ''; // Handle the chunked text

              if (responseContent === '') {
                // First chunk is the type
                responseType = chunk;
                responseMessage.type = responseType;
                responseContent += ' ';
              } else {
                // Append subsequent chunks to the content
                responseContent = chunk;
                if (
                  responseType === 'text' ||
                  responseType === 'code' ||
                  responseType === 'plan'
                ) {
                  // Remove the first four letters from the content (the "code" identifier)
                  responseContent = responseContent.slice(4).trim();
                }
                if (responseType === 'plan') {
                  responseContent = jsonrepair(responseContent);
                  const jsonData = JSON.parse(responseContent);
                  responseContent = JSON.stringify(jsonData, null, 2);
                }
                this.sessionsService.activeSession.history[
                  messageIndex
                ].content = responseContent;
              }
            }
          },
          error: (error) => {
            console.error('Error:', error);
            this.loading = false;
          },
          complete: () => {
            this.loading = false;
            this.getLatestPlanMessage();
            this.executeLatestCode();
          },
        });
    }
  }

  executeLatestCode(): void {
    const history = this.sessionsService.activeSession.history;
    const latestCodeMessage = history
      .slice()
      .reverse()
      .find((message) => message.type === 'code');
    if (latestCodeMessage) {
      let code = this.extractCode(latestCodeMessage.content);
      this.executeCode(code);
    }
  }



  getLatestPlanMessage(): void {
    const history = this.sessionsService.activeSession.history;
    // Find the last message of type 'plan'
    this.latestPlanMessage = history
      .slice()
      .reverse()
      .find((message) => message.type === 'plan');
    console.log(this.latestPlanMessage);
  }



  parseJson(content: string): any {
    try {
      return JSON.parse(content);
    } catch (e) {
      console.error('Error parsing JSON:', e);
      return null; // Return null if parsing fails
    }
  }



  continue() {
    this.loading = true;
    this.chatService
      .sendMessage(this.sessionsService.activeSession.history)
      .subscribe(
        (response: ChatMessage[]) => {
          this.sessionsService.activeSession.history = [
            ...this.sessionsService.activeSession.history,
            ...response,
          ];
          this.loading = false;
        },
        (error) => {
          console.error('Error:', error);
          this.loading = false;
        }
      );
  }

  executeCode(code: string) {
    this.executingCode = true;
    this.sandboxService.executeCode(code).subscribe(
      (result: any) => {
        console.log(result);
        // if result stdout is not empty, add it to the chat
        if (result.logs.stdout && result.logs.stdout.length > 0) {
          const stdoutContent = result.logs.stdout.join('\n');
          const codeResultMessage: ChatMessage = {
            type: 'result',
            role: 'assistant',
            content: stdoutContent,
          };
          this.sessionsService.addMessageToActiveSession(codeResultMessage);
        }
        if (result.results && result.results.length > 0) {
          if (result.results[0]['image/png']) {
            const base64Image = `data:image/png;base64,${result.results[0]['image/png']}`;
            const imageMessage: ChatMessage = {
              type: 'image',
              role: 'assistant',
              content: base64Image,
            };
            this.sessionsService.addMessageToActiveSession(imageMessage);
          }
        }
        if (result.error && result.error.length > 0) {
          const errorMessage: ChatMessage = {
            type: 'error',
            role: 'assistant',
            content: result.error,
          };
          this.sessionsService.addMessageToActiveSession(errorMessage);
        }
        this.executingCode = false;
      },
      (error) => {
        console.error('Error:', error);
        this.executingCode = false;
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
    this.addUserMessage();
    this.loading = true;
    this.chatService
      .requestPlan(this.sessionsService.activeSession.history)
      .subscribe(
        (response: any) => {
          //add the plan to message as a ChatMessage with type Plan
          let planMessage: ChatMessage = {
            type: 'plan',
            role: 'assistant',
            content: response['message'],
          };
          this.sessionsService.addMessageToActiveSession(planMessage);
        },
        (error) => {
          console.error('Error:', error);
          this.loading = false;
        }
      );
  }

  requestCode(withExecute: boolean = true) {
    if (!this.isConnected) {
      this.connectToSandBox();
    }
    this.addUserMessage();
    this.loading = true;
    this.chatService
      .requestCode(this.sessionsService.activeSession.history)
      .subscribe(
        (response: any) => {
          //add the code to message as a ChatMessage with type Code
          let code = this.extractCode(response['message']);
          let codeMessage: ChatMessage = {
            type: 'code',
            role: 'assistant',
            content: this.extractCode(response['message']), //TODO this is finicky because the agent is not always returning the code in the same format,
          };
          this.sessionsService.addMessageToActiveSession(codeMessage);
          console.log(this.sessionsService.activeSession.history);
          this.loading = false;
          if (withExecute) {
            this.executeCode(code);
          }
        },
        (error) => {
          console.error('Error:', error);
          this.loading = false;
        }
      );
  }


  convertNewlinesToBr(text: string): string {
    return text.replace(/\n/g, '<br>');
  }

  extractCode(response: string): string {
    // Static variables to store the state between function calls
    let isOpen = false;
    let buffer = '';

    // Regular expression to match complete and incomplete code blocks
    const codeRegex = /```(?:python\s*\n)?([\s\S]*?)(```|$)/g;
    let match;
    let codeParts: string[] = [];

    while ((match = codeRegex.exec(response)) !== null) {
      if (match[2] === '```') {
        // Complete code block found
        codeParts.push((buffer + match[1].trim()).trim());
        isOpen = false; // Reset state
        buffer = '';    // Clear the buffer
      } else {
        // Incomplete code block found
        isOpen = true;
        buffer += match[1].trim() + '\n'; // Append to buffer
      }
    }

    // If the stream is still open (no closing ```), return buffer for incomplete code
    if (isOpen) {
      return buffer;
    }

    // Return the complete code parts
    return codeParts.join('\n\n');
  }





  extractTextWithoutCode(response:string) {
    let isInCodeBlock = false; // To track whether we're inside a code block
    let result = ''; // To store the processed text
    const lines = response.split('\n'); // Split response by lines to process them one by one

    for (let line of lines) {
        if (line.startsWith('```')) {
            // Toggle the code block state when encountering ```
            isInCodeBlock = !isInCodeBlock;
            continue; // Skip the line containing ```
        }

        if (!isInCodeBlock) {
            // Process text lines outside of code blocks
            result += line + '\n';
        }
    }

    // Bold formatting for headers
    result = result.replace(/\*\*(.*?)\*\*/g, '<b>$1</b>'); // Bold formatting
    result = result.replace(/###\s*(.*?)(\n|$)/g, '<b>$1</b>'); // H3 style
    result = result.replace(/##\s*(.*?)(\n|$)/g, '<b>$1</b>');  // H2 style
    result = result.replace(/#\s*(.*?)(\n|$)/g, '<b>$1</b>');   // H1 style

    // Replace multiple newlines with <br/> for better line break handling in HTML
    result = result.replace(/\n{2,}/g, '<br/><br/>');

    return '<br/>' + result.trim();
}
}
