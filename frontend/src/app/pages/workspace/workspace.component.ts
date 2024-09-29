import { Component, OnDestroy, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient, HttpEventType } from '@angular/common/http';
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
import { firstValueFrom, Subscription } from 'rxjs';
import { delay } from '../../utils/time-utils';
import { CodeStream } from '../../models/code-stream';

import { adjustTextareaHeight } from '../../utils/adjust-textarea-height';

import { HlmButtonDirective } from '@spartan-ng/ui-button-helm';
import { HlmInputDirective } from '@spartan-ng/ui-input-helm';
import { HlmIconComponent } from '@spartan-ng/ui-icon-helm';
import {
  HlmLargeDirective,
  HlmMutedDirective,
  HlmSmallDirective,
  HlmPDirective,
  HlmH1Directive,
  HlmH2Directive,
  HlmH3Directive,
  HlmH4Directive,
  HlmUlDirective,
  HlmCodeDirective,
} from '@spartan-ng/ui-typography-helm';
import { HlmSeparatorDirective } from '@spartan-ng/ui-separator-helm';
import { HlmToggleDirective } from '@spartan-ng/ui-toggle-helm';
import { BrnToggleDirective } from '@spartan-ng/ui-toggle-brain';
import { BrnMenuTriggerDirective } from '@spartan-ng/ui-menu-brain';
import {
  HlmMenuComponent,
  HlmMenuItemDirective,
  HlmMenuItemIconDirective,
} from '@spartan-ng/ui-menu-helm';
import {
  BrnDialogCloseDirective,
  BrnDialogContentDirective,
  BrnDialogTriggerDirective,
} from '@spartan-ng/ui-dialog-brain';
import {
  HlmDialogComponent,
  HlmDialogContentComponent,
  HlmDialogDescriptionDirective,
  HlmDialogFooterComponent,
  HlmDialogHeaderComponent,
  HlmDialogTitleDirective,
} from '@spartan-ng/ui-dialog-helm';
import { HlmScrollAreaComponent } from '@spartan-ng/ui-scrollarea-helm';
import {
  HlmTabsComponent,
  HlmTabsContentDirective,
  HlmTabsListComponent,
  HlmTabsTriggerDirective,
} from '@spartan-ng/ui-tabs-helm';

import { provideIcons } from '@ng-icons/core';
import {
  lucidePlus,
  lucidePencil,
  lucideEllipsisVertical,
  lucideSettings,
  lucideTrash2,
  lucidePanelLeftDashed,
  lucideHouse,
  lucideCircleStop,
  lucideRotateCw,
  lucideSendHorizontal,
} from '@ng-icons/lucide';

import { Highlight, HighlightAuto } from 'ngx-highlightjs';
import { Router } from '@angular/router';
import { PlanMessagePipe } from '../../pipes/planmessage.pipe';
import { CodeMessagePipe } from '../../pipes/codemessage.pipe';
import { Session } from '../../models/session';

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    PlanMessagePipe,
    CodeMessagePipe,

    Highlight,
    HighlightAuto,

    HlmLargeDirective,
    HlmMutedDirective,
    HlmSmallDirective,
    HlmPDirective,
    HlmH1Directive,
    HlmH2Directive,
    HlmH3Directive,
    HlmH4Directive,
    HlmUlDirective,
    HlmCodeDirective,

    HlmButtonDirective,
    HlmIconComponent,
    HlmSeparatorDirective,
    HlmInputDirective,

    HlmToggleDirective,
    BrnToggleDirective,

    BrnMenuTriggerDirective,
    HlmMenuComponent,
    HlmMenuItemDirective,
    HlmMenuItemIconDirective,

    BrnDialogContentDirective,
    BrnDialogTriggerDirective,
    BrnDialogCloseDirective,
    HlmDialogComponent,
    HlmDialogContentComponent,
    HlmDialogDescriptionDirective,
    HlmDialogFooterComponent,
    HlmDialogHeaderComponent,
    HlmDialogTitleDirective,

    HlmScrollAreaComponent,
    HlmTabsComponent,
    HlmTabsContentDirective,
    HlmTabsListComponent,
    HlmTabsTriggerDirective,
  ],
  providers: [
    provideIcons({
      lucidePlus,
      lucidePencil,
      lucideEllipsisVertical,
      lucideSettings,
      lucideTrash2,
      lucidePanelLeftDashed,
      lucideHouse,
      lucideCircleStop,
      lucideRotateCw,
      lucideSendHorizontal,
    }),
  ],
  templateUrl: './workspace.component.html',
  styleUrls: ['./workspace.component.scss'],
})
export class WorkspaceComponent implements OnInit, OnDestroy {
  adjustTextareaHeight = adjustTextareaHeight;
  fileStorageSubscription?: Subscription;

  loading: boolean = false;
  executingCode: boolean = false;
  isConnected: boolean = false;
  userFiles: UserFile[] = [];
  selectedUploadFiles: File[] = [];
  uploadSubscription: Subscription | undefined;

  codeStream: CodeStream = {
    isOpen: false,
    buffer: '',
  };

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
  newSessionName: string = '';

  constructor(
    private http: HttpClient,
    public router: Router,
    private chatService: ChatService,
    private sandboxService: SandboxService,
    public sessionsService: SessionsService,
    public uploadService: UploadService,
    public userService: UserService,
    private fileStorageService: FileStorageService,
  ) {}

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content', '#den-sidebar'], {
      sizes: [20, 45, 35], // Initial sizes of the columns in percentage
      minSize: [200, 200, 300], // Minimum size of each column in pixels
      gutterSize: 12, // Size of the gutter (the draggable area between columns)
      snapOffset: 0,
    });
  }

  ngOnInit() {
    this.uploadSubscription = this.uploadService
      .getUploadProgress()
      .subscribe((uploads) => {
        console.log(uploads);
      });
    this.fileStorageSubscription = this.fileStorageService
      .getFiles()
      .subscribe((files) => {
        console.log(files);
        this.userFiles = files || [];
      });
    this.getLatestPlanMessage();
  }

  ngOnDestroy() {
    this.uploadSubscription?.unsubscribe();
    this.fileStorageSubscription?.unsubscribe();
  }

  renameSession(session: Session) {
    this.sessionsService.renameSession(session, this.newSessionName);
    this.newSessionName = '';
  }

  /**
   * Add a file to the e2b sandbox using the firebase storage link
   * @param fileUrl storage download link url of the file
   */
  async addFirebaseFileToSandbox(file: UserFile) {
    if (!this.isConnected) {
      await this.connectToSandbox();
    }

    const uploadText = `Uploaded ${file.name}`;
    const uploadMessage: ChatMessage = {
      type: 'text',
      role: 'user',
      content: uploadText,
    };
    this.sessionsService.addMessageToActiveSession(uploadMessage);
    this.loading = true;
    const filePath = file.storageLink;
    console.log('adding file to sandbox ' + filePath);
    this.sandboxService.addFirebaseFilesToSandBox([filePath]).subscribe(
      (response: any) => {
        console.log(response);
        this.loading = false;
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      },
    );
  }

  async connectToSandbox() {
    this.loading = true;
    try {
      if (this.sessionsService.activeSession.sandboxId) {
        await this.checkSandboxConnection();
      } else {
        await this.createSandbox();
      }
      // TODO find a better solution, but right now we give the e2b box an
      // extra second to get ready to recieve code
      await delay(1000);
    } catch (error) {
      console.error('Error connecting to sandbox:', error);
    } finally {
      this.loading = false;
    }
    window.addEventListener('unload', () => this.sandboxService.closeSandbox());
  }

  private async checkSandboxConnection(): Promise<void> {
    try {
      const response: any = await firstValueFrom(
        this.sandboxService.isSandboxConnected()
      );
      if (response.alive) {
        this.onSandboxConnected();
      } else {
        await this.createSandbox();
      }
    } catch (error) {
      console.error('Error:', error);
      this.loading = false;
    }
  }

  private onSandboxConnected() {
    this.isConnected = true;
    this.loading = false;
    if (this.sessionsService.activeSession.sandboxId) {
      this.sandboxService.setSandboxId(
        this.sessionsService.activeSession.sandboxId,
      );
    }
  }

  private async createSandbox() {
    try {
      const response: any = await firstValueFrom(
        this.sandboxService.createSandbox()
      );
      this.sandboxService.setSandboxId(response.sandboxId);
      this.onSandboxCreated();
    } catch (error) {
      console.error('Error:', error);
      this.loading = false;
    }
  }

  private onSandboxCreated() {
    this.isConnected = true;
    this.loading = false;
    console.log(this.sandboxService.getSandboxId());
  }

  /**
   * Add a new message to the active session from message bar
   **/
  async addUserMessage(message: string): Promise<void> {
    if (message.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: message,
        isLive: true,
      };
      await this.sessionsService.addMessageToActiveSession(userMessage);
    }
  }

  /**
   * Send a message to the generalist chat service
   */
  async sendMessage() {
    const message = this.newMessage;
    this.newMessage = '';

    if (!this.isConnected) {
      await this.connectToSandbox();
    }
    if (message.trim()) {
      console.log(this.sessionsService.activeSession.history);

      this.loading = true;
      await this.addUserMessage(message);

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
        isLive: true,
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
        },
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
      },
    );
  }

  processResponse(response: ChatMessage[]) {
    response.forEach((message) => {
      if (message.type === 'code') {
        this.executeCode(message.content);
      }
    });
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
        buffer = ''; // Clear the buffer
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

  extractTextWithoutCode(response: string) {
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
    result = result.replace(/##\s*(.*?)(\n|$)/g, '<b>$1</b>'); // H2 style
    result = result.replace(/#\s*(.*?)(\n|$)/g, '<b>$1</b>'); // H1 style

    // Replace multiple newlines with <br/> for better line break handling in HTML
    result = result.replace(/\n{2,}/g, '<br/><br/>');

    return '<br/>' + result.trim();
  }
}
