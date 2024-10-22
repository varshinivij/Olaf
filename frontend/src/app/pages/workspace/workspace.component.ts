import { AfterViewInit, Component, ElementRef, ViewChild } from '@angular/core';
import { CommonModule } from '@angular/common';
import { HttpEventType } from '@angular/common/http';
import { FormsModule } from '@angular/forms';
import { Router } from '@angular/router';
import { firstValueFrom } from 'rxjs';

import { jsonrepair } from 'jsonrepair';
import { Highlight } from 'ngx-highlightjs';
import Split from 'split.js';

import { ChatService } from '../../services/chat.service';
import { FileStorageService } from '../../services/file-storage.service';
import { SandboxService } from '../../services/sandbox.service';
import { SessionsService } from '../../services/sessions.service';

import { ChatMessage } from '../../models/chat-message';
import { getLucideIconFromType } from '../../models/extension-type';
import { Project } from '../../models/project';
import { Session } from '../../models/session';
import { UserFile } from '../../models/user-file';

import { CodeMessagePipe } from '../../pipes/codemessage.pipe';
import { PlanMessagePipe } from '../../pipes/planmessage.pipe';

import { adjustTextareaHeight } from '../../utils/adjust-textarea-height';
import { delay } from '../../utils/time-utils';

// icon imports
import { provideIcons } from '@ng-icons/core';
import {
  lucideArrowUpFromLine,
  lucideCheck,
  lucideCircleStop,
  lucideCode,
  lucideEllipsisVertical,
  lucideFileArchive,
  lucideFileChartColumn,
  lucideFileCode,
  lucideFileQuestion,
  lucideFileText,
  lucideFolder,
  lucideHouse,
  lucideLoaderCircle,
  lucidePanelLeftClose,
  lucidePanelLeftOpen,
  lucidePencil,
  lucidePlus,
  lucideRotateCw,
  lucideSendHorizontal,
  lucideSettings,
  lucideTrash2,
} from '@ng-icons/lucide';

// ui imports
import { HlmButtonDirective } from '@spartan-ng/ui-button-helm';
import {
  BrnDialogCloseDirective,
  BrnDialogContentDirective,
  BrnDialogTriggerDirective,
} from '@spartan-ng/ui-dialog-brain';
import {
  HlmDialogComponent,
  HlmDialogContentComponent,
  HlmDialogFooterComponent,
  HlmDialogHeaderComponent,
  HlmDialogTitleDirective,
} from '@spartan-ng/ui-dialog-helm';
import { HlmIconComponent } from '@spartan-ng/ui-icon-helm';
import { HlmInputDirective } from '@spartan-ng/ui-input-helm';
import { BrnMenuTriggerDirective } from '@spartan-ng/ui-menu-brain';
import {
  HlmMenuComponent,
  HlmMenuItemDirective,
  HlmMenuItemIconDirective,
} from '@spartan-ng/ui-menu-helm';
import { HlmScrollAreaComponent } from '@spartan-ng/ui-scrollarea-helm';
import { HlmSeparatorDirective } from '@spartan-ng/ui-separator-helm';
import { HlmSpinnerComponent } from '@spartan-ng/ui-spinner-helm';
import {
  HlmTabsComponent,
  HlmTabsContentDirective,
  HlmTabsListComponent,
  HlmTabsTriggerDirective,
} from '@spartan-ng/ui-tabs-helm';
import { BrnToggleDirective } from '@spartan-ng/ui-toggle-brain';
import { HlmToggleDirective } from '@spartan-ng/ui-toggle-helm';
import {
  HlmCodeDirective,
  HlmH2Directive,
  HlmH3Directive,
  HlmLargeDirective,
  HlmMutedDirective,
  HlmPDirective,
  HlmSmallDirective,
  HlmUlDirective,
} from '@spartan-ng/ui-typography-helm';

@Component({
  selector: 'app-chat',
  standalone: true,
  templateUrl: './workspace.component.html',
  styleUrls: ['./workspace.component.scss'],
  imports: [
    CommonModule,
    FormsModule,
    CodeMessagePipe,
    PlanMessagePipe,
    Highlight,

    HlmButtonDirective,

    BrnDialogContentDirective,
    BrnDialogCloseDirective,
    BrnDialogTriggerDirective,
    HlmDialogComponent,
    HlmDialogContentComponent,
    HlmDialogFooterComponent,
    HlmDialogHeaderComponent,
    HlmDialogTitleDirective,

    HlmIconComponent,
    HlmInputDirective,

    BrnMenuTriggerDirective,
    HlmMenuComponent,
    HlmMenuItemDirective,
    HlmMenuItemIconDirective,

    HlmScrollAreaComponent,
    HlmSeparatorDirective,
    HlmSpinnerComponent,

    HlmTabsComponent,
    HlmTabsContentDirective,
    HlmTabsListComponent,
    HlmTabsTriggerDirective,

    BrnToggleDirective,
    HlmToggleDirective,

    HlmCodeDirective,
    HlmH2Directive,
    HlmH3Directive,
    HlmLargeDirective,
    HlmMutedDirective,
    HlmPDirective,
    HlmSmallDirective,
    HlmUlDirective,
  ],
  providers: [
    provideIcons({
      lucideArrowUpFromLine,
      lucideCheck,
      lucideCircleStop,
      lucideCode,
      lucideEllipsisVertical,
      lucideFileArchive,
      lucideFileChartColumn,
      lucideFileCode,
      lucideFileQuestion,
      lucideFileText,
      lucideFolder,
      lucideHouse,
      lucideLoaderCircle,
      lucidePanelLeftOpen,
      lucidePanelLeftClose,
      lucidePencil,
      lucidePlus,
      lucideRotateCw,
      lucideSendHorizontal,
      lucideSettings,
      lucideTrash2,
    }),
  ],
})
export class WorkspaceComponent implements AfterViewInit {
  @ViewChild('sidebar') sidebar?: ElementRef;
  adjustTextareaHeight = adjustTextareaHeight;
  getLucideIconFromType = getLucideIconFromType;
  split?: Split.Instance;
  collapsed = false; // sidebar

  currentProject: Project;
  currentSession: Session;
  newMessage: string = ''; // ngModel variable
  newSessionName: string = ''; // ngModel variable
  responseLoading: boolean = false;
  executingCode: boolean = false;
  isConnected: boolean = false;

  // a bit ugly of a solution but we're going to rework files soon anyway.
  // this entire component needs serious refactoring on the chat/session side.
  // perhaps some logic should be moved to sessionservice/chatservice.
  uploadedFiles: Set<UserFile['id']> = new Set();

  constructor(
    public router: Router,
    private chatService: ChatService,
    public fileStorageService: FileStorageService,
    private sandboxService: SandboxService,
    public sessionsService: SessionsService,
  ) {
    this.currentProject =
      this.router.getCurrentNavigation()?.extras.state?.['project'];
    if (this.currentProject === undefined) {
      this.router.navigate(['dashboard']);
    }
    this.currentSession = sessionsService.blankSession(this.currentProject);
  }

  /*
    Sidebar UI methods
  */

  ngAfterViewInit() {
    this.split = Split(['#sidebar', '#main-content', '#den-sidebar'], {
      sizes: [20, 45, 35], // initial size of columns (default %)
      minSize: [60, 300, 300], // minimum size of each column (default px)
      gutterSize: 12, // size of the draggable slider (default px)
      snapOffset: 0, // snap to min size offset (default px)
      onDrag: () => {
        // collapsed is true if width smaller than a threshold (px)
        const width = this.sidebar?.nativeElement.offsetWidth;
        this.collapsed = width <= 85;
      },
      onDragStart: () => {
        this.sidebar?.nativeElement.classList.remove('transition-[width]');
      },
      onDragEnd: () => {
        this.sidebar?.nativeElement.classList.add('transition-[width]');
      },
    });
  }

  toggleSidebar() {
    if (!this.collapsed) {
      this.split?.collapse(0);
    } else {
      this.split?.setSizes([20, 45, 35]);
    }
    this.collapsed = !this.collapsed;
  }

  /*
    Session handling methods
  */

  newSession() {
    this.currentSession = this.sessionsService.blankSession(
      this.currentProject,
    );
  }

  changeSession(session: Session) {
    this.currentSession = session;
  }

  renameSession(session: Session) {
    this.sessionsService.renameSession(session, this.newSessionName);
    this.newSessionName = '';
  }

  deleteSession(session: Session) {
    if (session.id === this.currentSession.id) {
      this.newSession();
    }
    this.sessionsService.deleteSession(session);
  }

  /*
    Session file upload/sandbox updating methods
  */

  /**
   * Add a file to the e2b sandbox using the firebase storage link
   * @param fileUrl storage download link url of the file
   */
  async addFirebaseFileToSandbox(file: UserFile) {
    const uploadText = `Uploaded ${file.name}`;
    const uploadMessage: ChatMessage = {
      type: 'text',
      role: 'user',
      content: uploadText,
    };
    this.sessionsService.addMessageToSession(
      this.currentSession,
      uploadMessage,
    );
    this.responseLoading = true;

    if (!this.isConnected) {
      await this.connectToSandbox();
    }

    const filePath = file.storageLink;
    console.log('adding file to sandbox ' + filePath);
    this.sandboxService.addFirebaseFilesToSandBox([filePath]).subscribe(
      (response: any) => {
        console.log('add file response: ', response);
        this.uploadedFiles.add(file.id);
        this.responseLoading = false;
      },
      (error) => {
        console.error('Error:', error);
        this.responseLoading = false;
      },
    );
  }

  async connectToSandbox() {
    try {
      if (this.currentSession.sandboxId) {
        await this.checkSandboxConnection();
      } else {
        await this.createSandbox();
      }
      // TODO find a better solution, but right now we give the e2b box an
      // extra second to get ready to recieve code
      await delay(1000);
    } catch (error) {
      console.error('Error connecting to sandbox:', error);
    }
    window.addEventListener('unload', () => this.sandboxService.closeSandbox());
  }

  private async checkSandboxConnection(): Promise<void> {
    try {
      const response: any = await firstValueFrom(
        this.sandboxService.isSandboxConnected(),
      );
      if (response.alive) {
        this.onSandboxConnected();
      } else {
        this.uploadedFiles.clear();
        await this.createSandbox();
      }
    } catch (error) {
      console.error('Error:', error);
      this.responseLoading = false;
    }
  }

  private onSandboxConnected() {
    this.isConnected = true;
    if (this.currentSession.sandboxId) {
      this.sandboxService.setSandboxId(this.currentSession.sandboxId);
    }
  }

  private async createSandbox() {
    try {
      const response: any = await firstValueFrom(
        this.sandboxService.createSandbox(),
      );
      this.sandboxService.setSandboxId(response.sandboxId);
      this.onSandboxCreated();
    } catch (error) {
      console.error('Error:', error);
    }
  }

  private onSandboxCreated() {
    this.isConnected = true;
    console.log('created sandbox id:', this.sandboxService.getSandboxId());
  }

  /*
    Session chat/executing code methods
  */

  /**
   * Add a new message to the active session from message bar
   **/
  async addUserMessage(message: string): Promise<void> {
    if (message.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: message,
      };
      await this.sessionsService.addMessageToSession(
        this.currentSession,
        userMessage,
      );
    }
  }

  /**
   * Send a message to the generalist chat service
   */
  async sendMessage() {
    if (this.newMessage.trim() && !this.responseLoading) {
      console.log('curr session history: ', this.currentSession.history);

      const message = this.newMessage;
      this.newMessage = '';
      this.responseLoading = true;
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
      };
      // Add the placeholder message to the session and store a reference to it
      await this.sessionsService.addMessageToSession(
        this.currentSession,
        responseMessage,
      );
      const messageIndex = this.currentSession.history.length - 1; // Get the index of the added message

      this.chatService.sendMessage(this.currentSession.history).subscribe({
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
              this.currentSession.history[messageIndex].content =
                responseContent;
            }
          }
        },
        error: (error) => {
          console.error('Error:', error);
          this.responseLoading = false;
        },
        complete: () => {
          this.responseLoading = false;
          this.sessionsService.saveSession(this.currentSession);
          this.executeLatestCode();
        },
      });
    }
  }

  async executeLatestCode() {
    const history = this.currentSession.history;
    const latestCodeMessage = history
      .slice()
      .reverse()
      .find((message) => message.type === 'code');
    if (latestCodeMessage) {
      this.executingCode = true;

      if (!this.isConnected) {
        await this.connectToSandbox();
      }

      let code = this.extractCode(latestCodeMessage.content);
      this.executeCode(code);
    }
  }

  parseJson(content: string): any {
    try {
      return JSON.parse(content);
    } catch (e) {
      console.error('Error parsing JSON:', e);
      return null; // Return null if parsing fails
    }
  }

  executeCode(code: string) {
    this.sandboxService.executeCode(code).subscribe(
      (result: any) => {
        console.log('code result: ', result);
        // if result stdout is not empty, add it to the chat
        if (result.logs.stdout && result.logs.stdout.length > 0) {
          const stdoutContent = result.logs.stdout.join('\n');
          const codeResultMessage: ChatMessage = {
            type: 'result',
            role: 'assistant',
            content: stdoutContent,
          };
          this.sessionsService.addMessageToSession(
            this.currentSession,
            codeResultMessage,
          );
        }
        if (result.results && result.results.length > 0) {
          if (result.results[0]['image/png']) {
            const base64Image = `data:image/png;base64,${result.results[0]['image/png']}`;
            const imageMessage: ChatMessage = {
              type: 'image',
              role: 'assistant',
              content: base64Image,
            };
            this.sessionsService.addMessageToSession(
              this.currentSession,
              imageMessage,
            );
          }
        }
        if (result.error && result.error.length > 0) {
          const errorMessage: ChatMessage = {
            type: 'error',
            role: 'assistant',
            content: result.error,
          };
          this.sessionsService.addMessageToSession(
            this.currentSession,
            errorMessage,
          );
        }
        this.executingCode = false;
      },
      (error) => {
        console.error('Error:', error);
        this.executingCode = false;
      },
    );
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
