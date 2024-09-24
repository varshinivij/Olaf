import { Component, OnDestroy, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { Router } from '@angular/router';
import { Subscription } from 'rxjs';
import Split from 'split.js';

import { ChatService } from '../../services/chat.service';
import { SandboxService } from '../../services/sandbox.service';
import { FileStorageService } from '../../services/file-storage.service';
import { UploadService } from '../../services/upload.service';
import { UserService } from '../../services/user.service';
import { SessionsService } from '../../services/sessions.service';

import { ChatMessage } from '../../models/chat-message';
import { UserFile } from '../../models/user-file';

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

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,

    HlmLargeDirective,
    HlmMutedDirective,
    HlmSmallDirective,
    HlmPDirective,
    HlmH1Directive,
    HlmH2Directive,
    HlmH3Directive,
    HlmH4Directive,
    HlmUlDirective,

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
  loading: boolean = false;
  executingCode: boolean = false;
  isConnected: boolean = false;
  userFiles: UserFile[] = [];
  selectedUploadFiles: File[] = [];
  uploadSubscription: Subscription | undefined;

  adjustTextareaHeight = adjustTextareaHeight;

  selectedTab: string = 'planner'; // Default tab

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
  console = console;

  constructor(
    private http: HttpClient,
    private router: Router,
    private chatService: ChatService,
    private sandboxService: SandboxService,
    public sessionsService: SessionsService,
    public uploadService: UploadService,
    public userService: UserService,
    private fileStorageService: FileStorageService
  ) {}

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content', '#den-sidebar'], {
      sizes: [20, 45, 35], // Initial sizes of the columns in percentage
      minSize: [200, 200, 300], // Minimum size of each column in pixels
      gutterSize: 12, // Size of the gutter (the draggable area between columns)
      cursor: 'col-resize', // Cursor to show when hovering over the gutter
      snapOffset: 0,
    });
  }

  ngOnInit() {
    this.uploadSubscription = this.uploadService
      .getUploadProgress()
      .subscribe((uploads) => {
        console.log(uploads);
      });
    this.getUserFiles();
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
  addUserMessage(): void {
    if (this.newMessage.trim()) {
      const userMessage: ChatMessage = {
        type: 'text',
        role: 'user',
        content: this.newMessage,
      };
      this.sessionsService.addMessageToActiveSession(userMessage);
      this.newMessage = '';
    }
  }

  /**
   * Send a message to the generalist chat service
   */
  sendMessage() {
    if (this.newMessage.trim()) {
      this.addUserMessage();

      this.loading = true;
      this.chatService
        .sendMessage(this.sessionsService.activeSession.history)
        .subscribe(
          (responnse: any) => {
            let responseMessages: ChatMessage = {
              type: 'text',
              role: 'assistant',
              content: responnse['message'],
            };
            this.sessionsService.addMessageToActiveSession(responseMessages);
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
    console.log(this.sessionsService.activeSession.history);
    this.loading = true;
    this.chatService
      .requestCode(this.sessionsService.activeSession.history)
      .subscribe(
        (response: any) => {
          //add the code to message as a ChatMessage with type Code
          let code = this.getCode(response['message']);
          let codeMessage: ChatMessage = {
            type: 'code',
            role: 'assistant',
            content: this.getCode(response['message']), //TODO this is finicky because the agent is not always returning the code in the same format,
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

  getCode(message: string) {
    console.log(message);
    // Check if the message begins with ```python
    if (message.includes('```python')) {
      let code = message.split('```python\n');
      // If there's any code after splitting, it should end before the next ```
      if (code[1].includes('```')) {
        return code[1].split('```')[0];
      }
      return code[1];
    }
    // If the message contains ``` but not starting with ```python
    else if (message.includes('```')) {
      let code = message.split('```\n');
      // Return the code between the first and second ```
      return code[1];
    }
    return message;
  }
}
