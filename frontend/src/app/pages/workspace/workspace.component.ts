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

@Component({
  selector: 'app-chat',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,

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
      },
    );
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

  private createSandbox() {
    this.sandboxService.createSandbox().subscribe(
      (response: any) => {
        this.sandboxService.setSandboxId(response.sandboxId);
        this.onSandboxCreated();
      },
      (error) => {
        console.error('Error:', error);
        this.loading = false;
      },
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
          },
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
        },
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
        },
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

  code = `import numpy as np
import matplotlib.pyplot as plt

# Generate x values
x = np.linspace(0, 2 * np.pi, 1000)

# Compute y values (sine of x)
y = np.sin(x)

# Create the plot
plt.plot(x, y)

# Add title and labels
plt.title('Sine Wave')
plt.xlabel('x')
plt.ylabel('sin(x)')

# Show the plot
plt.grid(True)
plt.show()
`;

  error = `  Cell In[2], line 23
    pip install numpy matplotlib
        ^
SyntaxError: invalid syntax`;

  data = [
    'data:image/png;base64,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iOxBqKGJmzZtwvjx4xEVFQWZTIaVK1eaPS9JEl588UVERkbCy8sLo0aNwokTJ674vu+//z7i4uLg6emJAQMGYOfOnTb6DoiIqC29+uqrOH78OFQqFWpra/H3339zbywiInIIQiVitbW16NWrF95///2LPv/666/jnXfewZIlS7Bjxw74+Phg7NixqK+vv+R7fvvtt5g9ezZeeukl7N27F7169cLYsWNRVFRkq2+DiIiIiIhcnLCrJspkMvz888+YMGECAENvWFRUFObMmYOnnnoKgGEz0PDwcCxfvhx33XXXRd9nwIABuPbaa017ouj1esTExODxxx/Hs88+2ybfCxERERERuRanmSOWk5ODgoICsyEnAQEBGDBgALZt23bRRKyhoQF79uzB3LlzTcfkcjlGjRqFbdu2XfKz1Go11Gq16bFer0dZWRmCg4NNG4QSEREREZHrkSQJ1dXViIqKglx+6QGITpOIFRQUAADCw8PNjoeHh5ueu1BJSQl0Ot1FX3Ps2LFLftaCBQswb968qywxERERERE5qzNnzqB9+/aXfN5pErG2NHfuXMyePdv0uLKyEh06dEBOTg78/PzsWDJAo9Fg/fr1GDFiBNzd3e1aFpE9+tU+7Mgpx79v6Y5xPSNMx//5zX78nVWG52/sjNv7drBjCZ0T49c6/rfhJD7efAp39o3G3BsTTMc/3nQS/9t0CjcnhWPezT3sWELnxPi1ji3ZpZj59QEkRPjim4f6m47vyC7Go18fQlywF35+LMWOJXRejOGrV12vwdBFfwMAtj8zDEp3w36FDVo9BizcAABY90QK2vl52auITsuR4re6uhodO3a8Yl7gNIlYRIShslxYWIjIyEjT8cLCQvTu3fuirwkJCYFCoUBhYaHZ8cLCQtP7XYxSqYRSqWx2PCgoCP7+/q0ovfVoNBp4e3sjODjY7kEoMsndG3KlGmEhQQgODjYdV3r7Qa6sh69/oNlxsg7Gr3UovIogV3ojKKidWZz6+JVArvSGl68/49cGGL/W4ZGvgVzpDX//ALM49S/VQq70hruXD+PXRhjDV09XXQ+50hsyGRAZHmqastKg1UOu9AYABAUHIcjP257FdEqOFL/Gz7/SlCWhVk28nI4dOyIiIgJr1641HauqqsKOHTuQknLxljMPDw/07dvX7DV6vR5r16695GvINdRpdAAAz8aWLCKR1DU0xq+b09ziyYXUN95/vXj/JQHVN+gBAJ5uCq4bQFckVI9YTU0NsrKyTI9zcnKwf/9+BAUFoUOHDpg1axb+/e9/o0uXLujYsSNeeOEFREVFmVZWBICRI0fi1ltvxcyZMwEAs2fPxtSpU9GvXz/0798fb731Fmpra3H//fe39bdHDqSeiRgJzBi/SsYvCajp/suGBBJPvZbxSy0nVCK2e/dujBgxwvTYOE9r6tSpWL58Of71r3+htrYWDz/8MCoqKjB48GCsWrUKnp6eptdkZ2ejpKTE9HjSpEkoLi7Giy++iIKCAvTu3RurVq1qtoAHuZZ6TWOLFiuyJKB6rSF+2aNAImJDAomMPbpkCaESseHDh+Ny257JZDLMnz8f8+fPv+Q5p06danZs5syZph4yIoAtsiQ209BEVgRIQHWapqFdRKLh/ZcswVom0UWwRYtEpubQGBKY6f7rwfgl8RhHJLBHl1qCdzmiC0iSZLqRskWLRMSGBBKZaY4Ne8RIQE33X1ax6coYJUQX0Ogk6PSGIbCsCJCIuOoniayeQ7tIYFzsiyzBRIzoAsZKLAB4cmgMCYhzFEhkdaahiYxfEo/x/ssRCdQSrGUSXUDVoAUAuCtkULJHjASkaqwI+CqFWo+JCABQ2xi/PkzESECm+OX9l1qAiRjRBWrVbM0isdWqDY0J7FEgEaka49fbgxVZEk9T/PL+S1fGRIzoAsYeMbZmkaiMQ7t8lKwIkHiMPQrejF8SkKrx/suGBGoJJmJEFzAO62JrFomoQauHRmdYbIYVARJRnWloIuOXxGPsEWNDGLUEEzGiC7BHjERmjF+AjQkkptoGDu0icZl6dNmQQC3ARIzoAsY5YqwEkIiMlQAPNzncFbzFk3hUai52QOJqasxlHYKujH+liS5guomyNYsEZBoWw4YEEhR7xEhkTY25rEPQlTERI7qA6SbK1lgSkIrDYkhgkiSZYpg9YiSiOm6/QBZgIkZ0gaYeMd5ESTy1HBZDAlNr9dDpjYvNMIZJPKYeXTYkUAswESO6gHGODfdgIhEZ59d4sUeMBGTsDQPYq0ti4srLZAkmYkQXaJpjw0oAiaeWPbokMONm5Eo3ORRymZ1LQ2S5Wm7oTBZgIkZ0AW4mSiLjHDESGeeHkehU3AePLMBEjOgC3EyURFbLzURJYCqumEgCkyTpvDlijGG6MiZiRBfg0skksjr2iJHA2JtAIlNr9ZAMa80whqlFmIgRXaBpfDdvoiSeGs4RI4HVqNmbQOIyxi8AeLkzhunKmIgRXaC63nAj9fNkIkbiqaozxq+7nUtCZLmm+y/jl8RTVacBAPgp3SDnYjPUAkzEiC7ARIxEVl3fWBFg/JKATBVZxi8JiPUHshQTMaILNFVk2SJL4jFWBPy9GL8kHlP88v5LAuL9lyzFRIzoPJIkmcZ4s0WLRFTFHjESmLEhzJ/xSwLi/ZcsxUSM6DyqBh30jSse8UZKIuLQGBIZK7IkMo6oIUsxESM6j7ESq5DLuOIRCampR4EVARIPh3aRyJqG1rIhgVqGiRjReWrUhkqsr9INMhlXPCLxGFdNZCJGImKPLomsabEZ3n+pZZiIEZ2nipUAEphGp0edxrAhLmOYRGQa2qVkRZbEwzoEWcqpErG4uDjIZLJmXzNmzLjo+cuXL292rqenZxuXmhxJDfewIYEZ4xcAfFkRIAFVcWgiCYxDa8lSTvWXeteuXdDpdKbHhw8fxujRozFx4sRLvsbf3x+ZmZmmxxyO5tpMw2KUTvWrQS7CuNCBl7sC7gqnamcjF8F98EhkXGyGLOVUkRIaGmr2eOHChejUqROGDRt2ydfIZDJERETYumgkCOMcMd5ESURNrbGMXxITh3aRyLhqIlnKae90DQ0N+PLLLzF79uzL9nLV1NQgNjYWer0eycnJePXVV9GjR4/LvrdarYZarTY9rqqqAgBoNBpoNBrrfAOtZPx8e5dDVBW1hp+rt4f8otdQkvQAAJ1Ox2tsA4zfq1NeUw/AsNjMxa6hTm+IX71ez2tsA4zfq6PW6NCgNcSol6L5ddRpDUmaJEm8xjbCGL46xsU6vN2aX0NNY2wDgEaj5TW2AUeK35aWwWkTsZUrV6KiogLTpk275DkJCQn49NNPkZSUhMrKSixatAiDBg1CRkYG2rdvf8nXLViwAPPmzWt2fPXq1fD29rZG8a9aenq6vYsgpD25cgByVBSdQ1ra2WbPFxcbns/IOIK0kow2L5+rYPy2zoFSGQAFdHU1SEtLa/Z8Vp7h+by8PKSlnWnz8rkKxm/rVDUAgBtkkLBpXTrkF7ShZlUanq+trb1ofJP1MIZbp7BMAUCGw/t2oTbL/DlDHmaodq9fvx4cuGA7jhC/KpWqRec5bRgsXboUN954I6Kioi55TkpKClJSUkyPBw0ahG7duuHDDz/EK6+8csnXzZ07F7NnzzY9rqqqQkxMDMaMGQN/f3/rfAOtpNFokJ6ejtGjR8PdnV3jltrx2xEg7yx6JXZG6sjOzZ5fWboHKC9Fjx7dkdo/1g4ldG6M36tTtzcPOJ6BmIgQpKb2bfZ8zvosIPckoqOjkZra0w4ldG6M36uTU1IL7NkCH6U7bho3ttnzW08UAUf2w8fHB6mpg+1QQufHGL46L+1fD0CDMSOGokuYr9lzDVo95uxYAwAYMWIEgvy87FBC5+ZI8WscLXclTpmInT59GmvWrMFPP/1k0evc3d3Rp08fZGVlXfY8pVIJpVJ50dfb+wdv5EhlEUllvWGxl2A/z4teP5nMsACCQqHg9bUhxm/rqDQSACDA2+Oi108hN8SvXC7n9bUhxm/rNG6BB39Pt4vHr5uhyiKTyXh9bYwxbDlJklCtNgRxkK9Xs+snyZqGJrq7XzzGyTocIX5b+vlOuazWsmXLEBYWhnHjxln0Op1Oh0OHDiEyMtJGJSNHV6kyjOkN9OYNksRTxYniJDBj/HLpbxKRqkEHnd7QGMYFk6ilnC4R0+v1WLZsGaZOnQo3N/NfhClTpmDu3Lmmx/Pnz8fq1atx8uRJ7N27F/feey9Onz6Nhx56qK2LTQ6iXNUAAAj09rBzSYgsZ1o1kSvOkYCquWIiCcwYvwq5DF7uCjuXhkThdHe7NWvWIDc3Fw888ECz53JzcyGXN+We5eXlmD59OgoKCtCuXTv07dsXW7duRffu3duyyORAKow9YmyRJQEZ45c9CiQiY/wGMH5JQBV1hobcAC937klLLeZ0idiYMWMgSdJFn9uwYYPZ48WLF2Px4sVtUCoSRWXj0rPt2CNGAjL26Ab7MH5JPMb4DWL8koDKahm/ZDmnG5pI1FoNWj1qGifaco4YichYEWjHigAJiPFLIiuvNTTkBrEhlyzARIyokbE3TCbjYgckJvYokMjKjT0KrMiSgMpUxoYE1h+o5ZiIETWqUDWN71ZcuJMokQBMPQqsyJKAmiqyjF8STzmHJlIrMBEjalRRx4U6SFwand60ahcrAiQi9oiRyNgQRq3BRIyokWnFRN5ESUDGYYlyGVedIzGxR4xExqHh1BpMxIgaNe0hxkosicc4UTzQ24NDa0lIpsUOWJElAbFHjFqDiRhRo0oVl64ncTVVAtiQQOJRa3WmVWs5NJFExB4xag0mYkSNys5brININKwEkMiMQ8MVchn8PJ1ui1NyAcYeXQ6tJUswESNqVFKtBgCE+intXBIiy3FYDIns/PiVc2gtCcgYw8FMxMgCTMSIGhXXMBEjcXHpZBJZU/xyRAKJp65BhzqNDgB7xMgyTMSIGhUbe8R8mYiReMo4NJEEZloxkT26JCDj0HAPhRw+Hgo7l4ZEwkSMqFEJe8RIYKU1TMRIXIxfEpkxftv5uEMm49BaajkmYkQA9HoJJY030hD2iJGAiqrrAbAhgcRkjN8wxi8JqCl+Pe1cEhINEzEiABV1Guj0EgAg2JctsiSeosahtawIkIiKqhrj15/xS+Jpuv+yIYEsw0SMCE3zw4J8POCu4K8FicdYkQ33Z0WAxFPIiiwJrLCqsUeMDQlkIdY4idCUiIWwN4wEVKvWmjbDZUWARFTEiiwJjD1i1FpMxIjAhTpIbMZKgI+HAr5KboZL4jHGMHt0SUTGhoRwNiSQhZiIEYFL15PY2JtAImvQ6k2b4XKOI4mIPWLUWkzEiNC0mTNXTCQRcX4Nicx4/3VXyNDOmxs6k3gK2SNGrcREjAhAfqXhJhoRwJsoiYc9YiQy00IHfp7cg4mEo9NLplE1YRxaSxZiIkYE4FxFHQAgKtDLziUhspxpfg17xEhATUvXM35JPKW1auglQC4DgrkhOVmIiRgRgPzGRCySPWIkoKalk1mRJfFwM2cSmbEhIdhXCTduf0MWYsSQy9Pq9ChorMhGs0eMBGQcWsv5CSQi09Bwxi8JiPFLV4OJGLm8omrDsAJ3hYyLdZCQ8soNPbrt23nbuSREljvL+CWBnS1XAQDat2NDLlmOiRi5POP8sHB/T8jlnChOYjm/RzeGFQESUB4rsiSwpoYwxi9ZjokYubxzjcMKuFAHiSi/sh46vQQPhZw9uiQk9oiRyIzxy6kN1BpMxMjlmVZM5EIdJCBTJaCdF3t0STj1Gp1p1c9o9iiQgM5WGHt02ZBAlmMiRi4vn0vXk8DyKtgaS+IyLnTg7aHgZs4kpLzzGsOILOVUidjLL78MmUxm9pWYmHjZ13z//fdITEyEp6cnevbsibS0tDYqLTmK3DJDaxZvoiQiThQnkZ0fv9zMmURTq9aiXKUBwDoEtY5TJWIA0KNHD+Tn55u+Nm/efMlzt27dismTJ+PBBx/Evn37MGHCBEyYMAGHDx9uwxKTvZ0uNVQEOgb72LkkRJbL4/wEEhjjl0RmHJHg7+kGf0/26JLlnC4Rc3NzQ0REhOkrJCTkkue+/fbbuOGGG/D000+jW7dueOWVV5CcnIz33nuvDUtM9qTV6XGmsUU2NoSJGInHGL9sjSURMX5JZGdMI2o4P4xax83eBbC2EydOICoqCp6enkhJScGCBQvQoUOHi567bds2zJ492+zY2LFjsXLlyst+hlqthlqtNj2uqqoCAGg0Gmg0mqv7Bq6S8fPtXQ5RnClXQaOT4OEmR4iX4orXTZL0AACdTsdrbAOMX8vllNQCANoHKq943XR6Q/zq9XpeYxtg/FruZFENAKBDO68rx69WCwCQJInX2EYYw5bJLqoGAMQGXTl+NVp90/81Wl5jG3Ck+G1pGZwqERswYACWL1+OhIQE5OfnY968eRgyZAgOHz4MPz+/ZucXFBQgPDzc7Fh4eDgKCgou+zkLFizAvHnzmh1fvXo1vL0do1UkPT3d3kUQwrEKGQAF2rnrsGrVn1c8v7hYDkCOjIwjSCvJsHn5XBXjt2XUOqCwynAbz9q7FfmHLn9+Vp4h3vPy8pCWdsb2BXRRjN+WO5CjACBDSc4RpFVc/p6aVQkAbqitreV8bhtjDLfMxpOGOoGuIh9paXmXPdeQhxnu1+vXr4eXU9XAHYsjxK9KpWrReU4VBjfeeKPp/0lJSRgwYABiY2Px3Xff4cEHH7Ta58ydO9esJ62qqgoxMTEYM2YM/P39rfY5raHRaJCeno7Ro0fD3Z3jla+kfOcZ4OhR9OgQhtTUPlc8f2XpHqC8FD16dEdq/9g2KKFrYfxa5kh+FbBzO9p5u2PiLWOueH7O+iwg9ySio6ORmtqzDUroWhi/ltHrJTyzey0APe68YRhigy/fkLn1RBFwZD98fHyQmjq4bQrpYhjDlvlm2W4AZRg9IAmpfaIue26DVo85O9YAAEaMGIEgPw7HtTZHil/jaLkrcapE7EKBgYHo2rUrsrKyLvp8REQECgsLzY4VFhYiIiLisu+rVCqhVDbfONXd3d3uP3gjRyqLIztbblg6uWOob4uul0xmmFapUCh4fW2I8dsyueWGIdLxLYxfhdwQv3K5nNfXhhi/LXOuog71Gj3c5DLEhfrBTXH5aesKN0OVRSaT8fraGGO4ZXJKDL0enSP8r3i9JFnT0ER3dzdeXxtyhPht6ec73WId56upqUF2djYiIyMv+nxKSgrWrl1rdiw9PR0pKSltUTxyAKdKDfNr4q7QEkvkiIzzwzpyoRkSkDF+OwR7XzEJI3I0tWotCqoMjbnxvAdTKznVne+pp57Cxo0bcerUKWzduhW33norFAoFJk+eDACYMmUK5s6dazr/iSeewKpVq/DGG2/g2LFjePnll7F7927MnDnTXt8CtbHjhYaJ4p3Dms8hJHJ0xopsfCgrASSek8WG+298iK+dS0JkOWNDbpCPBwK9PexcGhKVUw1NPHv2LCZPnozS0lKEhoZi8ODB2L59O0JDQwEAubm5kMubcs9BgwZhxYoVeP755/Hcc8+hS5cuWLlyJa655hp7fQvUhuoadKalk7uGsyJA4jFWZLkHHokou9jYo8sRCSSek8UcUUNXz6kSsW+++eayz2/YsKHZsYkTJ2LixIk2KhE5sqyiGkgSEOzjgWDf5nP+iByZXi+ZenS7hLNHl8RzvNCw9Dfjl0RkjN+ujF+6Ck41NJHIEk2VAPaGkXhyy1So0+igdJOzRZaEI0kSjhUY7sHdIuy72jBRaxzNN8RvYgQTMWo9JmLkso4XsTWLxHWswLA0bpdwXy50QMIprlGjrLYBchkbw0hMxntwAhsS6Crwrze5rBMc1kUCM/YmJLISQAI61tibEBfiA093hZ1LQ2SZ6noNzpbXAWCPGF0dJmLksg7nVQIAukfyJkriOcZhMSSwzALGL4nLOLUh3F+Jdj5cMZFaj4kYuaSiqnoUVashlwHdItmjQOI5ahoWw4osiedofmP8hvP+S+I50tgQxmGJdLWYiJFLOnzO0BvWKdQX3h5OtXgouYDy2gacLjVsvdAzOsDOpSGy3P6zFQCApPaMXxLPgTMVAIAk3n/pKjERI5d06KyhNZaVWBLRgcZKbMcQH24kSsKpVGlMezD1igm0b2GIWsGYiPVm/NJVYiJGLulQ4/ywa5iIkYD2sxJAAjuYVwEAiA32RhDn15Bgqus1yCo2LPbFhgS6WkzEyOVIkoT9Z8oBcFgMicmYiPVi/JKA9udWAAB6tQ+0azmIWuPg2UpIEhAd6IVQP6W9i0OCYyJGLudUqQolNQ3wcJOjJyuyJBhJkpqGxXRoZ9/CELUCe3RJZKb47RBo13KQc2AiRi5nV04ZAKB3+0Ao3bh/DYklt0yFcpUGHgo5unHrBRKMYURCBQBWZElMxvjtw4YEsgImYuRydp0yJGL94tibQOLZdcowrLZ7lD8bEkg4OSW1KK1tgIdCju7cOoQEo9dL2HPacA9mjy5ZAxMxcjnGROzajkF2LgmR5bZmlwAABnUKtnNJiCy3NbsUAJAcGwhPdzYkkFgyC6tRVtsAbw8FkjjHkayAiRi5lKLqepwqVUEmA5I5v4YEI0kStjVWZAd1CrFzaYgsx/glkRkbEq6NC4KHG6vQdPUYReRStmYZbqLdIvwR4OVu59IQWSanpBb5lfXwUMg5tJaEo9dLph7d6zqzR5fEszWL8UvWxUSMXMqGzCIAwLCEUDuXhMhyHNZFIjtWUI1ylYbDukhIWp0eOxoX+2KPLlkLEzFyGXq9hE0nDK1Zw7syESPxbD5hnB/GSgCJZ0tjb0L/jkFwV7D6QWI5cLYCNWotArzc0Y0LzZCV8E5ILuNQXiXKahvgp3RDciyHdZFY6jU6bDpRDAAYxoYEElD60UIAjF8SU/oRw4iaIV1CoJDL7FwachZMxMhlbMg0VGKv6xzC1lgSzrbsUqgadAj3V6JnNDciJ7GU1TZgd+OKtaO7h9u5NESWSz9SAIDxS9bF2ii5jL8yDDfREYlsjSXxrD5i6E0Y1S0ccrbGkmDWHSuCXgK6RfqjfTtvexeHyCIni2uQXVwLN7kMwxPC7F0cciJMxMgl5JTU4kh+FRRyGcZ0j7B3cYgsotdLWNs4rIutsSQi9iaQyNIbG8IGxgdzxWWyKiZi5BLSDuUDMGyC287Hw86lIbLM7tPlKKpWw1fphhRu5EyCqa7XmIaGj2EiRgIy1iHG9GD8knUxESOX8MdBw030pqRIO5eEyHI/7zsLALjxmggo3bhsPYll1eECqLV6xIf6oEcUV5sjsWQX1+DA2Uoo5DKk9mQdgqyLiRg5vWMFVTiSXwU3DkskAdVrdPi9sSHh1uRoO5eGyHI/78sDANzWJxoyGec3klhWNsbvsK6hCPFV2rk05GyYiJHT+3bXGQCGRQ44LJFEs/5YEarrtYgM8MTAjhyWSGLJr6zDtpOGjchv6c2GBBKLXi9h5X5DIjahD+OXrI+JGDk1tVZnao2d1D/GzqUhstyKnbkADJVYrpZIovl21xlIEtA/LggxQVwtkcSyOasEZ8rq4Kd0w+hunB9G1sdEjJzaqsMFqFBpEBngiaFduGw9ieVkcQ3+PlECmQy4u38HexeHyCIanR4rdhgaEu4ZyPgl8Xy+7TQA4Pa+7eHlwfm5ZH1OlYgtWLAA1157Lfz8/BAWFoYJEyYgMzPzsq9Zvnw5ZDKZ2Zenp2cblZhsSZIkfLrlFADgzn4xULA3gQTz5XZDJXZ411B0CGZvAolldUYhiqrVCPFV4sZruMgBieVsuQrrjhmWrb+XDQlkI06ViG3cuBEzZszA9u3bkZ6eDo1GgzFjxqC2tvayr/P390d+fr7p6/Tp021UYrKlPafLceBMBTzc5LgvJdbexSGySGWdBt/vNsxvnJISZ9/CEFlIkiQs3XwSADC5fww83JyqukEuYPmWU9BLhm1vOof52bs45KTc7F0Aa1q1apXZ4+XLlyMsLAx79uzB0KFDL/k6mUyGiAiupudsPtpkqATc1ieaKx2RcL7YdgrVai26hPliWFcOqyWxbDtZir25jQ1hA9kQRmIpq23AV43DaqcPibdzaciZOVUidqHKykoAQFBQ0GXPq6mpQWxsLPR6PZKTk/Hqq6+iR48elzxfrVZDrVabHldVVQEANBoNNBqNFUreesbPt3c57C3jXBVWHymETAZMGRhjteshSXoAgE6nc/lrbAuMX4NatRZLN+cAAB4d2hE6nRY63dW/r05viF+9Xu/y19gWGL9N3l17AgAwqW802nkprHJNdFotAENvG6+xbTCGDT7elIU6jQ49ovxwXXygVa6HRqtv+r9G6/LX2BYcKX5bWgaZJEmSjctiF3q9HjfffDMqKiqwefPmS563bds2nDhxAklJSaisrMSiRYuwadMmZGRkoH379hd9zcsvv4x58+Y1O75ixQp4e3MehyP48KgcRyrk6Buix5Qu+iu/oIU+PibH4XI57orXISXcKX91yAGsOyfDL6cVCPGU8FxvHRRWmt6YnifD77kKDAzTY3In6/1eEJ0vpxp467AbFDIJz/fRIchKAxKyKoF3j7gh3Mvwe0FkCyotMG+vAvU6GR7oqkOvYOv8rdfqgTk7DP0fC6/Vwsupu0JIpVLh7rvvRmVlJfz9L72RvdOGwYwZM3D48OHLJmEAkJKSgpSUFNPjQYMGoVu3bvjwww/xyiuvXPQ1c+fOxezZs02Pq6qqEBMTgzFjxlz2YrcFjUaD9PR0jB49Gu7u7nYti73sy63AkW07oZDLsOCeIegY4mO1915ZugcoL0WPHt2R2p/DbayN8QtUqDR48a2/AWgx+4ZrML6v9fauyVmfBeSeRHR0NFJTe1rtfcmA8Wvorbrrk10AKnB7cnvcO+HSo0sstfVEEXBkP3x8fJCaOthq70tNGMPAgj8zUa87ja5hvnjmnhSrbRvSoNVjzo41AIARI0YgyM/LKu9LTRwpfo2j5a7EKROxmTNn4vfff8emTZsu2at1Ke7u7ujTpw+ysrIueY5SqYRS2byJz93d3e4/eCNHKktb0uslLPzrOADg9uRodI0MtOr7y2SGCecKhcIlr29bcdX4BYD3Nx5HZZ0WiRF+mNQ/1qqrfSrkhviVy+Uue33bgivH728HzmFvbgW83BWYMzbRqtdB4WaosshkMpe9vm3FVWM4p6QWXzTODXtuXDcolR5We29J1jQKwd3dzSWvb1txhPht6ec71TJGkiRh5syZ+Pnnn7Fu3Tp07NjR4vfQ6XQ4dOgQIiO51K6Iftx7FntzK+DtocCTo7vauzhEFskursGX2w2rtj4/rju3XCCh1Gt0WPjnMQDAo8M6IdyfW8GQWBakHYVGJ2FY11AMTwizd3HIBThVj9iMGTOwYsUK/PLLL/Dz80NBQQEAICAgAF5ehi7gKVOmIDo6GgsWLAAAzJ8/HwMHDkTnzp1RUVGB//73vzh9+jQeeughu30f1DqVdRpTJeCJkV0QGcBufxKHXi9h7k+HoNVLGJkYhsFdQuxdJCKLvLXmBPIq6hAZ4ImHh3KlORLL6owCrD5SCIVchufHdbN3cchFOFUi9sEHHwAAhg8fbnZ82bJlmDZtGgAgNzcXcnlTR2B5eTmmT5+OgoICtGvXDn379sXWrVvRvXv3tio2Wcn8346gtLYBnUJ9cP91lveGEtnTVztzsTOnDN4eCrx8s/Xm1RC1hUNnK/Hx34YtQ+bfcg28PBR2LhFRy1XWafD8ysMAgIeHxqNLOPcNo7bRqkSsoqICP/zwA7Kzs/H0008jKCgIe/fuRXh4OKKjrTex3FItWQByw4YNZo8XL16MxYsX26hE1Fb+yijAj3vPQiYDXrs9iZuHklDOlKmwMO0oAODpsQmICeLqqySOeo0OT/9wADq9hJuSIjG6e7i9i0Rkkfm/HUFRtRrxIT54YmQXexeHXIjFidjBgwcxatQoBAQE4NSpU5g+fTqCgoLw008/ITc3F59//rktykl0SYVV9Xjup0MAgEeGdkK/uMvvG0fkSBq0esxcsRe1DTr0i22HKSlx9i4SkUX+/ccRHCuoRrCPB3tzSTg/7jmLH/eehVwGLLw9CZ7u7M2ltmNxt8Hs2bMxbdo0nDhxAp6eTRNxU1NTsWnTJqsWjuhKGrR6PPblHpTWNiAxwg9PjmZLFollwZ9HceBsJQK83PH25D5coIOE8tuBc/hyey5kMmDxpN4I8bXSpmFEbSCrqNo0JPGJkV3RvyMbcqltWZyI7dq1C4888kiz49HR0abFMYjayr//OIK9uRXw83TDknv7QunGliwSxy/787BsyykAwBsTeyE6kAvMkDiO5lfh2R8PAgD+MbwThnYNtXOJiFqusk6DR7/cizqNDoM6BWPm9Z3tXSRyQRYnYkql8qKblB0/fhyhobwJU9v5bOspfL7NsNT3W5N6I86KGzcT2drOnDI8/b2hEvvIsHiM4rwaEkhhVT0eWL4LtQ06pMQH48lR3C6ExGEcTZNVVIMIf0+8dVdvjkYgu7A4Ebv55psxf/58aDQaAIaNFXNzc/HMM8/g9ttvt3oBiS7m94Pn8PJvGQCAp8Z0xchurMSSOLKKavDwF7vRoNPjhh4ReGZsor2LRNRiVfUaPLB8F/Ir69Ep1AdL7u0LNwUXSCIxGLcK2ZpdCh8PBT6ddi3C/LjnHdmHxXfON954AzU1NQgLC0NdXR2GDRuGzp07w8/PD//5z39sUUYiM+szizD72wOQJOC+gbGYMYLDCUgc2cU1mPzxdlSoNOgdE4jFk3pDzpZYEkR1vQZTP92JjHNVCPbxwLJp/RHg7W7vYhG1iF4v4f9WHsKPe89CIZfhvXuS0T3K397FIhdm8aqJAQEBSE9Px+bNm3Hw4EHU1NQgOTkZo0aNskX5iMz8lVGAmSv2QqOTkNozAi/f3AMyGSuxJIbs4hpM/mg7iqvVSIzww6fTruV+SySMGrUWUz/diX25FQjwcsdnD/RHh2ButUBi0OslPP/LYXy98wzkMuDNO3thREKYvYtFLq7VGzoPHjwYgwcPtmZZiC7rl/15mP2dYa+acT0jOaabhHLwbAUeWL4bJTWGJGzF9IEI8vGwd7GIWqS4Wo0HP9uFg40rfH710ABcEx1g72IRtUiDVo9nfjyIn/flNSZhvXFLb/vte0tk1KJE7J133mnxG/7zn/9sdWGILkaSJLy/PguLVh8HANzWJxqv35HEOQkkjHXHCjHjq32o0+jQLdIfXz7Yn0kYCSOrqAbTlu3E2fI6BPl44LP7+zMJI2FU1Wvw6Bd7sDW7FAq5DIsmJmFCHyZh5BhalIgtXrzY7HFxcTFUKhUCAwMBABUVFfD29kZYWBgTMbKqeo0Oz/10CD/tywMA3H9dHJ4f1509YSQESZLw8d8nsfDPY9BLwJAuIfjfPcnw8+ScGhLDumOFePLbA6is0yAu2BvL7+/PFWpJGFlF1Xj0y73IKqqBj4cC/7u3L4ZxmwVyIC1KxHJyckz/X7FiBf73v/9h6dKlSEhIAABkZmZi+vTpF91fjKi1sotrMOOrvThWUA2FXIZ5N/fAvQNj7V0sohaprNPg6e8PYPWRQgDAnf3a4z+39oQ7e3JJADq9hMXpx/He+iwAQJ8OgfhkSj8Ec8NmEsSvB87h2R8PQtWgQ7i/Ep9OuxY9otiTS47F4jliL7zwAn744QdTEgYACQkJWLx4Me644w7cc889Vi0guR5JkvDzvjw8v/IwVA06hPh64K1JfTC4S4i9i0bUIrtOlWHOdweQW6aCh0KOF8d3xz0DOnBhGRLCmTIVnv7hALafLAMATEmJxf+N6walGxeWIcdXo9biP38cxdc7cwEAgzoF4+27+iDUj40I5HgsTsTy8/Oh1WqbHdfpdCgsLLRKoch1FVTW4/mVh7HmqCGWUuKD8fZdvRHmzz0+yPHVa3RY9Fcmlm7JgSQB0YFe+ODeZCS1D7R30YiuSJIkrNiZi1f/OIraBh28PRRYcFtPLmpAwtiaVYKnfziIvIo6AMCMEZ0we3QCpzOQw7I4ERs5ciQeeeQRfPLJJ0hOTgYA7NmzB4899hiXsKdW0+slfLf7DP7zx1FUq7VwV8jwz+u74B8jOvMGSkLYklWCF1YexsmSWgCGoYjP39Qd/pwPRgLILq7BS79kYHNWCQDg2rh2+O8dvTgfjIRQXtuA/67OxIodhl6w9u288PodSRjUiSNpyLFZnIh9+umnmDp1Kvr16wd3d0MFQ6vVYuzYsfjkk0+sXkByfrtPlWH+70dw8GwlAKBXTCD+e0cSuob72blkRFd2pkyF//xxFKsyCgAA4f5KLLwtCSMSuT8NOb7qeg3eXZeFTzfnQKuXoHST4183JGLaoDg2gpHD0+klrNhxGotWH0dlnQYAcM+ADpib2g2+ylbv0ETUZiyO0tDQUKSlpeH48eM4duwYACAxMRFdu3a1euHIuZ0pU+H1vzLx24FzAABfpRtmjeqC+6/ryAoAObzKOg0+3nQSH/99EmqtHgq5DPcNjMWTo7oiwJu9YOTYNDo9fthzFm+mH0dxtRoAMDIxDC/c1J29YOTwJEnCxuPFeG1VJo7mVwEAEiP88NL4HkjpFGzn0hG1XKubC7p27crki1olr6IO763Lwg97zkCjkyCTAZP6xWDOmAROpiWHV6PWYtnmHHz890lU1Rvmy6bEB+Olm7sjMcLfzqUjujydXsLKfXl4e+0J5JapAAAdQ3zw4k3d2YtLQtiWXYo3Vmdi9+lyAECAlzvmjOmKu/t34P6iJByLE7EHHnjgss9/+umnrS4MObczZSos2ZiN73YbEjAAGNw5BHNTE7mkLDm8ClUDvtqRi0/+PolylWEITEK4H54c3RVje4RzRURyaGqtDr/uP4clG7ORXWyYxxjiq8Q/hnfCvQNj4eHGCiw5LkmS8PeJEizZmI2t2aUAAKWbHFNSYvHY8M4I8vGwcwmJWsfiRKy8vNzssUajweHDh1FRUYHrr7/eagUj5yBJEvbmluOTv3PwV0YB9Ib8C4M6BWPWqK7o3zHIvgUkuoLcUhWWbj6J73afRZ1GBwCID/HBrNFdcVPPSMg5jJYcWKVKgy93nMZnW0+hqHEIYqC3Ox4d1glTUmLh7cF5NOS4GrR6/HrgHD75+ySOFVQDANwVMtzdvwP+MaIzwrmiMgnO4jvwzz//3OyYXq/HY489hk6dOlmlUCQ+tVaHvzIKsXRzDg6cqTAdH9IlBDNHdMaAeI7hJsel10vYdrIUX24/bdaA0C3SHw8P7YjxSVEcAkMO7VhBFb7ekYvv95yFqsHQgBDh74lp18XhngEd4MfVPMmBFVbV47tdZ/DljtMorDI0IHh7KDDp2hg8NCQe0YFedi4hkXVYpSlMLpdj9uzZGD58OP71r39Z4y1JUMcLq/HtrjP4ae9Z0/AtDzc5bu0djQcGd0RCBFdCJMdVVF2PH/acxbe7zuB0qcp0fFjXUDw8NB6DOgVzCCI5LFWDFr8fzMfXO3OxL7fCdLxbpD+mD+mIm5KiOASRHJZOL2HT8WKs2JmLdceKoGtsAQv3V2LaoI64u38HLoRETsdqYxKys7MvutEzOb/qeg3+PFSAb3blYu95f/wj/D1xV/8Y3DswFiG+XISDHFODVo9Nx4vx496zSD9SCG3jH38/pRtu6ROF+wbGsQGBHJZh+HcFft53Fr/sO4dqteHvsJtchtHdw3HPgFhc15kNCOS4ckpqsXJfHr7ffQbnKutNx6+Na4fJ/TuwAYGcmsWJ2OzZs80eS5KE/Px8/PHHH5g6darVCkaOrV6jw7pjRfh1/zmsyyxCg1YPAFDIZRiZGIa7+sdgaJdQDt8ih6TXS9h5qgy/7D+HtEP5pv1nAKBPh8DGP/6RnD9DDutEYTVW7s/DrwfO4UxZnel4bLA3Jl0bgzv6tkeYH+fPkGMqqqrHbwfz8ev+PBxo3EMUMMxfvK1Pe0zuH4Mu3EuUXIDFtYx9+/aZPZbL5QgNDcUbb7xxxRUVSWwanR6bT5Tg1wPnsDqjALWN8w4AID7UBxP7xuD2vtH8408OSa+XcDCvEn8eysevB84h/7yW11A/JcYnReHOa9tzCXpyWKdKavHn4QL8euCcae8kwDB3ZmyPCNye3B6DOgVzARlySKU1aqw5WohfD5zDtuxS09xbhVyGwZ1DcFtyNMb2iICnu8K+BSVqQxYnYuvXr7dFOchB1ai12HS8GKszCrDuWJFp3yQAiA70wk29InFzryh0j/Tn0BdyOBqdHjtOlmH1kQKszihEQVVT8uXn6YYbr4nALb2jMTA+mJuIk8ORJAkZ56qwOqMAf2UUIrOw2vScu0KGYV3DcEvvKIzqFg4vD1ZeyfGcKVNh9ZFC/JVRgN2nykzJFwD0jW2HW3pHIbVnJKcvkMuyOBG7/vrr8dNPPyEwMNDseFVVFSZMmIB169ZZq2xkJ0XV9Vh7tAirMwqwJasUDTq96bkQXw+M6xmJ8b2ikNyhHVteyeHUNeiwsbHxYO2xIrNhhz4eCgxPDMP4pCgMTwhlyys5HK1Oj92ny7E6w1B5zatoGnboJpdhQHwQbkqKwo3XRCDQm3snkWORJAmZhdWm+M04V2X2fI8of6T2NDTgxgR526mURI7D4kRsw4YNaGhoaHa8vr4ef//9t1UKRW1Lrze0um48XoR1x4qw70wFpPNareKCvTGmRwRGdw9Hcod27DkghyJJErKLa7AhsxgbjxdjR06Zac4iAAT7eGBUt3CMvSYcgzqFMPkih1NYVY+Nxw3xu/lEiVnjgae7HMO6hmJsjwiMTAznqnHkcKrrNdiSVYqNx4ux6XixWeOBXAZcGxeEsY11CCZfROZanIgdPHjQ9P8jR46goKDA9Fin02HVqlWIjo62bunIZkpr1NicVYKNmcXYdKIYJTXmyXWvmECM6R6OMd3D0TnMl8MOyaFc7g8/ALRv54Ux3SMwtkc4+sUFsfGAHEqDVo/dp8sMyVdmsWmjWqMAL3eM7BaGsT0iMLRLKIcdkkPR6yUcya8yNR7sPV1uWm0WMGxZM6RziKHxoFsYgjnskOiSWpyI9e7dGzKZDDKZDNdff32z5728vPDuu+9atXCt9f777+O///0vCgoK0KtXL7z77rvo37//Jc///vvv8cILL+DUqVPo0qULXnvtNaSmprZhiW1Pq9PjwNkKbGzsNTiYV2nW6+XjocB1nUMwLCEUo7qFc7d6cigt+cM/oGMQhnUNxfCEUHQKZeMBOZYzZSpsaEy8tmWXmC12JJMBSe0DMaxrKIZ1DUWv9gFccZYcSlltA/4+UdzY+FWCkhq12fPxIT4Y2jUUwxJCMbBjMBsPiFqoxYlYTk4OJElCfHw8du7cidDQUNNzHh4eCAsLg0Jh/1+8b7/9FrNnz8aSJUswYMAAvPXWWxg7diwyMzMRFhbW7PytW7di8uTJWLBgAW666SasWLECEyZMwN69e3HNNdfY4TuwnsKqelPi9feJYrOFNgCge6Q/hiUY/vAnd2jHfTrIofAPP4msrkGH7TmlhlEHx4txsqTW7PkQXyWGdg3BsK6hGNIlFEE+nO9FjkOnl3DwdLmp8evgWfMpC94eCgzqZGi8HdYlFB2COeSQqDVanIjFxsYCAPR6/RXOtK8333wT06dPx/333w8AWLJkCf744w98+umnePbZZ5ud//bbb+OGG27A008/DQB45ZVXkJ6ejvfeew9Llixp07JfrQatYYW4X0/L8cF7W3GssMbs+QAvdwzpEmJqdQ1jrxc5EJ1ewr7cCqSdkWPph9txKK+Kf/hJGIa5irXYkC/DD5/twc5T5VCfN1dRIZehb2w70/23e6Q/Fzsih1JUVY91Rwvw7XE5Xty/HpV15o23iRF+psbbfrFBbLwlh1LXoMOWE8WobL6MhUNrUSL266+/4sYbb4S7uzt+/fXXy5578803W6VgrdHQ0IA9e/Zg7ty5pmNyuRyjRo3Ctm3bLvqabdu2NdukeuzYsVi5cuUlP0etVkOtbmqdr6oyrAqk0Wig0Wgu9TKb25JVivs/2wNADqDGMNwlOgBDuwRjSJcQJEUHmM2VsWdZRSVJhoqVTqfj9bOComo1Np0owd8nSrAlu7TxD78cgOF3KjHcF0O6hGBolxAkdwg0+8PP6285XWNDml6v5/Wzgup6LbafLMOmLEMM51XUA1AAKAUARAZ4YmiXEAztEoyU+CD4eTYttKHTaaHTXfx96eJ0WkNiIEkS49cKGrR67DtTgU0nSrDpROl5cxXlALTw93TD4M6G+sPgzsGIOL/xVtJBo2EAW0JzXsOMRqNlDF8lSZJwrKAGm7NLsPlEKXadLodGJ+GOjjLc4QDXtqU/3xYlYhMmTEBBQQHCwsIwYcKES54nk8mgs+NflpKSEuh0OoSHh5sdDw8Px7Fjxy76moKCgouef/5iJBdasGAB5s2b1+z46tWr4e1tv1Z6jR4IVioQ7y+hW6CEhAAJvu6lQH0p8g8dR/4huxXNaRQXywHIkZFxBGklGfYujnB0eiCnGjhSIcexChnyVOY9Al4KCYmBjV8BEgKVFYCuAmXHsrDm4r/CZIGsPBkABfLy8pCWdsbexRGOJAHnVMDRChmOVshxshrQS00x7CaT0NnfEL/dAiWEe9VAJquB9tQp/H3KfuV2FlmVAOCG2tpapKWl2bs4QipTN8ZvuQzHK2VQ65viVwYJHXzRGL96dPDVQiHLAwrysPfSVSJqIUMeZqh2r1+/Hl4Wr1tOKi2QWSHDkQoZjlXIUKUxr0O08zAMo0lPT7dH8cyoVKoWndeiMDh/OKKjD01sC3PnzjXrRauqqkJMTAzGjBkDf39/O5YMuHFMA9asWYPRo0fD3Z3LHFvbytI9QHkpevTojtT+sfYujhDyK+sbW1xLsDW7DDXqpuEuMhnQM9ofQzsber0Sw72wfu1axq+N5KzPAnJPIjo6GqmpPe1dHCFU1WmwJbsUm06U4u8TJSisNp+rGBvkjaFdDb1efaL9sGXjOsavjWw9UQQc2Q8fHx+kpg62d3GEoNYa9qXbdNxwD84qNp+rGOzjgSGNvV7XdQ6Gv4cM6enpjGEbaNDqMWfHGgDAiBEjEOTnZecSOT7DQl3V2NhYh9h/psJsU3Avd8NCXUO6hGBwp2C0D3B3mDqwcbTclThVPh4SEgKFQoHCwkKz44WFhYiIiLjoayIiIiw6HwCUSiWUyubLsbq7u9v9B2/kSGVxJjKZYWicQqHg9b2EBq0eu08ZlubekFmMzELzpbmDfTwMi2x0DcWQLiFmSxsbu/IZv7ahkBviVy6X8/pegiQZ91UsxobMIuzNrYDuvL/8nu5yDOoUguGNc2Vig31MzzF+bUvhZqiyyGQyXt/LOFOmwobMImw8Xoyt2aVQnbdCp1wGJHdoh+EJoRieENZsriJj2HYkWVNHhru7G6/vJZTXNmDTieJLbq/UNdwXwxPCDHMV49pB6da0UJcjxW9LP79Fidg777zT4g/+5z//2eJzrc3DwwN9+/bF2rVrTUMo9Xo91q5di5kzZ170NSkpKVi7di1mzZplOpaeno6UlJQ2KDGRc6is02BDZhFWHynExsxis14vuQzoHRNounH2jA7gIgXkUNRaHbZllyL9SCHWHC1EYZV5r1enUB8MTwjD8IRQXBsXxE3ByaHo9RIO5lUi/UgBVmcU4kSR+UJdYX7Kxq09wjC4cwg3BSeHk1NSi/QjBUg/Uog9p8vNer18lW64rnMwhnUNw7CEUEQHOldPYosSscWLF7fozWQymV0TMQCYPXs2pk6din79+qF///546623UFtba1pFccqUKYiOjsaCBQsAAE888QSGDRuGN954A+PGjcM333yD3bt346OPPrLnt0Hk8PIq6pCeUYD0o4XYcbLMbF+vEF9Dr9fwhDAM6RyCdlyamxxMpUqD9ZlFSD9SiI3HzRsPjCt0Gnu9YoK4Qic5FrVWh63GxoMjhSg6b8js+St0Dk8wrNDJfRXJkej1EvafrUD6kUKkHylE1gWNB8YVOod3DUPfWOfeXqlFiVhOTo6ty2E1kyZNQnFxMV588UUUFBSgd+/eWLVqlWlBjtzcXMjlTT/QQYMGYcWKFXj++efx3HPPoUuXLli5cqXwe4gR2cLp0lr8fjAfaYfykXHOfPxzlzBfjO4ejtHdw9GrfSB7vcjhlNU2YNXhAvxx6FyzxoMwPyVGNcbvoE7BZsNdiBxBvUaHDZlF+O1gPjYcKzLbFNzHQ4HhCWEY0yMcw7uGsdeLHI5OL2HXqTL8fvAc/sooRPF5jQduchkGxgdjdPdwjOoe7nS9XpdzVXPEpMZNfhytpWXmzJmXHIq4YcOGZscmTpyIiRMn2rhURGI6V1GHPw7m47eD53DwbKXpuFwG9IsNMiVfcSE+l3kXIvuorNNgdUYBfj+Yj81ZJWbzvbqGGxsPIpDEIbPkgBq0evx9ohi/HTiH9COFZslXhL8nRnUPw6hu4Uhh4wE5IEmSsDe3Ar8fPIc/Duab9dz6Kd0wLCEUo7uHY3hCGAK8XLPxoFWJ2NKlS7F48WKcOHECANClSxfMmjULDz30kFULR0T2UVWvwe8H8vHT3rPYfbrcdFwuA67rHIKbkiIxqlu42UIbRI6iQavHumNF+HHvWWzMLEaDrmmSfI8of9yUFIXUnhFmC20QOQpJkrDndDl+2HMWaYfyUVXfNGw2OtALNyVFIrVnJOfbksM6WVyDH/acxS/7zyGvos503N/TDWN7RGBcUiQGdQpx6iGHLWVxIvbiiy/izTffxOOPP25a0GLbtm148sknkZubi/nz51u9kERke3q9hO0nS/H9nrP483A+6jWGyqtMBvSPC8JNvaJw4zURCGHyRQ4qs6Aa3+0+g5X78lBa27TSVtdwX4xPisK4pEjEh/rasYREl1ZYVY8f957FD7vP4mRJ0zLzYX5KjEuKxE1JUUjuEOhwo5CIAKBGrUXawXx8v+cMdp1qasD18VBgTI8I3JQUiSFdQpl8XcDiROyDDz7Axx9/jMmTJ5uO3XzzzUhKSsLjjz/ORIxIMKU1anyz6wy+3pmLs+VNLVddwnxxZ78Y3Nw7CuH+nnYsIdGl1Wt0+HX/OXy14zQOnDd0NtRPiduT2+O25Gh0DfezYwmJLk2vl7DxRDG+2HYaGzKLTKvFeXsokNozErcnt0f/jkFQsOeLHNThvEp8vu0Ufj+Yb9oqQS4DhieE4Y6+7XF9YhhXmr0MixMxjUaDfv36NTvet29faLXai7yCiBxRxrlKLN9yCr8cOIcGraH3y0/phvG9o3Bnvxj0ah/AlldyWPmVdfhi22l8vTMX5SrD3jFuchlGdgvDnf1iMKxrKNwUbHklx1Sj1uKH3Wfw2bbTyDmv96tfbDvc2S8GqUmR8FU61Vav5ES0Oj1WHynE8i2nsPNUmel4xxAfTOzXHrcnt2cDbgtZ/Ft+33334YMPPsCbb75pdvyjjz7CPffcY7WCEZH1SZKEDceL8cH6bLObZ8/oAEwbFIfUnpHw8mDLFTmuo/lVeH99Fv48XGBaeCM60Av3DozFxH7tOXSWHFpRVT0+2nQS3+w6Y9oywU/phon9YnDPwA7oxKGz5MBUDVqs2JGLTzfn4FxlPQBDA9iNPSNx38BYXBvXjg24Fmr1Yh2rV6/GwIEDAQA7duxAbm4upkyZgtmzZ5vOuzBZIyL70OslpB8txHvrsnAozzB8y3jznDYojvMOyOEdOFOBd9dlYc3RQtOxgfFBmDYoDqO6hbP3ixxaXkUdlmzIxre7z5hGIMSH+uD+QXG4Lbk9fNj7RQ6sul6Dz7edxtLNOShrnH8b7OOBuwd0wD0DYhERwN6v1rL4N//w4cNITk4GAGRnZwMAQkJCEBISgsOHD5vOY6WOyP4kScKGzGK8tuoYjhVUAwC83BW4d2AHPDg4njdPcnhH86vw2qpj2JBZDMCweMy4npH4x/DO6B7lb+fSEV1eUXU93l5zAt/uOmPat65fbDvMuL4zhnUJ5aqH5NDqNTp8uiUHSzZkm1bv7BDkjX8M74QJfaI598sKLE7E1q9fb4tyEJGVHc6rxKtpR7E1uxQA4Kt0w9RBsXhwcDyCfDzsXDqiyyuorMcbqzPxw96zkCRAIZfhlt5RmDGiM4dvkcNTNWjx8aYcfLgp27SAwaBOwXj8+i4YGB/ExmpyaHq9hJ/25eGN1ZnIbxyC2CnUBzOv74zxSVEcgWBF7AsncjJltQ1YkHYU3+85CwDwUMgx7bo4zBjeGQHerrlhIomjQavHx3+fxLvrTpi2UEjtGYGnxyaiIzcNJwcnSRJ+P5iPV34/Ytq8tldMIJ67MRED4oPtXDqiKztwpgL/t/IQDudVAQCiAjwxZ0wCJvSJ5uqdNmBxIlZfX493330X69evR1FREfR6vdnze/futVrhiKjlJEnCj3vz8J8/jphWkZvQOwpzxiQgJsjbzqUjurJdp8rw3E+HcKKoBoBhCNdz47ohuUM7O5eM6MrOlKnw/MrD2HjcMIw2JsgL/xqbiJuSItkDRg6vul6DRX9l4vPtpyFJgJ+nG2aM6Ixpg+I4BNGGLE7EHnzwQaxevRp33HEH+vfvz5sLkQPIr6zD098fxOasEgBAYoQf/nNrT/SNZQWWHF9dgw6vph3FF9tPAwBCfD3w/LjuuKV3FP/GkMOTJAmfbT2FhauOoV6jh4dCjhkjOuPR4fFQurECS47v7xPFeOr7AyisMvTi3tYnGs+N68ZVaNuAxYnY77//jrS0NFx33XW2KA8RWSjtUD7m/nQIlXUaeLrLMWtUVzw4uCPcOYabBHA4rxJPfLMP2cWGvZQm94/BMzckItCb8xjJ8RVV1+Pp7w+aesFS4oPx71uv4TxGEkK9RofXV2Xi0y05AIC4YG/859aeuK5ziJ1L5josTsSio6Ph5+dni7IQkQUatHq8/FsGVuzIBQAktQ/AW5N6I54VABKAJEn4YvtpvPL7EWh0EsL8lHjjzl4Y0iXU3kUjapGt2SV4fMU+lNY2QOkmx3Op3TAlJZa9uCSE3FIVHvlyD47mG+aCTUmJxdwbu3Ev0TZmcSL2xhtv4JlnnsGSJUsQGxtrizIR0RUUVdfjsS/3Ys/pcshkwD+Gd8KsUV3ZC0ZCUGt1eHFlBr7dfQYAMLZHOBbeloR2XM2TBCBJEpZvPYV//3EUOr2EbpH+ePuu3ugazkZqEsOWrBLMWLEXFSoNgn088N+JSbg+MdzexXJJFidi/fr1Q319PeLj4+Ht7Q13d/NV2MrKyqxWOCJq7nBeJR76bDcKqurh5+mGdyf3wfCEMHsXi6hFSmrUePjz3dibWwG5DHj2xkRMHxLPXgQSgkanx//9fAjf7TasSntbn2i8eltPLmZAwli+JQevNDYi9GofgA/v68c9Re3I4kRs8uTJyMvLw6uvvorw8HD+8SRqQztOluLBz3ajRq1F5zBffDylH5f0JmHkVdThvk924GRJLfw93fDu3ckY1pVDEUkM9Rod/vHVXqw7VgS5DHgutRseHNyR9SASgiRJeDP9ON5dlwUAuD25Pf5z6zVsRLAzixOxrVu3Ytu2bejVq5ctykNEl7D2aCH+8dVeqLV6DOgYhI+n9oO/J/cFIzFkF9fgvk924FxlPaIDvfD5g/25oAEJo6peg4c+242dOWVQusnxwb3JHMpFwtDrJbz0a4ZpZdqnxnTFjBGd2YjgACxOxBITE1FXV2eLshDRJazPLMKjX+6BRidhVLdwvHd3H7ZikTBOldTiro+2o7hajU6hPvjyoQGIDPCyd7GIWkTVoMX9y3Zhz+ly+Hm64dNp1+LauCB7F4uoRSRJwgu/HMZXO3IhkwGv3HIN7h3INR4chcUz+xcuXIg5c+Zgw4YNKC0tRVVVldkXEVnXjpOlePQLQxJ2U1IkPrg3mUkYCSO/sg73fLIDxdVqJEb44ftHBzEJI2HUa3R4+PM92HO6HP6ebvh6+kAmYSSU11ZlmpKwxXf2ZhLmYCzuEbvhhhsAACNHjjQ7LkkSZDIZdDqddUpGRDheWI2HPtsNtVaPkYlhWDypN1dGJGFU1WswZelO5FXUoWOID754cACCuDIiCUKSJMz5/gA2Z5XA20OB5Q/0xzXRAfYuFlGLLd2cgyUbswEAr97aExP6RNu5RHQhixOx9evXX/K5Q4cOXVVhiKhJeW0DHvpsN6rVWvTvGIT370lmEkbC0OklPPH1PpwoqkG4vxJfPjQAoX5KexeLqMXeWZuFPw7mw10hw8dT+iG5Qzt7F4moxTYeL8Z//jgCwLA67eT+HexcIroYixOxYcOGmT2urq7G119/jU8++QR79uzBzJkzrVY4Ilel0enx2Fd7kFumQkyQF5bc25fDEUkor686hvWZxVC6yfHxlH6IDuRwRBLHqsP5WLzmOADDnJrrOofYuURELZddXIOZK/ZCLwF39muPR4bG27tIdAmtbl7ftGkTpk6disjISCxatAjXX389tm/fbs2yEbmst9ecwPaTZfDxUGDp1Gs5nIuEsvZoIT7cdBIAsGhiLyS1D7RvgYgscKZMhae/PwgAuP+6ONzFngQSSL1Ghxlf7UV1vRbXxrXDKxOu4eqIDsyiHrGCggIsX74cS5cuRVVVFe68806o1WqsXLkS3bt3t1UZiVzKjpOleH+DYZ+P1+/oha7hfnYuEVHLFVXX418/GCqxD1zXEeN7Rdm5REQtp9Xp8eS3+1Gt1qJvbDv8X2o3exeJyCL//SsTxwqqEezjgffvSYbSjaNpHFmLe8TGjx+PhIQEHDx4EG+99RbOnTuHd99915ZlI3I5lXUazP7uACQJmNi3PcYlRdq7SEQtJkkSnv7+IEprG5AY4Yd/3ZBg7yIRWeT99dnYfbocfko3vDWpN9w4L5cEsul4MZZuzgEAvH5HEsL8PO1cIrqSFveI/fnnn/jnP/+Jxx57DF26dLFlmYhc1murjiGvog5xwd54+eYe9i4OkUV+3JuHjccN88Lemcy97kgsxwur8e66EwCAf996DWKCvO1cIqKWUzVoMfcnw6J5U1JiMbIbNxwXQYubejZv3ozq6mr07dsXAwYMwHvvvYeSkhJblo3Ipew5XY4VO3IBAK/dngQfpcVr6RDZTYWqAa+mHQUAzBrVlUNqSSh6vYT/+/kQtHoJo7uH45beXOabxPL22hPIq6hDdKAXnr0x0d7FoRZqcSI2cOBAfPzxx8jPz8cjjzyCb775BlFRUdDr9UhPT0d1dbUty3lFp06dwoMPPoiOHTvCy8sLnTp1wksvvYSGhobLvm748OGQyWRmX48++mgblZrIQKvT4/mVhwEAd/RtjwHxwXYuEZFlXlt1DGW1Dega7ouHhnS0d3GILPLD3rPYdaocXu4KjkYg4WQWVGPp34YhifNv6QFvDzbkisLiwc8+Pj544IEHsHnzZhw6dAhz5szBwoULERYWhptvvtkWZWyRY8eOQa/X48MPP0RGRgYWL16MJUuW4Lnnnrvia6dPn478/HzT1+uvv94GJSZq8uPesziaX4UAL3fMZUsWCSazoBrf7DoDAPj3hJ7c746EUtegw6K/MgEAs0Z14VYLJJwFfx6FVi9hTPdwDkkUzFX9tUxISMDrr7+Os2fP4uuvv7ZWmVrlhhtuwLJlyzBmzBjEx8fj5ptvxlNPPYWffvrpiq/19vZGRESE6cvf378NSkxkUK/R4a01hnkJj1/fGcG+3PSWxPLfvzIhScCN10Sgf8cgexeHyCLLtuagqFqN9u28MO26OHsXh8gi20+WYkNmMdzkMjzHVT6FY5W+S4VCgQkTJmDChAnWeDurqaysRFDQlSsFX331Fb788ktERERg/PjxeOGFF+DtfelJumq1Gmq12vS4qqoKAKDRaKDRaK6+4FfB+Pn2LoezkiQ9AECn01ntGi/fcgr5lfWIDPDEXX2jXPpnx/i1LZ3eEL96vd5q13hvbgXWHC2EQi7DrOs7ufTPjvFrWzqtFoBhdU5rXePKOg2WbMgGADxxfSfIJT00Gr1V3ltEjGHb0Wib4kqj0VrlGkuShNf+NMzNvbNfNKIDPFz6Z+dI8dvSMjjtINKsrCy8++67WLRo0WXPu/vuuxEbG4uoqCgcPHgQzzzzDDIzMy/bk7ZgwQLMmzev2fHVq1dfNoFrS+np6fYuglMqLpYDkCMj4wjSSjKu+v0adMC7exUAZBgeUou16X9d9Xs6A8avbWTlyQAokJeXh7S0M1Z5z/ePGH4n+ofocGzXRhyzyruKjfFrG1mVAOCG2tpapKWlWeU9087IUVUvR6S3BLe8/Ug7t98q7ys6xrD1GfIwQ7V7/fr18LJCDfxouQz7zijgIZeQqDuFtLRTV/+mTsAR4lelUrXoPIdPxJ599lm89tprlz3n6NGjSExsmleTl5eHG264ARMnTsT06dMv+9qHH37Y9P+ePXsiMjISI0eORHZ2Njp16nTR18ydOxezZ882Pa6qqkJMTAzGjBlj92GNGo0G6enpGD16NNzd3e1aFme0snQPUF6KHj26I7V/7FW/34qdZ1CrPYr2gZ548b7BLr9nDePXtnLWZwG5JxEdHY3U1J5X/X5H8qtwfNt2KOQyLLhvmMvPrWH82tbWE0XAkf3w8fFBaurgq36/ugYdXn5jEwAN5o7vhRuvibj6QgqOMWw7DVo95uxYAwAYMWIEgvyu/n757bLdAMpwz8A4TL6R+zY6UvwaR8tdicMnYnPmzMG0adMue058fLzp/+fOncOIESMwaNAgfPTRRxZ/3oABAwAYetQulYgplUoolc3n8bi7u9v9B2/kSGVxJjKZIVFSKBRXfX11egnLtp4GADw0JB5enpwbZsT4tQ2F3BC/crncKtd3+TZDr9q4npGIC+XcWiPGr20o3AxVFplMZpXr+83uPJSrNIgJ8sK4Xu2hkMuu+j2dBWPY+iRZ09BEd3e3q76+h/MqsfVkGRRyGR4a2ok/r/M4Qvy29PMdPhELDQ1FaGhoi87Ny8vDiBEj0LdvXyxbtgxyueW9C/v37wcAREZGWvxaIkukHynEqVIVArzcMbFfjL2LQ2SRcxV1+O3AOQDA9CHxVzibyLHo9BI+2WxY7vuhwfFMwkg4H/99EgBwU1Kky49GEJnTjIPKy8vD8OHD0aFDByxatAjFxcUoKChAQUGB2TmJiYnYuXMnACA7OxuvvPIK9uzZg1OnTuHXX3/FlClTMHToUCQlJdnrWyEX8dnWUwCAewd24ObNJJyvdpyGVi9hYHwQerYPsHdxiCyy8XgRTpsawtrbuzhEFimqrscfB/MBsCFMdE5T+0tPT0dWVhaysrLQvr35TVWSJACGsaOZmZmmCXQeHh5Ys2YN3nrrLdTW1iImJga33347nn/++TYvP7mW06W12HayFDIZcM+Aq59rRtSWtDo9vt99FgAwJSXOvoUhaoVvdhqG1U7s256b35JwftqbB61eQp8Ogbgmmg1hInOau8+0adOuOJcsLi7OlJQBQExMDDZu3GjjkhE1991uQyVgaJdQRHFIAQlm4/FiFFWrEeTjgVHcPJQEU1Rdj3XHigAAk67lsHASiyRJ+G6XoQ4xidMahOc0QxOJRHF+b8JdrASQgL5prATc1icaHm78M0JiMfYmJHcIRJdwP3sXh8giu06V42RJLbw9FLipV5S9i0NXiX9BidrYzpwyFFWr0c7bHSPZm0CCqazTYEOmoTeBi8yQiH7db1hkhvFLIvr1QB4AILVnJHw5v1x4TMSI2tjvhwwTbMf2iGBvAgkn/UghNDoJXcN9kRDB3gQSS05JLY7kV0Ehl+GGHtw3jMSi00tYddiwCN149oY5BdYCidqQVqfHX4030XFJ3CKBxJPW2JCQ2pPxS+Ixxu91nUPQzsfDzqUhssyOnFKU1DQg0NsdgzoF27s4ZAVMxIja0I6cMpTWNqCdtztS4nkTJbFU1mnw94liAIZNnIlE83vjkt/jerI3jMRjXLJ+bPcIuCtYhXcG/CkStaH0I4UAgDHdI+DGmygJZkNmETQ6CV3CfLnIAQnnbLkKRxuHJY7pzkSMxCJJEtYcNdQhbmRDgtNgTZCoDW08buhNuL5bmJ1LQmQ5xi+JzBi/yR0COSyRhHOsoBqFVWp4uSswkCNqnAYTMaI2crq0FjkltXCTyzi2m4Sj10vY1FiRHdY11M6lIbLchkzGL4nLGL8pnYLh6a6wc2nIWpiIEbURY2tsv7h28PN0t3NpiCxzJL8KJTUN8PFQoF9skL2LQ2SRBq0eW7NKAADDE9ijS+LZeNywbcjwBDYkOBMmYkRtpKk1lpUAEo9x77BBnUO47QIJZ/fpMtQ26BDi64Hukf72Lg6RRarrNdh9qhwAe3SdDf+aErUBrU6PHSdLAQBDu4bYuTREltuSZYxfVgJIPFsb43dIl1DI5TI7l4bIMrtOlUGrlxAb7I3YYB97F4esiIkYURs4ml+N2gYd/D3d0C2CrbEkFo1Oj31nDK2xAztyWCKJZ9epMgDAAMYvCWhXY28Y49f5MBEjagM7GysB/eKC2BpLwsk4V4V6jR7tvN3ROczX3sUhsohaq8P+MxUAgGtZkSUB7cox1CGujWP8OhsmYkRtYLcpEWtn55IQWc5YCegbGwSZjA0JJJbDeVVQa/UI8vFAfAiHdZFY6jU6HDxbCYCJmDNiIkZkY5IkmYbF8CZKImqKXzYkkHiM8dsvth0bEkg4B89WokGnR4ivErHB3vYuDlkZEzEiGztVqkJJTQM83ORIah9g7+IQWUSSJOw+bZifwGFdJCLjiIT+jF8S0C5T/LIhwRkxESOysX25hkpsz+gAKN24CSOJJbdMhbJaQ0PCNVFsSCCxSJKEfbkVAIDkWPboknhM8duB8euMmIgR2djhvCoAhkSMSDTG+E2M8OP+YSScgqp6lNY2QCGXcf8wElLGOcP8MNYhnBP/qhLZmPEm2iOKlQAST1P8shJA4slobEjoHOoLT3eOSCCxlNU2IL+yHgDQnXUIp8REjMiG9HoJR84ZKgKsyJKIMkzxy0oAiYfxSyIzNoTFBXvDz9PdzqUhW2AiRmRDZ8pVqFZr4aGQo0s4918isUiSxB5dEtrhxvhlbwKJKIMNuU6PiRiRDRlvogkRfnBX8NeNxFJUrUZJjWF+TTfOryEBGUckXMP5NSQgUyIWzfuvs2LNkMiGDuexN4HEZYzfTqE+nF9DwimvbUBeRR0A9oiRmDLyOEfX2TERI7Kho/mcn0DiaopfVgJIPMb47RDkDX/OryHBqBq0yCmtBcA6hDNjIkZkQyeKagAAXcP97FwSIssxfklkjF8S2cniWkgSEOzjgRBfpb2LQzbCRIzIRuoadKZhMZ3CuFAHiSersSLbKdTHziUhspwpfsMYvySepvsv6w/OjIkYkY2cLKmBJAGB3u4I9vGwd3GILKLXSzhZbBgW05kNCSSg7GJDRbYzK7IkIGP8siHXuTERI7IRY2tW51BfyGQyO5eGyDLnKutQp9HBXSFDhyBvexeHyGKmezArsiQgjkhwDU6ViMXFxUEmk5l9LVy48LKvqa+vx4wZMxAcHAxfX1/cfvvtKCwsbKMSkzPLZm8CCcwYv3HBPnDj1gskmKp6DYqq1QDYo0BiMvXoMn6dmtP9dZ0/fz7y8/NNX48//vhlz3/yySfx22+/4fvvv8fGjRtx7tw53HbbbW1UWnJm2RzfTQJjbwKJzHj/DfNTcsVEEo5Wp0dOCRtzXYGbvQtgbX5+foiIiGjRuZWVlVi6dClWrFiB66+/HgCwbNkydOvWDdu3b8fAgQNtWVRycicbb6KcKE4iyikxVGTjOSyGBGSsxLIhjESUV1EHjU6C0k2OqAAvexeHbMjpErGFCxfilVdeQYcOHXD33XfjySefhJvbxb/NPXv2QKPRYNSoUaZjiYmJ6NChA7Zt23bJREytVkOtVpseV1UZ9irRaDTQaDRW/G4sZ/x8e5fDWUmSHgCg0+kue40lSUJumaEiEOWv5M+jhRi/tqXTG+JXr9df8RqfbqzIRgd48ufRQoxf29JptQAM99crXeNTjcO6Ytoxfi3BGLYdjVbf9H+N9rLX+GRRNQAgpp0XdDotdDqbF88pOFL8trQMTpWI/fOf/0RycjKCgoKwdetWzJ07F/n5+XjzzTcven5BQQE8PDwQGBhodjw8PBwFBQWX/JwFCxZg3rx5zY6vXr0a3t6OMak9PT3d3kVwSsXFcgByZGQcQVpJxiXPq9UAtWrDr9fhHRuR6XSDgG2L8WsbWXkyAArk5eUhLe3MZc89dlYBQIZzxw8irfBAm5TPWTB+bSOrEgDcUFtbi7S0tMueuz3LcK+uLcpFWtrptiieU2EMW58hDzPUC9avXw+vy9TAtxQa7tUemuorxjo15wjxq1KpWnSewydizz77LF577bXLnnP06FEkJiZi9uzZpmNJSUnw8PDAI488ggULFkCptN5meHPnzjX7rKqqKsTExGDMmDHw97fv7ucajQbp6ekYPXo03N05Lt7aVpbuAcpL0aNHd6T2j73keQfPVgK7dyDcT4lbbhrThiUUG+PXtnLWZwG5JxEdHY3U1J6XPE+nl/DUzjUAJNxx43BEB3JoTEswfm1r64ki4Mh++Pj4IDV18GXP/eKTnQAqMGpgb6QmRbZNAZ0AY9h2GrR6zNmxBgAwYsQIBPld+r6asfo4cPIU+ibGITU1sa2KKDxHil/jaLkrcfhEbM6cOZg2bdplz4mPj7/o8QEDBkCr1eLUqVNISEho9nxERAQaGhpQUVFh1itWWFh42XlmSqXyoomdu7u73X/wRo5UFmcikxm6thQKxWWv77mqBgBAh2Bv/hxagfFrGwq5IX7lcvllr29x4/wEN7kMMcF+UMi5/YIlGL+2oWicZiCTya54fc+W1wMAOob582fRCoxh65NkTUMT3d3dLnt98yoM01/iQnz5c2gFR4jfln6+wydioaGhCA0NbdVr9+/fD7lcjrCwsIs+37dvX7i7u2Pt2rW4/fbbAQCZmZnIzc1FSkpKq8tMdKbc0CUdw/2XSEC5ZYb4bd/Oi0kYCadeo0NhtSER4x54JCLjPZjx6/wcPhFrqW3btmHHjh0YMWIE/Pz8sG3bNjz55JO499570a5dOwBAXl4eRo4cic8//xz9+/dHQEAAHnzwQcyePRtBQUHw9/fH448/jpSUFK6YSFflTONNNKYdb6IkHmMlgA0JJKK8ijpIEuDjoUA7b/YmkHia7sEcFu7snCYRUyqV+Oabb/Dyyy9DrVajY8eOePLJJ83mcmk0GmRmZppNoFu8eDHkcjluv/12qNVqjB07Fv/73//s8S2QE2FrFonsDBMxEtj5DQkyGXt0SSyVdRpU1hlW3GNjrvNzmkQsOTkZ27dvv+w5cXFxkCTJ7Jinpyfef/99vP/++7YsHrmYvPI6AIahXUSiYfySyJril5VYEo8xfoN8POCjdJpqOl0CF9UmsjJJkpBfaZifEMXV5khApvjlRqIkoALT/dfTziUhslxBlSERiwxg/LoCJmJEVlah0kDduHFjmL/1tk0gaiuFVYaKbLg/KwIkHmNDAuOXRFRQaVgxMYLx6xKYiBFZWUFjJTbYxwNKN4WdS0NkmfN7dNkiSyIyNiQwfklEBZWGHrEIxq9LYCJGZGUFbI0lgVXVa1Gn0QFgRYDElG+syPIeTAIyNuYyfl0DEzEiKytgaywJzNiQEOjtDk939uiSeAqrGod28R5MAjKOSGD8ugYmYkRWZpqfwJsoCYitsSSy6noNatRaAKzIkpiMQ2sZv66BiRiRlRUa59ewIksC4vwEEpmxEuvv6QZvDy79TeLhHF3XwkSMyMryq9gjRuLiil0kMg7rIpHVqrWorjf06HKeuWtgIkZkZcYeBbZmkYiMe9iwIksiakrEuAceicc4NNxX6QY/T3c7l4baAhMxIisrrjb0KIT5sSJL4mH8ksia4pd7OJJ4GL+uh4kYkRVpdXqUqzQAgBBfDzuXhshyJTUNABi/JKZSU/yyIkviYfy6HiZiRFZUpjLcROUyINCbFVkST2mtoUU2mBUBEpAxftmQQCJquv8yfl0FEzEiKzK2ZrXz9oBCLrNzaYgsZ4zhYB9WBEg8xvgNYvySgEoYvy6HiRiRFZkqsWzNIgGpGrRQNegAMIZJTCU17NElcZUyfl0OEzEiKzINK/DhTZTEY2xI8HCTw1fJPZhIPKW17NElcZVyjq7LYSJGZEXsESORlTVWYkN8PCCTcWgtiUWvl1Bey8UOSFxlpoYExq+rYCJGZEVNE8V5EyXxcKEOEllVvQZavQSAc2xITCVcrMPlMBEjsiJOFCeRcaI4icwYv36ebvBwY/WGxMPFklwP71REVlTCoYkkMA6tJZEZFzrgiAQSUYNWj8o6wz6kHJXgOpiIEVkRF+sgkbEiSyLjQh0ksvLz9yH1crdzaaitMBEjsiKueEQiY0WWRNa09Dfjl8Rj3HohyEcJOfchdRlMxIisyLjiEefYkIhKGb8ksKb4ZY8uiaeMDWEuiYkYkZVodXrUqLUAgEBv3khJPMb5CYxfElFT/HJYF4nHGL8BjF+XwkSMyEqq6rWm//t7cjNcEk+VsSLA+QkkoErGLwmM8euamIgRWYnxJurjoYCbgr9aJB5WBEhkxoYEf0/GL4mH91/XxNoikZXwJkoikySJMUxCY/ySyBi/romJGJGVGG+i/ryJkoBqG3TQ6SUAgL8Xh9aSeFiRJZGxR9c1OU0itmHDBshksot+7dq165KvGz58eLPzH3300TYsOTkLzq8hkRnj110hg5e7ws6lIbJcVZ1hni7vwSSipvhlQ5grcZqf9qBBg5Cfn2927IUXXsDatWvRr1+/y752+vTpmD9/vumxt7e3TcpIzo2tsSSy8+NXJuMeNiQe3oNJZFw10TU5TSLm4eGBiIgI02ONRoNffvkFjz/++BUrFd7e3mavJWoNVgJIZBxaSyJr0OpRp9EB4D2YxMQ6hGtymkTsQr/++itKS0tx//33X/Hcr776Cl9++SUiIiIwfvx4vPDCC5ftFVOr1VCr1abHVVVVAAzJn0ajufrCXwXj59u7HM5KkvQAAJ1O1+wal9caYsJXqeD1byXGr23p9Ib41ev1za5xWXU9AMPWC7z+rcP4tS2d1jB0S5KkZte4tKbpb7JS0fx5ahnGsO1otPqm/2u0za5xhcqwobOPu5zXv5UcKX5bWganTcSWLl2KsWPHon379pc97+6770ZsbCyioqJw8OBBPPPMM8jMzMRPP/10ydcsWLAA8+bNa3Z89erVDjOsMT093d5FcErFxXIAcmRkHEFaSYbZcxnZhucKz5xEWlq2XcrnLBi/tpGVJwOgQF5eHtLSzpg9t73I8Jy6uhxpaWl2KZ+zYPzaRlYlALihtra2WYwW1hme81JI+GvVn3YonXNhDFufIQ8zVLvXr1+PC6eClVYpAMiwf+dWFBxu69I5F0eIX5VK1aLzHD4Re/bZZ/Haa69d9pyjR48iMTHR9Pjs2bP466+/8N13313x/R9++GHT/3v27InIyEiMHDkS2dnZ6NSp00VfM3fuXMyePdv0uKqqCjExMRgzZgz8/f2v+Jm2pNFokJ6ejtGjR8Pdnd3b1raydA9QXooePbojtX+s2XNpX+8HiorQL6kHUgd2sE8BBcf4ta2c9VlA7klER0cjNbWn2XMFW04B2cfRuUMUUlOT7FNAwTF+bWvriSLgyH74+PggNXWw2XP7ciuA/TsR7OeF1NSh9imgE2AM206DVo85O9YAAEaMGIEgPy/Tc3q9hFnbDcnD+LEjEeqntEsZRedI8WscLXclDp+IzZkzB9OmTbvsOfHx8WaPly1bhuDgYNx8880Wf96AAQMAAFlZWZdMxJRKJZTK5r8k7u7udv/BGzlSWZyJTGZYaFShUDS7vjUNhvkJQb6evPZXifFrGwq5IX7lcnmz61vbYBg2E+it5LW/Soxf21C4GaosMpmsefxqDVsvBHh78NpbAWPY+iRZ09BEd3c3s+tbVa+BZAhhBPl5wZ0r114VR4jfln6+wydioaGhCA0NbfH5kiRh2bJlmDJlSqt+CPv37wcAREZGWvxacm2caEsiY/ySyLh9CImsUmWIX6WbHJ5MwlyK0+wjZrRu3Trk5OTgoYceavZcXl4eEhMTsXPnTgBAdnY2XnnlFezZswenTp3Cr7/+iilTpmDo0KFISuLQHLJM06pzDt++QdQMEzESWSU3wyWB8f7rupyuxrh06VIMGjTIbM6YkUajQWZmpmkCnYeHB9asWYO33noLtbW1iImJwe23347nn3++rYtNTsDYosUbKYmIFQESGe+/JDL26Loup0vEVqxYccnn4uLiIBkH4QKIiYnBxo0b26JY5OT0egnVasPSymyRJRFV1TfGL3t0SUBV9RyRQOJqil/WH1yN0w1NJLKHOo3ONNHW15MVARJPbWNDgq+SFQEST43asFgS45dE1BS/rD+4GiZiRFZgrMTKZIAXJ9qSgGoaY9hHyfgl8dQyfklgTQ1hTMRcDRMxIiuobVy63sfDDTKZzM6lIbKcyhjDrAiQgFQNxkSM8UviqW2MX28PNiS4GiZiRFbA1lgSXVOPGCuyJB7GL4mslvHrspiIEVmB6SbqwZsoiUej06NBa9hs1IctsiQgU48u45cEVKs2jkhg/LoaJmJEVlDLYTEkMFVjJQAAvNmYQAJijxiJjD1irouJGJEVGFuzOL6bRGRsSPBQyOHhxj8LJB5jYwJHJZCIVA2MX1fFv7hEVsAVj0hkxvj15rAYEhTn6ZLIjD26bMx1PUzEiKzAuGqiNxMxElAtW2NJYJIkcXg4Cc246icbc10PEzEiK2jqEWNrFomHvQkksnqNHnrJ8H8mYiQi44bObMx1PUzEiKygaQ8Q3kRJPJwoTiIz3n8BwNudjQkknqYeMcavq2EiRmQFrMiSyEzDutiQQAKqPW9+jVwus3NpiCzXFMO8B7saJmJEVtC0Yhdbs0g83MOGRNa0ai0rsSQmYwxzjpjrYSJGZAWcKE4iU7FHjATGYV0kMp1eQp2GW+C4KiZiRFbAHgUSWY0pfpmIkXhqOKyLBKY6b44j78Guh4kYkRVwjg2JTMV9xEhgxs1wOayLRGSMX4VcBqUbq+Wuhj9xIivgYh0kMjYkkMhq2JBAAjt/M2eZjIvNuBomYkRWUMuhXSQwxi+JTMWGMBKYigt1uDQmYkRW0NSjwBZZEo+pR5fxSwKqbeCqtSSu83vEyPUwESOyAhV7FEhgXPWTRMah4SSyplU/Gb+uiIkY0VVq0OrRoNMD4BwbEhNX/SSRNfXo8v5L4uGqn66NiRjRVaprHBYDcLI4icnYIsuKAInIuOoc778kImMdgg1hromJGNFVMm7EqJDL4K7grxSJp15j6NH1dGdFgMRjvAd7ujF+STym+OX91yWx1kh0leobb6JevImSoBjDJDJjQ4IXFzsgAbEhzLUxESO6Sk2tWfx1IjExhklk9YxfEhjvv66NP3Wiq2SsBCg5LIYEJEnSeRVZxjCJp55DE0lgao5IcGlMxIiuEofFkMg0Ogl6yfB/JmIkonptYyLGezAJiA1hro2JGNFV4rAYEplxWAzAGCYxGVedY48YiYiLdbg2Yf7q/uc//8GgQYPg7e2NwMDAi56Tm5uLcePGwdvbG2FhYXj66aeh1Wov+75lZWW455574O/vj8DAQDz44IOoqamxwXdAzorDYkhkxmExMhngwVU/SUBNix0wfkk8XKzDtQlz12poaMDEiRPx2GOPXfR5nU6HcePGoaGhAVu3bsVnn32G5cuX48UXX7zs+95zzz3IyMhAeno6fv/9d2zatAkPP/ywLb4FclLGYTEcmkgiMg2tdVdAJpPZuTREllPzHkwC46q1rk2Y3TvnzZsHAFi+fPlFn1+9ejWOHDmCNWvWIDw8HL1798Yrr7yCZ555Bi+//DI8PDyavebo0aNYtWoVdu3ahX79+gEA3n33XaSmpmLRokWIioqy2fdDzqOuwVCR5WIdJCIOiyHRcWgiiYyrJro2YRKxK9m2bRt69uyJ8PBw07GxY8fiscceQ0ZGBvr06XPR1wQGBpqSMAAYNWoU5HI5duzYgVtvvfWin6VWq6FWq02PKysrARiGOWo0Gmt9S62i0WigUqlQWloKd3d3u5bFGalV1dCrVaipqkBpaSkAoKSsFHq1CrKGWtMxah3Gr23VVldBr1ahrqbKFKsFRYZjCqWO8XuVGL+2VVVRDr1aBU2dZBarqtoa6PUSVDUVKJXV27GE4mMM206DVg+9WgUAKCstg9RQBwCoqqyEXq1CQ20V78FXyZHit7q6GoBhZeLLcZpErKCgwCwJA2B6XFBQcMnXhIWFmR1zc3NDUFDQJV8DAAsWLDD10J2vY8eOlhabBPXIW8AjFxz7BMAnHNVKAni/8et8ZwCEzLdDYYgsdAZAyAvNj/d4q61LQtQ6Xd5qfmzyRY6R+KqrqxEQEHDJ5+2aiD377LN47bXXLnvO0aNHkZiY2EYlapm5c+di9uzZpsd6vR5lZWUIDg62+xyLqqoqxMTE4MyZM/D397drWZwRr69t8fraFq+vbfH62havr+3xGtsWr69tOdL1lSQJ1dXVV5zmZNdEbM6cOZg2bdplz4mPj2/Re0VERGDnzp1mxwoLC03PXeo1RUVFZse0Wi3Kysou+RoAUCqVUCqVZscutZKjvfj7+9s9CJ0Zr69t8fraFq+vbfH62havr+3xGtsWr69tOcr1vVxPmJFdE7HQ0FCEhoZa5b1SUlLwn//8B0VFRabhhunp6fD390f37t0v+ZqKigrs2bMHffv2BQCsW7cOer0eAwYMsEq5iIiIiIiILiTMEi25ubnYv38/cnNzodPpsH//fuzfv9+059eYMWPQvXt33HfffThw4AD++usvPP/885gxY4ap92rnzp1ITExEXl4eAKBbt2644YYbMH36dOzcuRNbtmzBzJkzcdddd3HFRCIiIiIishlhFut48cUX8dlnn5keG1dBXL9+PYYPHw6FQoHff/8djz32GFJSUuDj44OpU6di/vym2ecqlQqZmZlmKxt+9dVXmDlzJkaOHAm5XI7bb78d77zzTtt9Y1amVCrx0ksvNRs6SdbB62tbvL62xetrW7y+tsXra3u8xrbF62tbIl5fmXSldRWJiIiIiIjIqoQZmkhEREREROQsmIgRERERERG1MSZiREREREREbYyJGBERERERURtjIubg3n//fcTFxcHT0xMDBgxotmn1hb7//nskJibC09MTPXv2RFpamtnzkiThxRdfRGRkJLy8vDBq1CicOHHClt+Cw7PkGn/88ccYMmQI2rVrh3bt2mHUqFHNzp82bRpkMpnZ1w033GDrb8NhWXJ9ly9f3uzaeXp6mp3DGDZnyfUdPnx4s+srk8kwbtw40zmM3yabNm3C+PHjERUVBZlMhpUrV17xNRs2bEBycjKUSiU6d+6M5cuXNzvH0vu6s7L0+v70008YPXo0QkND4e/vj5SUFPz1119m57z88svN4jcxMdGG34XjsvT6btiw4aL3h4KCArPzGL8Gll7fi91bZTIZevToYTqH8WuwYMECXHvttfDz80NYWBgmTJiAzMzMK75OxDowEzEH9u2332L27Nl46aWXsHfvXvTq1Qtjx45FUVHRRc/funUrJk+ejAcffBD79u3DhAkTMGHCBBw+fNh0zuuvv4533nkHS5YswY4dO+Dj44OxY8eivr6+rb4th2LpNd6wYQMmT56M9evXY9u2bYiJicGYMWNMe9MZ3XDDDcjPzzd9ff31123x7TgcS68vAPj7+5tdu9OnT5s9zxhuYun1/emnn8yu7eHDh6FQKDBx4kSz8xi/BrW1tejVqxfef//9Fp2fk5ODcePGYcSIEdi/fz9mzZqFhx56yCxZaM3vhLOy9Ppu2rQJo0ePRlpaGvbs2YMRI0Zg/Pjx2Ldvn9l5PXr0MIvfzZs326L4Ds/S62uUmZlpdv3CwsJMzzF+m1h6fd9++22z63rmzBkEBQU1u/8yfoGNGzdixowZ2L59O9LT06HRaDBmzBjU1tZe8jXC1oElclj9+/eXZsyYYXqs0+mkqKgoacGCBRc9/84775TGjRtndmzAgAHSI488IkmSJOn1eikiIkL673//a3q+oqJCUiqV0tdff22D78DxWXqNL6TVaiU/Pz/ps88+Mx2bOnWqdMstt1i7qEKy9PouW7ZMCggIuOT7MYbNXW38Ll68WPLz85NqampMxxi/FwdA+vnnny97zr/+9S+pR48eZscmTZokjR071vT4an9mzqol1/diunfvLs2bN8/0+KWXXpJ69eplvYI5iZZc3/Xr10sApPLy8kuew/i9uNbE788//yzJZDLp1KlTpmOM34srKiqSAEgbN2685Dmi1oHZI+agGhoasGfPHowaNcp0TC6XY9SoUdi2bdtFX7Nt2zaz8wFg7NixpvNzcnJQUFBgdk5AQAAGDBhwyfd0Zq25xhdSqVTQaDQICgoyO75hwwaEhYUhISEBjz32GEpLS61adhG09vrW1NQgNjYWMTExuOWWW5CRkWF6jjHcxBrxu3TpUtx1113w8fExO874bZ0r3YOt8TOjJnq9HtXV1c3uvydOnEBUVBTi4+Nxzz33IDc3104lFFPv3r0RGRmJ0aNHY8uWLabjjF/rWrp0KUaNGoXY2Fiz44zf5iorKwGg2e/6+UStAzMRc1AlJSXQ6XQIDw83Ox4eHt5svLZRQUHBZc83/mvJezqz1lzjCz3zzDOIiooy+8W+4YYb8Pnnn2Pt2rV47bXXsHHjRtx4443Q6XRWLb+ja831TUhIwKeffopffvkFX375JfR6PQYNGoSzZ88CYAyf72rjd+fOnTh8+DAeeughs+OM39a71D24qqoKdXV1VrnnUJNFixahpqYGd955p+nYgAEDsHz5cqxatQoffPABcnJyMGTIEFRXV9uxpGKIjIzEkiVL8OOPP+LHH39ETEwMhg8fjr179wKwzt9MMjh37hz+/PPPZvdfxm9zer0es2bNwnXXXYdrrrnmkueJWgd2s9snEwlu4cKF+Oabb7BhwwazBSXuuusu0/979uyJpKQkdOrUCRs2bMDIkSPtUVRhpKSkICUlxfR40KBB6NatGz788EO88sordiyZ81m6dCl69uyJ/v37mx1n/JIIVqxYgXnz5uGXX34xm8N04403mv6flJSEAQMGIDY2Ft999x0efPBBexRVGAkJCUhISDA9HjRoELKzs7F48WJ88cUXdiyZ8/nss88QGBiICRMmmB1n/DY3Y8YMHD582GnnyrFHzEGFhIRAoVCgsLDQ7HhhYSEiIiIu+pqIiIjLnm/815L3dGatucZGixYtwsKFC7F69WokJSVd9tz4+HiEhIQgKyvrqssskqu5vkbu7u7o06eP6doxhptczfWtra3FN99806I/7K4av61xqXuwv78/vLy8rPI7QcA333yDhx56CN99912zoUgXCgwMRNeuXRm/rdS/f3/TtWP8WockSfj0009x3333wcPD47Lnunr8zpw5E7///jvWr1+P9u3bX/ZcUevATMQclIeHB/r27Yu1a9eajun1eqxdu9asx+B8KSkpZucDQHp6uun8jh07IiIiwuycqqoq7Nix45Lv6cxac40Bw6o7r7zyClatWoV+/fpd8XPOnj2L0tJSREZGWqXcomjt9T2fTqfDoUOHTNeOMdzkaq7v999/D7VajXvvvfeKn+Oq8dsaV7oHW+N3wtV9/fXXuP/++/H111+bbbtwKTU1NcjOzmb8ttL+/ftN147xax0bN25EVlZWixrCXDV+JUnCzJkz8fPPP2PdunXo2LHjFV8jbB3YbsuE0BV98803klKplJYvXy4dOXJEevjhh6XAwECpoKBAkiRJuu+++6Rnn33WdP6WLVskNzc3adGiRdLRo0ell156SXJ3d5cOHTpkOmfhwoVSYGCg9Msvv0gHDx6UbrnlFqljx45SXV1dm39/jsDSa7xw4ULJw8ND+uGHH6T8/HzTV3V1tSRJklRdXS099dRT0rZt26ScnBxpzZo1UnJystSlSxepvr7eLt+jPVl6fefNmyf99ddfUnZ2trRnzx7prrvukjw9PaWMjAzTOYzhJpZeX6PBgwdLkyZNanac8Wuuurpa2rdvn7Rv3z4JgPTmm29K+/btk06fPi1JkiQ9++yz0n333Wc6/+TJk5K3t7f09NNPS0ePHpXef/99SaFQSKtWrTKdc6WfmSux9Pp+9dVXkpubm/T++++b3X8rKipM58yZM0fasGGDlJOTI23ZskUaNWqUFBISIhUVFbX592dvll7fxYsXSytXrpROnDghHTp0SHriiSckuVwurVmzxnQO47eJpdfX6N5775UGDBhw0fdk/Bo89thjUkBAgLRhwwaz33WVSmU6x1nqwEzEHNy7774rdejQQfLw8JD69+8vbd++3fTcsGHDpKlTp5qd/91330ldu3aVPDw8pB49ekh//PGH2fN6vV564YUXpPDwcEmpVEojR46UMjMz2+JbcViWXOPY2FgJQLOvl156SZIkSVKpVNKYMWOk0NBQyd3dXYqNjZWmT5/ukn+kjCy5vrNmzTKdGx4eLqWmpkp79+41ez/GsDlL7xHHjh2TAEirV69u9l6MX3PG5bwv/DJe06lTp0rDhg1r9prevXtLHh4eUnx8vLRs2bJm73u5n5krsfT6Dhs27LLnS5Jhu4DIyEjJw8NDio6OliZNmiRlZWW17TfmICy9vq+99prUqVMnydPTUwoKCpKGDx8urVu3rtn7Mn4NWnN/qKiokLy8vKSPPvroou/J+DW42HUFYHY/dZY6sEySJMlm3W1ERERERETUDOeIERERERERtTEmYkRERERERG2MiRgREREREVEbYyJGRERERETUxpiIERERERERtTEmYkRERERERG2MiRgREREREVEbYyJGRERERETUxpiIERGRS5s2bRomTJjQ5p+7fPlyyGQyyGQyzJo1q0WvmTZtmuk1K1eutGn5iIjIttzsXQAiIiJbkclkl33+pZdewttvvw1JktqoROb8/f2RmZkJHx+fFp3/9ttvY+HChYiMjLRxyYiIyNaYiBERkdPKz883/f/bb7/Fiy++iMzMTNMxX19f+Pr62qNoAAyJYkRERIvPDwgIQEBAgA1LREREbYVDE4mIyGlFRESYvgICAkyJj/HL19e32dDE4cOH4/HHH8esWbPQrl07hIeH4+OPP0ZtbS3uv/9++Pn5oXPnzvjzzz/NPuvw4cO48cYb4evri/DwcNx3330oKSmxuMz/+9//0KVLF3h6eiI8PBx33HHH1V4GIiJyQEzEiIiILvDZZ58hJCQEO3fuxOOPP47HHnsMEydOxKBBg7B3716MGTMG9913H1QqFQCgoqIC119/Pfr06YPdu3dj1apVKCwsxJ133mnR5+7evRv//Oc/MX/+fGRmZmLVqlUYOnSoLb5FIiKyMw5NJCIiukCvXr3w/PPPAwDmzp2LhQsXIiQkBNOnTwcAvPjii/jggw9w8OBBDBw4EO+99x769OmDV1991fQen376KWJiYnD8+HF07dq1RZ+bm5sLHx8f3HTTTfDz80NsbCz69Olj/W+QiIjsjj1iREREF0hKSjL9X6FQIDg4GD179jQdCw8PBwAUFRUBAA4cOID169eb5pz5+voiMTERAJCdnd3izx09ejRiY2MRHx+P++67D1999ZWp142IiJwLEzEiIqILuLu7mz2WyWRmx4yrMer1egBATU0Nxo8fj/3795t9nThxwqKhhX5+fti7dy++/vprREZG4sUXX0SvXr1QUVFx9d8UERE5FA5NJCIiukrJycn48ccfERcXBze3q/vT6ubmhlGjRmHUqFF46aWXEBgYiHXr1uG2226zUmmJiMgRsEeMiIjoKs2YMQNlZWWYPHkydu3ahezsbPz111+4//77odPpWvw+v//+O9555x3s378fp0+fxueffw69Xo+EhAQblp6IiOyBiRgREdFVioqKwpYtW6DT6TBmzBj07NkTs2bNQmBgIOTylv+pDQwMxE8//YTrr78e3bp1w5IlS/D111+jR48eNiw9ERHZg0ySJMnehSAiInI1y5cvx6xZs1o1/0smk+Hnn3822/+MiIjEwh4xIiIiO6msrISvry+eeeaZFp3/6KOPwtfX18alIiKitsAeMSIiIjuorq5GYWEhAMOQxJCQkCu+pqioCFVVVQCAyMhI+Pj42LSMRERkO0zEiIiIiIiI2hiHJhIREREREbUxJmJERERERERtjIkYERERERFRG2MiRkRERERE1MaYiBEREREREbUxJmJERERERERtjIkYERERERFRG2MiRkRERERE1Mb+HxsSHtENPSyGAAAAAElFTkSuQmCC',
  ];
}
