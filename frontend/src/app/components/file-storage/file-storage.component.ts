/* TODO:
1) maybe some sort of toast message on upload success/fail (spartan sonner)
*/
import {
  ChangeDetectionStrategy,
  Component,
  OnDestroy,
  OnInit,
} from '@angular/core';
import { CommonModule } from '@angular/common';
import {
  FormBuilder,
  FormGroup,
  FormsModule,
  NgModel,
  ReactiveFormsModule,
  Validators,
} from '@angular/forms';
import { firstValueFrom, Observable, Subscription } from 'rxjs';

import { FileStorageService } from '../../services/file-storage.service';
import { UploadService } from '../../services/upload.service';

import {
  ExtensionType,
  getLucideIconFromType,
} from '../../models/extension-type';
import { UserFile } from '../../models/user-file';
import { UserUploadTask } from '../../models/user-upload-task';

import { formatBytes } from '../../utils/format-bytes';

import {
  BrnProgressComponent,
  BrnProgressIndicatorComponent,
} from '@spartan-ng/ui-progress-brain';
import {
  HlmProgressDirective,
  HlmProgressIndicatorDirective,
} from '@spartan-ng/ui-progress-helm';

import { DecimalPipe, TitleCasePipe } from '@angular/common';
import { HlmButtonModule } from '@spartan-ng/ui-button-helm';
import {
  HlmCheckboxCheckIconComponent,
  HlmCheckboxComponent,
} from '@spartan-ng/ui-checkbox-helm';
import { HlmInputDirective } from '@spartan-ng/ui-input-helm';
import { BrnMenuTriggerDirective } from '@spartan-ng/ui-menu-brain';
import { HlmIconComponent, provideIcons } from '@spartan-ng/ui-icon-helm';

import { HlmMenuModule } from '@spartan-ng/ui-menu-helm';
import {
  HlmTableModule,
  HlmTdComponent,
  HlmTrowComponent,
} from '@spartan-ng/ui-table-helm';
import { BrnSelectImports } from '@spartan-ng/ui-select-brain';
import { HlmSelectImports } from '@spartan-ng/ui-select-helm';
import {
  HlmH1Directive,
  HlmH2Directive,
  HlmH4Directive,
  HlmLargeDirective,
  HlmMutedDirective,
} from '@spartan-ng/ui-typography-helm';
import {
  HlmPaginationContentDirective,
  HlmPaginationDirective,
  HlmPaginationEllipsisComponent,
  HlmPaginationItemDirective,
  HlmPaginationLinkDirective,
  HlmPaginationNextComponent,
  HlmPaginationPreviousComponent,
} from '@spartan-ng/ui-pagination-helm';
import { HlmNumberedPaginationComponent } from '@spartan-ng/ui-pagination-helm';
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
import {
  HlmMenuComponent,
  HlmMenuItemDirective,
  HlmMenuItemIconDirective,
} from '@spartan-ng/ui-menu-helm';

import {
  lucideArrowDownUp,
  lucideArrowDownWideNarrow,
  lucideArrowUpNarrowWide,
  lucideChevronRight,
  lucideCircleEllipsis,
  lucideEllipsis,
  lucideFile,
  lucideFileArchive,
  lucideFileChartColumn,
  lucideFileCode,
  lucideFileQuestion,
  lucideFileText,
  lucideFileUp,
  lucideFolderOpen,
  lucideFolderPlus,
  lucideFolderUp,
  lucideSearch,
  lucideTrash2,
  lucideUpload,
} from '@ng-icons/lucide';
import { HlmSpinnerComponent } from '@spartan-ng/ui-spinner-helm';

import { toast } from 'ngx-sonner';
import { HlmToasterComponent } from '@spartan-ng/ui-sonner-helm';

interface FilterOption {
  name: string;
  type: ExtensionType;
}

@Component({
  selector: 'app-file-storage',
  changeDetection: ChangeDetectionStrategy.OnPush,
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    ReactiveFormsModule,

    BrnProgressComponent,
    BrnProgressIndicatorComponent,
    HlmProgressDirective,
    HlmProgressIndicatorDirective,

    BrnMenuTriggerDirective,
    HlmMenuModule,

    HlmPaginationContentDirective,
    HlmPaginationDirective,
    HlmPaginationEllipsisComponent,
    HlmPaginationItemDirective,
    HlmPaginationLinkDirective,
    HlmPaginationNextComponent,
    HlmPaginationPreviousComponent,

    HlmNumberedPaginationComponent,

    HlmTableModule,

    HlmButtonModule,

    DecimalPipe,
    TitleCasePipe,
    HlmIconComponent,
    HlmInputDirective,

    HlmCheckboxCheckIconComponent,
    HlmCheckboxComponent,

    BrnSelectImports,
    HlmSelectImports,

    HlmMutedDirective,
    HlmLargeDirective,
    HlmH2Directive,
    HlmH1Directive,
    HlmH4Directive,

    BrnDialogContentDirective,
    BrnDialogCloseDirective,
    BrnDialogTriggerDirective,
    HlmDialogComponent,
    HlmDialogContentComponent,
    HlmDialogFooterComponent,
    HlmDialogHeaderComponent,
    HlmDialogTitleDirective,

    BrnMenuTriggerDirective,
    HlmMenuComponent,
    HlmMenuItemDirective,
    HlmMenuItemIconDirective,

    HlmSpinnerComponent,
    HlmToasterComponent,
  ],
  providers: [
    provideIcons({
      lucideArrowDownUp,
      lucideArrowDownWideNarrow,
      lucideArrowUpNarrowWide,
      lucideChevronRight,
      lucideCircleEllipsis,
      lucideEllipsis,
      lucideFile,
      lucideFileArchive,
      lucideFileChartColumn,
      lucideFileCode,
      lucideFileQuestion,
      lucideFileText,
      lucideFileUp,
      lucideFolderOpen,
      lucideFolderPlus,
      lucideFolderUp,
      lucideSearch,
      lucideTrash2,
      lucideUpload,
    }),
  ],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent {
  createFolderForm: FormGroup;

  filterOptions: FilterOption[] = [
    { name: 'Code', type: 'code' } as const,
    { name: 'Dataset', type: 'dataset' } as const,
    { name: 'Document', type: 'document' } as const,
    { name: 'Folder', type: 'folder' } as const,
    { name: 'Model', type: 'model' } as const,
    { name: 'Unknown', type: 'unknown' } as const,
  ];

  // make imported util functions available to template
  formatBytes = formatBytes;
  getLucideIconFromType = getLucideIconFromType;

  constructor(
    private formBuilder: FormBuilder,
    public fileStorageService: FileStorageService,
    public uploadService: UploadService
  ) {
    this.createFolderForm = this.formBuilder.group({
      name: ['', [Validators.required]],
    });
  }

  async onFileUploadSelected(event: Event) {
    const target = event.target as HTMLInputElement;
    const files = target.files as FileList;
    const path = await firstValueFrom(this.fileStorageService.getPathFilter());

    Array.from(files).forEach((file: File) => {
      this.uploadService.uploadFile(
        file,
        path,
        (completedUpload) => {
          this.uploadService.removeUpload(completedUpload);
        },
        (errorUpload) => {
          this.uploadService.removeUpload(errorUpload);
          toast.error(`File "${file.name}" upload failed, try again later.`);
        }
      );
    });

    // reset selected files, else files remain selected
    target.value = '';
  }

  async createFolder() {
    const path = await firstValueFrom(this.fileStorageService.getPathFilter());

    this.uploadService.createNewFolder(
      this.createFolderForm.value.name,
      path,
      (completedUpload) => {
        this.uploadService.removeUpload(completedUpload);
      },
      (errorUpload) => {
        this.uploadService.removeUpload(errorUpload);
        toast.error(`Folder creation failed, try again later.`);
      }
    );

    this.createFolderForm.reset();
  }

  async deleteItem(
    file: UserFile,
    fileRow: HlmTrowComponent,
    fileData: HlmTdComponent
  ) {
    try {
      fileRow.muted.set(true);
      fileData.fileStorageComponentDeleting.set(true);
      await this.fileStorageService.deletePath([file.path, file.name]);
    } catch (error) {
      fileRow.muted.set(false);
      fileData.fileStorageComponentDeleting.set(false);
      toast.error(`File "${file.name}" failed to delete, try again later.`);
    }
  }

  async downloadItem(file: UserFile) {
    // TODO: doesn't work right now. see FileStorageService
    if (file.isFolder) {
      // this.fileStorageService.downloadFolder(file);
    } else {
      // this.fileStorageService.downloadFile(file);
    }
  }
}
