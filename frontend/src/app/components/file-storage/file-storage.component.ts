/* TODO:
1) maybe some sort of toast message on upload success/fail (spartan sonner)
*/
import { Component, OnDestroy, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { OrderByDirection } from '@angular/fire/firestore';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { Observable, Subscription } from 'rxjs';

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
import { HlmTableModule } from '@spartan-ng/ui-table-helm';
import { BrnSelectModule } from '@spartan-ng/ui-select-brain';
import { HlmSelectModule } from '@spartan-ng/ui-select-helm';
import { HlmH1Directive, HlmH2Directive, HlmH4Directive, HlmLargeDirective, HlmMutedDirective } from '@spartan-ng/ui-typography-helm';
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
  lucideCircleEllipsis,
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
  lucideUpload,
} from '@ng-icons/lucide';

interface FilterButton {
  name: string;
  type: ExtensionType;
}

@Component({
  selector: 'app-file-storage',
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

    BrnSelectModule,
    HlmSelectModule,

    HlmMutedDirective,
    HlmLargeDirective,
    HlmH2Directive,
    HlmH1Directive,
    HlmH4Directive
  ],
  providers: [
    provideIcons({
      lucideCircleEllipsis,
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
      lucideUpload,
    }),
  ],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent implements OnInit, OnDestroy {
  private uploadSubscription?: Subscription;

  filterButtons: FilterButton[] = [
    { name: 'Folder', type: 'folder' } as const,
    { name: 'Dataset', type: 'dataset' } as const,
    { name: 'Code', type: 'code' } as const,
    { name: 'Model', type: 'model' } as const,
  ];

  userUploads?: UserUploadTask[];

  // make imported util functions available to template
  formatBytes = formatBytes;
  getLucideIconFromType = getLucideIconFromType;

  constructor(
    public fileStorageService: FileStorageService,
    public uploadService: UploadService
  ) {}

  ngOnInit() {

  }

  ngOnDestroy(): void {
    this.uploadSubscription?.unsubscribe();
  }

  /*
    UPLOAD/DOWNLOAD FILE/FOLDER UTILITY METHODS
  */

  onUploadFilesSelected(event: Event) {
    // const target = event.target as HTMLInputElement;
    // const files = target.files as FileList;
    // Array.from(files).forEach((file: File) => {
    //   this.uploadService.uploadFile(
    //     file,
    //     this.currentPathString,
    //     (completedUpload) => {
    //       console.log('Upload succeeded:', completedUpload);
    //       this.uploadService.removeUpload(completedUpload);
    //     },
    //     (errorUpload) => {
    //       console.error('Upload failed:', errorUpload);
    //       this.uploadService.removeUpload(errorUpload);
    //     }
    //   );
    // }
    // );
    // reset selected files, else you can't select the same files again
    // target.value = '';
  }

  createFolder() {
    // this.createFolderName already updated through ngModel
    // this.uploadService.createNewFolder(
    //   this.createFolderName,
    //   this.currentPathString,
    //   (completedUpload) => {
    //     console.log('Upload succeeded:', completedUpload);
    //     this.uploadService.removeUpload(completedUpload);
    //   },
    //   (errorUpload) => {
    //     console.error('Upload failed:', errorUpload);
    //     this.uploadService.removeUpload(errorUpload);
    //   }
    // );
    // this.createFolderName = '';
  }

  deleteFileOrFolder(event: Event, file: UserFile) {
    /*
      Currently uses direct DOM manipulation. Another way would be
      to store the "deleted/error" state in the actual data model,
      but that doesn't make too much sense. Maybe an entirely separate
      layer that represents each table row, but this is good enough for now.

      If we did an entirely separate layer, then maybe separate all the
      filtering/pagination to a utility class. But since we're loading
      every file into memory at the moment this simple is enough.
    */
    // const buttonElement = event.target as HTMLButtonElement;
    // const fileFullPath = posix.join(file.path, file.name);
    // buttonElement.textContent = 'progress_activity';
    // buttonElement.classList.toggle('animate-spin');
    // this.fileStorageService
    //   .deletePath(fileFullPath)
    //   .then(() => {
    //     buttonElement.classList.toggle('animate-spin');
    //     buttonElement.textContent = 'delete';
    //   })
    //   .catch((error) => {
    //     buttonElement.classList.toggle('animate-spin');
    //     buttonElement.textContent = 'error';
    //   });
  }

  async downloadFileOrFolder(file: UserFile) {
    // TODO: doesn't work right now. see FileStorageService
    if (file.isFolder) {
      // this.fileStorageService.downloadFolder(file);
    } else {
      // this.fileStorageService.downloadFile(file);
    }
  }
}
