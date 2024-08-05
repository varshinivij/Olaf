import { Component, OnDestroy, OnInit } from '@angular/core';
import { AsyncPipe, CommonModule } from '@angular/common';
import { OrderByDirection } from '@angular/fire/firestore';
import { FormsModule } from '@angular/forms';
import { Subscription } from 'rxjs';

import { posix } from 'path-browserify';

import { formatBytes } from '../../utils/format-bytes';
import { titleCase } from '../../utils/title-case';

import { FileStorageService } from '../../services/file-storage.service';
import { FileUploadService } from '../../services/file-upload.service';

import {
  UserFile,
  UserFileSortKey,
  UserFileType,
} from '../../models/user-file';
import { UserFileUpload } from '../../models/user-file-upload';

interface FilterButton {
  name: string;
  type: UserFileType;
}

type FileTypeIcons = {
  [K in UserFileType]: string;
};

@Component({
  selector: 'app-file-storage',
  standalone: true,
  imports: [AsyncPipe, CommonModule, FormsModule],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent implements OnInit, OnDestroy {
  currentPath = ['/'];
  createFolderName = '';
  searchFilter = '';
  typeFilter = null as UserFileType | null;
  sortFilter = 'name' as UserFileSortKey;
  sortDirectionFilter = 'asc' as OrderByDirection;

  // make imported function available to template
  formatBytes = formatBytes;
  titleCase = titleCase;

  private fileUploadSubscription?: Subscription;

  filterButtons: FilterButton[] = [
    { name: 'Folders', type: 'folder' },
    { name: 'Datasets', type: 'dataset' },
    { name: 'Code', type: 'code' },
    { name: 'Models', type: 'model' },
  ];

  fileTypeIcons: FileTypeIcons = {
    code: 'assets/file-type-code.svg',
    dataset: 'assets/file-type-dataset.svg',
    folder: 'assets/file-type-folder.svg',
    model: 'assets/file-type-model.svg',
    text: 'assets/file-type-text.svg',
    unknown: 'assets/file-type-unknown.svg',
  };

  constructor(
    public fileStorageService: FileStorageService,
    public fileUploadService: FileUploadService
  ) {}

  ngOnInit() {
    this.fileStorageService.setPageSize(10);
    this.fileUploadSubscription = this.fileUploadService
      .getUploadProgress()
      .subscribe((uploadProgress: UserFileUpload[]) => {
        uploadProgress.forEach((upload) => {
          if (upload.status == 'completed') {
            // probably make a toast message or something here
            console.log('Upload succeeded:', upload);
            this.fileUploadService.removeUpload(upload);
          } else if (upload.status == 'error') {
            // probably make a toast message or something here
            console.error('Upload failed:', upload);
            this.fileUploadService.removeUpload(upload);
          }
        });
      });
  }

  ngOnDestroy(): void {
    this.fileUploadSubscription?.unsubscribe();
  }

  onUploadFilesSelected(event: Event) {
    const target = event.target as HTMLInputElement;
    const files = target.files as FileList;
    this.fileUploadService.uploadFiles(
      Array.from(files),
      posix.join(...this.currentPath)
    );
  }

  toPreviousDirectory(pathIndex: number) {
    if (pathIndex < 0) {
      return;
    }

    this.currentPath.splice(pathIndex + 1);
    this.fileStorageService.setPath(posix.join(...this.currentPath));
  }

  toNewDirectory(directory: string) {
    this.currentPath.push(directory);
    this.fileStorageService.setPath(posix.join(...this.currentPath));
  }

  downloadFileOrFolder(file: UserFile) {
    // TODO: doesn't work right now. see FileStorageService
    if (file.isFolder) {
      this.fileStorageService.downloadFolder(file);
    } else {
      this.fileStorageService.downloadFile(file);
    }
  }

  deleteFileOrFolder(file: UserFile) {
    const fileFullPath = posix.join(file.path, file.name);
    this.fileStorageService.deletePath(fileFullPath);
  }

  createFolder() {
    const currentDir = posix.join(...this.currentPath);
    console.log(currentDir, this.createFolderName);
    this.fileStorageService.createFolder(this.createFolderName, currentDir);
    this.createFolderName = '';
  }

  applySearchFilter() {
    this.fileStorageService.setSearchFilter(this.searchFilter);
  }

  applyTypeFilter(type: UserFileType) {
    if (type == this.typeFilter) {
      this.typeFilter = null;
    } else {
      this.typeFilter = type;
    }

    this.fileStorageService.setTypeFilter(this.typeFilter);
  }

  applySortFilter(sort: UserFileSortKey) {
    if (sort == this.sortFilter) {
      this.sortDirectionFilter == 'asc'
        ? (this.sortDirectionFilter = 'desc')
        : (this.sortDirectionFilter = 'asc');
    } else {
      this.sortFilter = sort;
      this.sortDirectionFilter = 'asc';
    }

    this.fileStorageService.setSortField(
      this.sortFilter,
      this.sortDirectionFilter
    );
  }

  getSortDirectionMatIcon(sort: UserFileSortKey) {
    if (sort == this.sortFilter) {
      if (this.sortDirectionFilter == 'asc') {
        return 'expand_less'
      } else {
        return 'expand_more'
      }
    }
    return 'unfold_more'
  }
}
