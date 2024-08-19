/* TODO:
1) maybe some sort of toast message on upload success/not
2) progress indicators on folder creation/file deletion
  (tricky since folder creation is currently a cloud function.
  maybe finish loading based on when the request returns,
  but this would be different from the uploadService entirely.
  tricky to design this well.
3) a lot of logic for filtering data on the client side is here.
   perhaps moving it out to a separate utility class just like
   firestore-paginator is ideal.
*/
import { Component, OnDestroy, OnInit } from '@angular/core';
import { CommonModule } from '@angular/common';
import { OrderByDirection } from '@angular/fire/firestore';
import { FormsModule } from '@angular/forms';
import { Subscription } from 'rxjs';

import { posix } from 'path-browserify';

import { formatBytes } from '../../utils/format-bytes';
import { titleCase } from '../../utils/title-case';

import { FileStorageService } from '../../services/file-storage.service';
import { UploadService } from '../../services/upload.service';

import { ExtensionType } from '../../models/extension-type';
import { UserFileUpload } from '../../models/user-file-upload';
import { UserFile } from '../../models/user-file';

interface FilterButton {
  name: string;
  type: ExtensionType;
}

type FileTypeIcons = {
  [K in ExtensionType]: string;
};

@Component({
  selector: 'app-file-storage',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent implements OnInit, OnDestroy {
  private fileStorageSubscription?: Subscription;
  private uploadSubscription?: Subscription;

  pageSize = 10;  // setting
  pageNumber = 1;
  totalPages = 1;  // updated automatically

  currentPath = ['/'];
  createFolderName = ''; // ngModel variable
  searchFilter = ''; // ngModel variable
  typeFilter = null as ExtensionType | null;
  sortFilter = 'name' as keyof UserFile;
  sortDirectionFilter = 'asc' as OrderByDirection;

  private userFiles: UserFile[] = [];
  displayedFiles: UserFile[] = []; // will store userFiles after sort/filtering
  uploadQueue: UserFileUpload[] = [];

  // make imported util functions available to template
  // (is there really no better way to do this?)
  formatBytes = formatBytes;
  titleCase = titleCase;

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
    public uploadService: UploadService
  ) {}

  ngOnInit() {
    this.fileStorageSubscription = this.fileStorageService
      .getFiles()
      .subscribe((files: UserFile[] | null) => {
        if (files) {
          this.userFiles = files;
        } else {
          this.userFiles = [];
        }
        this.sortAndFilterFiles();
      });

    this.uploadSubscription = this.uploadService
      .getUploadProgress()
      .subscribe((uploadProgress: UserFileUpload[]) => {
        this.uploadQueue = uploadProgress;

        uploadProgress.forEach((upload) => {
          if (upload.status == 'completed') {
            // probably make a toast message or something here
            console.log('Upload succeeded:', upload);
            this.uploadService.removeUpload(upload);
          } else if (upload.status == 'error') {
            // probably make a toast message or something here
            console.error('Upload failed:', upload);
            this.uploadService.removeUpload(upload);
          }
        });
      });
  }

  ngOnDestroy(): void {
    this.fileStorageSubscription?.unsubscribe();
    this.uploadSubscription?.unsubscribe();
  }

  /*
    PAGE NAVIGATOR METHODS
  */
  previousPage() {
    if (this.pageNumber > 1) {
      this.pageNumber--;
      this.sortAndFilterFiles();
    }
  }

  nextPage() {
    if (this.pageNumber < this.totalPages) {
      this.pageNumber++;
      this.sortAndFilterFiles();
    }
  }

  setPage(page: number) {
    if (page >= 1 && page <= this.totalPages) {
      this.pageNumber = page;
      this.sortAndFilterFiles();
    }
  }

  /*
    DOWNLOAD/CREATE FOLDER UTILITY METHODS
  */
  currentPathToString() {
    return posix.join(...this.currentPath);
  }

  onUploadFilesSelected(event: Event) {
    const target = event.target as HTMLInputElement;
    const files = target.files as FileList;
    this.uploadService.uploadFiles(
      Array.from(files),
      posix.join(...this.currentPath)
    );
    // reset input file selection (if u don't, can't select same file again)
    target.value = '';
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
    // this.createFolderName already updated through ngModel
    const currentDir = posix.join(...this.currentPath);
    console.log(currentDir, this.createFolderName);
    this.fileStorageService.createFolder(this.createFolderName, currentDir);
    this.createFolderName = '';
  }

  /*
    SORT/FILTER UTILITY METHODS
  */

  /*
   * Updates displayedFiles and totalPages to match current
   * page number/sort/filter settings.
   * */
  private sortAndFilterFiles() {
    this.userFiles.sort((a, b) => {
      const valueA: any = a[this.sortFilter];
      const valueB: any = b[this.sortFilter];
      let cmp;

      if (typeof valueA === 'string') {
        cmp = valueA.localeCompare(valueB);
      } else {
        cmp = valueA - valueB;
      }

      if (this.sortDirectionFilter === 'desc') {
        cmp = -cmp;
      }

      return cmp;
    });

    const pageIndexStart = this.pageSize * (this.pageNumber - 1);
    const pageIndexLast = pageIndexStart + this.pageSize - 1;
    const currentPath = posix.join(...this.currentPath);
    let totalFilteredFiles = 0;

    this.displayedFiles = this.userFiles
      .filter((file: UserFile) => {
        if (this.searchFilter !== '') {
          if (
            // file name should start with the search query, non-case-sensitive
            !file.name.toLowerCase().startsWith(this.searchFilter.toLowerCase())
          ) {
            return false;
          }
        } else {
          // if no search filter, filter based on path
          if (file.path !== currentPath) {
            return false;
          }
        }

        // finally, filter based on type
        if (this.typeFilter !== null && file.type !== this.typeFilter) {
          return false;
        }

        totalFilteredFiles++;
        return true;
      })
      .filter((file: UserFile, index: number) => {
        // after filtering, display only the pagesize amount of results
        return index >= pageIndexStart && index <= pageIndexLast;
      });

    this.totalPages = Math.max(
      1,
      Math.ceil(totalFilteredFiles / this.pageSize)
    );
  }

  toPreviousDirectory(pathIndex: number) {
    if (pathIndex < 0) {
      return;
    }
    this.currentPath.splice(pathIndex + 1);
    this.sortAndFilterFiles();
  }

  toNextDirectory(directory: string) {
    this.currentPath.push(directory);
    this.sortAndFilterFiles();
  }

  applySearchFilter() {
    // this.searchFilter is already updated through ngModel
    // maybe restructure it to be more readable?
    this.sortAndFilterFiles();
  }

  applyTypeFilter(type: ExtensionType) {
    if (type == this.typeFilter) {
      this.typeFilter = null;
    } else {
      this.typeFilter = type;
    }
    this.sortAndFilterFiles();
  }

  applySortFilter(sort: keyof UserFile) {
    if (sort == this.sortFilter) {
      this.sortDirectionFilter == 'asc'
        ? (this.sortDirectionFilter = 'desc')
        : (this.sortDirectionFilter = 'asc');
    } else {
      this.sortFilter = sort;
      this.sortDirectionFilter = 'asc';
    }
    this.sortAndFilterFiles();
  }

  getSortDirectionMatIcon(sort: keyof UserFile) {
    if (sort == this.sortFilter) {
      if (this.sortDirectionFilter == 'asc') {
        return 'expand_less';
      } else {
        return 'expand_more';
      }
    }
    return 'unfold_more';
  }
}
