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

import {
  ExtensionType,
  getImageUrlFromType,
} from '../../models/extension-type';
import { UserFile } from '../../models/user-file';

interface FilterButton {
  name: string;
  type: ExtensionType;
}

@Component({
  selector: 'app-file-storage',
  standalone: true,
  imports: [CommonModule, FormsModule],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent implements OnInit, OnDestroy {
  private fileStorageSubscription?: Subscription;

  // Settings/default values, can be modified
  pageSize = 10; // setting
  typeFilter = null as ExtensionType | null;
  sortFilter = 'name' as keyof UserFile;
  sortDirectionFilter = 'asc' as OrderByDirection;
  filterButtons: FilterButton[] = [
    { name: 'Folders', type: 'folder' } as const,
    { name: 'Datasets', type: 'dataset' } as const,
    { name: 'Code', type: 'code' } as const,
    { name: 'Models', type: 'model' } as const,
  ];

  // Do not modify
  pageNumber = 1;
  totalPages = 1; // updated automatically
  currentPathArray = ['/'];
  currentPathString = '/';
  createFolderName = ''; // ngModel variable
  searchFilter = ''; // ngModel variable
  private userFiles: UserFile[] = [];
  displayedFiles: UserFile[] = []; // will store userFiles after sort/filtering

  // make imported util functions available to template
  formatBytes = formatBytes;
  getImageUrlFromType = getImageUrlFromType;
  titleCase = titleCase;

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
  }

  ngOnDestroy(): void {
    this.fileStorageSubscription?.unsubscribe();
  }

  /*
    UPLOAD/DOWNLOAD FILE/FOLDER UTILITY METHODS
  */

  onUploadFilesSelected(event: Event) {
    const target = event.target as HTMLInputElement;
    const files = target.files as FileList;

    Array.from(files).forEach((file: File) => {
      this.uploadService.uploadFile(
        file,
        this.currentPathString,
        (completedUpload) => {
          console.log('Upload succeeded:', completedUpload);
          this.uploadService.removeUpload(completedUpload);
        },
        (errorUpload) => {
          console.error('Upload failed:', errorUpload);
          this.uploadService.removeUpload(errorUpload);
        }
      );
    });

    // reset selected files, else you can't select the same files again
    target.value = '';
  }

  createFolder() {
    // this.createFolderName already updated through ngModel
    this.uploadService.createNewFolder(
      this.createFolderName,
      this.currentPathString,
      (completedUpload) => {
        console.log('Upload succeeded:', completedUpload);
        this.uploadService.removeUpload(completedUpload);
      },
      (errorUpload) => {
        console.error('Upload failed:', errorUpload);
        this.uploadService.removeUpload(errorUpload);
      }
    );
    this.createFolderName = '';
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
    const buttonElement = event.target as HTMLButtonElement;
    const fileFullPath = posix.join(file.path, file.name);

    buttonElement.textContent = 'progress_activity';
    buttonElement.classList.toggle('animate-spin');

    this.fileStorageService
      .deletePath(fileFullPath)
      .then(() => {
        buttonElement.classList.toggle('animate-spin');
        buttonElement.textContent = 'delete';
      })
      .catch((error) => {
        buttonElement.classList.toggle('animate-spin');
        buttonElement.textContent = 'error';
      });
  }

  async downloadFileOrFolder(file: UserFile) {
    // TODO: doesn't work right now. see FileStorageService
    if (file.isFolder) {
      this.fileStorageService.downloadFolder(file);
    } else {
      this.fileStorageService.downloadFile(file);
    }
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

      const cmp =
        typeof valueA === 'string'
          ? valueA.localeCompare(valueB)
          : valueA - valueB;

      if (this.sortDirectionFilter === 'desc') {
        return -cmp;
      }
      return cmp;
    });

    const pageIndexStart = this.pageSize * (this.pageNumber - 1);
    const pageIndexLast = pageIndexStart + this.pageSize - 1;
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
          if (file.path !== this.currentPathString) {
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
    this.currentPathArray.splice(pathIndex + 1);
    this.currentPathString = posix.join(...this.currentPathArray);
    this.pageNumber = 1;
    this.sortAndFilterFiles();
  }

  toNextDirectory(directory: string) {
    this.currentPathArray.push(directory);
    this.currentPathString = posix.join(...this.currentPathArray);
    this.pageNumber = 1;
    this.sortAndFilterFiles();
  }

  applySearchFilter() {
    // this.searchFilter is already updated through ngModel
    // maybe restructure to be more readable? unfamiliar w/ angular design patterns.
    this.pageNumber = 1;
    this.toPreviousDirectory(0);
    this.sortAndFilterFiles();
  }

  applyTypeFilter(type: ExtensionType) {
    if (type == this.typeFilter) {
      this.typeFilter = null;
    } else {
      this.typeFilter = type;
    }
    this.pageNumber = 1;
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
    this.pageNumber = 1;
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
