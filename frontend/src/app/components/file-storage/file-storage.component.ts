import { Component, OnDestroy, OnInit } from '@angular/core';
import { AsyncPipe } from '@angular/common';

import { formatBytes } from '../../utils/format-bytes'
import { titleCase } from '../../utils/title-case'

import { FileStorageService } from '../../services/file-storage.service';

import { UserFile, UserFileType } from '../../models/user-file';

interface FilterButton {
  name: string,
  type: UserFileType
}

@Component({
  selector: 'app-file-storage',
  standalone: true,
  imports: [AsyncPipe],
  templateUrl: './file-storage.component.html',
  styleUrl: './file-storage.component.scss',
})
export class FileStorageComponent implements OnInit, OnDestroy {
  private currentPath = ['/'];

  // make imported function available to template
  formatBytes = formatBytes;
  titleCase = titleCase;

  filterButtons: FilterButton[] = [
    { name: 'Folders', type: 'folder' },
    { name: 'Datasets', type: 'dataset' },
    { name: 'Code', type: 'code' },
    { name: 'Models', type: 'model' },
  ];

  constructor(public fileStorageService: FileStorageService) {}

  ngOnInit() {
  }

  ngOnDestroy(): void {
  }

  activateButton() {
    // this.fileStorageService.setSearchFilter('test');
  }
}
