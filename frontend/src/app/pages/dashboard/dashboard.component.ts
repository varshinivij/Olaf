import { AfterViewInit, Component } from '@angular/core';

import { FileStorageComponent } from '../../components/file-storage/file-storage.component';
import { SidebarComponent, PageName } from '../../components/sidebar/sidebar.component';

import Split from 'split.js';

@Component({
  selector: 'app-dashboard',
  standalone: true,
  imports: [FileStorageComponent, SidebarComponent],
  templateUrl: './dashboard.component.html',
  styleUrl: './dashboard.component.scss',
})
export class DashboardComponent implements AfterViewInit {
  selectedPage: PageName = 'data';

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content'], {
      sizes: [10, 90], // Initial sizes of the columns in percentage
      minSize: [200, 700], // Minimum size of each column in pixels
      gutterSize: 17, // Size of the gutter (the draggable area between columns)
      cursor: 'col-resize', // Cursor to show when hovering over the gutter
    });
  }

  selectPage(page: PageName) {
    this.selectedPage = page;
  }
}
