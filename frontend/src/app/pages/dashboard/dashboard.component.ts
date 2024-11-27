import {
  AfterViewInit,
  ChangeDetectionStrategy,
  Component,
} from '@angular/core';

import Split from 'split.js';

import { FileStorageComponent } from '../../components/file-storage/file-storage.component';
import { ProjectsComponent } from '../../components/projects/projects.component';
import {
  SidebarComponent,
  PageName,
} from '../../components/sidebar/sidebar.component';

import { HlmH2Directive, HlmH3Directive, HlmLargeDirective } from '@spartan-ng/ui-typography-helm';

@Component({
  selector: 'app-dashboard',
  changeDetection: ChangeDetectionStrategy.OnPush,
  standalone: true,
  imports: [
    FileStorageComponent,
    SidebarComponent,
    ProjectsComponent,

    HlmH2Directive,
    HlmH3Directive,
    HlmLargeDirective,
  ],
  templateUrl: './dashboard.component.html',
  styleUrl: './dashboard.component.scss',
})
export class DashboardComponent implements AfterViewInit {
  page: PageName = 'projects';

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content'], {
      sizes: [17, 83], // Initial sizes of the columns in percentage
      minSize: [150, 700], // Minimum size of each column in pixels
      gutterSize: 12, // Size of the gutter (the draggable area between columns)
      snapOffset: 0,
    });
  }
}
