import { AfterViewInit, Component, ViewChild } from '@angular/core';
import Split from 'split.js';

import { FileStorageComponent } from '../../components/file-storage/file-storage.component';
import { ProjectsComponent } from '../../components/projects/projects.component';
import { SettingsComponent } from '../../components/settings/settings.component';

import {
  SidebarComponent,
  PageName,
} from '../../components/sidebar/sidebar.component';

import {
  HlmH2Directive,
  HlmH3Directive,
  HlmLargeDirective,
} from '@spartan-ng/ui-typography-helm';


@Component({
  selector: 'app-dashboard',
  standalone: true,
  // 2. Add it to the "imports" array
  imports: [
    FileStorageComponent,
    SidebarComponent,
    ProjectsComponent,
    HlmH2Directive,
    HlmH3Directive,
    HlmLargeDirective,
    SettingsComponent, // <--- Add this
  ],
  templateUrl: './dashboard.component.html',
  styleUrls: ['./dashboard.component.scss'],
})
export class DashboardComponent implements AfterViewInit {
  page: PageName = 'projects';

  @ViewChild(SettingsComponent) settingsComponent!: SettingsComponent;

  ngAfterViewInit() {
    Split(['#sidebar', '#main-content'], {
      sizes: [17, 83], // Initial split sizes
      minSize: [150, 700],
      gutterSize: 12,
      snapOffset: 0,
    });
  }
  openSettingsModal() {
    this.settingsComponent.openModal();
  }
}