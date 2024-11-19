import { CommonModule } from '@angular/common';
import { Component, EventEmitter, Input, Output } from '@angular/core';

import { UserService } from '../../services/user.service';

import {
  HlmAvatarImageDirective,
  HlmAvatarComponent,
  HlmAvatarFallbackDirective,
} from '@spartan-ng/ui-avatar-helm';
import { HlmButtonDirective } from '@spartan-ng/ui-button-helm';
import { HlmIconComponent, provideIcons } from '@spartan-ng/ui-icon-helm';
import { HlmSeparatorDirective } from '@spartan-ng/ui-separator-helm';
import { HlmH4Directive } from '@spartan-ng/ui-typography-helm';

import {
  lucideClipboardCheck,
  lucideFiles,
  lucideSettings,
} from '@ng-icons/lucide';

export type PageName = 'projects' | 'data' | 'settings';

interface Page {
  name: PageName;
  lucideIcon: string;
}

@Component({
  selector: 'app-sidebar',
  standalone: true,
  imports: [
    CommonModule,

    HlmAvatarImageDirective,
    HlmAvatarComponent,
    HlmAvatarFallbackDirective,

    HlmButtonDirective,
    HlmIconComponent,
    HlmSeparatorDirective,

    HlmH4Directive,
  ],
  providers: [
    provideIcons({
      lucideClipboardCheck,
      lucideFiles,
      lucideSettings,
    }),
  ],
  templateUrl: './sidebar.component.html',
  styleUrl: './sidebar.component.scss',
})
export class SidebarComponent {
  // emits a PageName whenever a new page is selected. event can be
  // two-way binded from parent component.
  @Input() page!: PageName;
  @Output() pageChange = new EventEmitter<PageName>();

  updatePage(page: PageName): void {
    this.page = page;
    this.pageChange.emit(page);
  }

  pages: Page[] = [
    {
      name: 'projects',
      lucideIcon: 'lucideClipboardCheck',
    },
    {
      name: 'data',
      lucideIcon: 'lucideFiles',
    },
    {
      name: 'settings',
      lucideIcon: 'lucideSettings',
    },
  ];

  constructor(public userService: UserService) {}
}
