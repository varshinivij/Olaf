import { AsyncPipe, CommonModule } from '@angular/common';
import { Component, Output, EventEmitter } from '@angular/core';

import { UserService } from '../../services/user.service';

import { titleCase } from '../../utils/title-case';

export type PageName = "projects" | "data" | "settings"

interface Page {
  name: PageName,
  mat_i_outlined_icon: string
}

@Component({
  selector: 'app-sidebar',
  standalone: true,
  imports: [AsyncPipe, CommonModule],
  templateUrl: './sidebar.component.html',
  styleUrl: './sidebar.component.scss',
})
export class SidebarComponent {
  // emits a PageName whenever a new page is selected. event can be
  // listented to from parent component.
  @Output() pageSelected = new EventEmitter<PageName>();

  pages: Page[] = [
    {
      name: 'projects',
      mat_i_outlined_icon: 'inventory',
    },
    {
      name: 'data',
      mat_i_outlined_icon: 'file_copy',
    },
    {
      name: 'settings',
      mat_i_outlined_icon: 'settings',
    },
  ];

  selectedPage = this.pages[1];  // default selected page

  titleCase = titleCase;

  constructor(public userService: UserService) {}

  selectPage(page: Page) {
    if (this.selectedPage === page) {
      return;
    }

    this.selectedPage = page;
    this.pageSelected.emit(page.name);
  }
}
