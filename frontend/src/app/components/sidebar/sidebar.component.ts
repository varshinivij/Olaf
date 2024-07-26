import { AsyncPipe, CommonModule } from '@angular/common';
import { Component } from '@angular/core';

import { UserService } from '../../services/user.service';

interface Page {
  name: string,
  icon: string
}

@Component({
  selector: 'app-sidebar',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './sidebar.component.html',
  styleUrl: './sidebar.component.scss',
})
export class SidebarComponent {
  pages: Page[] = [
    {
      name: 'Projects',
      icon: 'inventory',
    },
    {
      name: 'Data',
      icon: 'file_copy',
    },
    {
      name: 'Settings',
      icon: 'settings',
    },
  ];

  selectedPageIndex = 0;

  constructor(public userService: UserService) {}
}
