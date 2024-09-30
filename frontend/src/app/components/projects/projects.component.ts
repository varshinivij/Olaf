import { CommonModule } from '@angular/common';
import { Component } from '@angular/core';
import {
  FormBuilder,
  FormGroup,
  FormsModule,
  ReactiveFormsModule,
  Validators,
} from '@angular/forms';
import { Router } from '@angular/router';

import { UserService } from '../../services/user.service';

import { Project } from '../../models/project';

@Component({
  selector: 'app-projects',
  standalone: true,
  imports: [CommonModule, FormsModule, ReactiveFormsModule],
  templateUrl: './projects.component.html',
  styleUrl: './projects.component.scss',
})
export class ProjectsComponent {
  // Dummy data to be replaced by a ProjectService implementation.
  pageNumber = 1;
  totalPages = 1; // updated automatically
  searchFilter = ''; // ngModel variable
  createProjectForm: FormGroup;
  displayedProjects: Project[] = [
    {
      id: 'abcdefg',
      name: 'Sample Project 1',
      language: 'Python',
      model: 'GPT-4o',
      updatedAt: new Date(),
    },
  ];

  constructor(
    private formBuilder: FormBuilder,
    private router: Router,
    public userService: UserService
  ) {
    this.createProjectForm = this.formBuilder.group({
      name: ['', [Validators.required]],
      language: ['', [Validators.required]],
      model: ['', [Validators.required]],
    });
  }

  goToWorkspace() {
    this.router.navigate(['/workspace']);
  }

  createProject() {
    this.displayedProjects.push({
      id: crypto.randomUUID(),
      name: this.createProjectForm.value.name,
      language: this.createProjectForm.value.language,
      model: this.createProjectForm.value.model,
      updatedAt: new Date(),
    });

    this.createProjectForm.reset();
  }

  /*
    PAGE NAVIGATOR METHODS
  */

  previousPage() {
    if (this.pageNumber > 1) {
      this.pageNumber--;
    }
  }

  nextPage() {
    if (this.pageNumber < this.totalPages) {
      this.pageNumber++;
    }
  }

  setPage(page: number) {
    if (page >= 1 && page <= this.totalPages) {
      this.pageNumber = page;
    }
  }
}
