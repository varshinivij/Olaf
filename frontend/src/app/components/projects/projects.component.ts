import { CommonModule } from '@angular/common';
import { Component } from '@angular/core';
import {
  FormBuilder,
  FormGroup,
  FormsModule,
  ReactiveFormsModule,
  Validators,
} from '@angular/forms';
import { NavigationExtras, Router } from '@angular/router';

import { ProjectService } from '../../services/project.service';
import { UserService } from '../../services/user.service';

import { Project, ProjectLanguage, ProjectModel } from '../../models/project';
import { Subscription } from 'rxjs';

@Component({
  selector: 'app-projects',
  standalone: true,
  imports: [CommonModule, FormsModule, ReactiveFormsModule],
  templateUrl: './projects.component.html',
  styleUrl: './projects.component.scss',
})
export class ProjectsComponent {
  projectSubscription?: Subscription;
  pageNumberSubscription?: Subscription;
  totalPagesSubscription?: Subscription;
  createProjectForm: FormGroup;

  projects?: Project[];
  pageNumber = 1;
  totalPages = 1;

  constructor(
    private formBuilder: FormBuilder,
    public router: Router,
    public projectService: ProjectService,
    public userService: UserService,
  ) {
    this.createProjectForm = this.formBuilder.group({
      name: ['', [Validators.required]],
      language: ['', [Validators.required]],
      model: ['', [Validators.required]],
    });
  }

  ngOnInit() {
    this.projectSubscription = this.projectService
      .getProjects()
      .subscribe((projects) => {
        this.projects = projects;
      });
    this.pageNumberSubscription = this.projectService
      .getPageNumber()
      .subscribe((number) => {
        this.pageNumber = number;
      });
    this.totalPagesSubscription = this.projectService
      .getTotalPages()
      .subscribe((total) => {
        this.totalPages = total;
      });
  }

  ngOnDestroy() {
    this.projectSubscription?.unsubscribe();
    this.pageNumberSubscription?.unsubscribe();
    this.totalPagesSubscription?.unsubscribe();
  }

  async createProject() {
    await this.projectService.createProject(
      this.createProjectForm.value.name,
      this.createProjectForm.value.language as ProjectLanguage,
      this.createProjectForm.value.model as ProjectModel,
    );
    this.createProjectForm.reset();
  }

  goToWorkspace(project: Project) {
    const navigationExtras: NavigationExtras = {
      state: { project: project },
    };
    this.router.navigate(['workspace'], navigationExtras);
  }
}
