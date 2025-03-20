import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormBuilder, FormGroup, ReactiveFormsModule } from '@angular/forms';
import { ThemeService } from '../../services/theme.service';

@Component({
  selector: 'app-settings',
  standalone: true,
  imports: [
    CommonModule,
    ReactiveFormsModule,
  ],
  templateUrl: './settings.component.html',
  styleUrls: ['./settings.component.scss'],
})
export class SettingsComponent implements OnInit {
  @Input() showModal: boolean = false;
  @Output() closeModalEvent = new EventEmitter<void>();

  // Active tab state
  activeTab: 'general' | 'code' = 'general';

  // Collapsible sections
  pythonOpen = true;
  rOpen = false;
  shellOpen = false;
  systemOpen = false;

  // Example data
  pythonPackages: string[] = ['numpy', 'pandas', 'scipy', 'sklearn', 'matplotlib', 'seaborn', 'scanpy', 'anndata'];
  rPackages: string[] = ['dplyr', 'ggplot2', 'tidyr', 'shiny'];
  shellPackages: string[] = ['curl', 'wget', 'git'];
  systemInfo: string[] = ['OS: Ubuntu 20.04', 'Memory: 16GB', 'CPU: Intel i7'];

  // Form for general settings
  generalSettingsForm!: FormGroup;

  constructor(
    private fb: FormBuilder,
    private themeService: ThemeService
  ) {}

  ngOnInit(): void {
    // Initialize form, using localStorage theme or default to 'Light'
    this.generalSettingsForm = this.fb.group({
      projectName: [''],
      description: [''],
      agentSpecialization: ['scRNA-seq'],
      theme: [localStorage.getItem('theme-preference') || 'Light']
    });

    // Listen for theme changes in the form
    this.generalSettingsForm.get('theme')?.valueChanges.subscribe(value => {
      this.themeChanged(value);
    });
  }

  // Switch between tabs
  switchTab(tab: 'general' | 'code'): void {
    this.activeTab = tab;
  }

  // Toggle collapsible sections
  togglePython() { this.pythonOpen = !this.pythonOpen; }
  toggleR() { this.rOpen = !this.rOpen; }
  toggleShell() { this.shellOpen = !this.shellOpen; }
  toggleSystem() { this.systemOpen = !this.systemOpen; }

  // Save changes and close the modal
  onSaveChanges(): void {
    console.log('General Settings:', this.generalSettingsForm.value);
    this.closeModal();
  }

  // Dynamically apply theme
  themeChanged(theme: string): void {
    if (theme === 'Light' || theme === 'Dark') {
      this.themeService.setTheme(theme);
    }
  }

  // Open the modal
  openModal(): void {
    this.showModal = true;
  }

  // Close the modal
  closeModal(): void {
    this.closeModalEvent.emit();
  }
}

