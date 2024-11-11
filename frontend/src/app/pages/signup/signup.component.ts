import { CommonModule } from '@angular/common';
import { Component, OnDestroy, OnInit } from '@angular/core';
import {
  ReactiveFormsModule,
  FormBuilder,
  FormGroup,
  Validators,
} from '@angular/forms';
import { Router } from '@angular/router';
import { Subscription } from 'rxjs';

import { UserService } from '../../services/user.service';

@Component({
  selector: 'app-login',
  standalone: true,
  imports: [CommonModule, ReactiveFormsModule],
  templateUrl: './signup.component.html',
  styleUrl: './signup.component.scss',
})
export class SignupComponent implements OnInit, OnDestroy {
  formSubmitted = false;
  signupForm: FormGroup;
  errorMessage: string | null = null;
  private subscription: Subscription | null = null;

  constructor(
    private formBuilder: FormBuilder,
    private router: Router,
    private userService: UserService
  ) {
    this.signupForm = this.formBuilder.group({
      email: ['', [Validators.required, Validators.email]],
      password: ['', [Validators.required, Validators.minLength(6)]],
    });
  }

  ngOnInit(): void {
    this.subscription = this.userService.getCurrentUser().subscribe({
      next: (user) => {
        if (user && this.formSubmitted) {
          if (user.name !== null) {
            this.navigateToDashboard();
          } else {
            this.navigateToOnboarding();
          }
        }
      },
      error: (error) => {
        this.errorMessage = 'Failed to fetch user, try again later';
        console.error('Error retrieving user data: ', error);
      },
    });
  }

  ngOnDestroy() {
    this.subscription?.unsubscribe();
  }

  navigateToLogin() {
    this.router.navigate(['/login']);
  }

  navigateToDashboard() {
    this.router.navigate(['/dashboard'])
  }

  navigateToOnboarding() {
    this.router.navigate(['/onboarding']);
  }

  async loginWithGoogle() {
    try {
      this.formSubmitted = true;
      await this.userService.loginWithGoogle();
    } catch (error) {
      this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error logging in with Google: ', error);
    }
  }

  async signupWithEmail() {
    try {
      this.formSubmitted = true;
      await this.userService.signupWithEmail(
        this.signupForm.value.email,
        this.signupForm.value.password
      );
    } catch (error) {
      this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error signing up with email: ', error);
    }
  }
}
