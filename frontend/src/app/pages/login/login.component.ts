import { CommonModule, isPlatformServer } from '@angular/common';
import {
  Component,
  Inject,
  OnInit,
  OnDestroy,
  PLATFORM_ID,
} from '@angular/core';
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
  templateUrl: './login.component.html',
  styleUrl: './login.component.scss',
})
export class LoginComponent implements OnInit, OnDestroy {
  isServer = false;
  formSubmitted = false;
  loginForm: FormGroup;
  errorMessage: string | null = null;
  private subscription: Subscription | null = null;

  constructor(
    private formBuilder: FormBuilder,
    private router: Router,
    private userService: UserService,
    @Inject(PLATFORM_ID) platformId: Object,
  ) {
    this.isServer = isPlatformServer(platformId);
    this.loginForm = this.formBuilder.group({
      email: ['', [Validators.required, Validators.email]],
      password: ['', [Validators.required]],
    });
  }

  ngOnInit(): void {
    this.subscription = this.userService.getCurrentUser().subscribe({
      next: (user) => {
        if (user && this.formSubmitted) {
          console.log('Logged in: ', user);
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

  navigateToSignup() {
    this.router.navigate(['/signup']);
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

  async loginWithEmail() {
    try {
      this.formSubmitted = true;
      await this.userService.loginWithEmail(
        this.loginForm.value.email,
        this.loginForm.value.password,
      );
    } catch (error) {
      this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error logging in with email: ', error);
    }
  }
}
