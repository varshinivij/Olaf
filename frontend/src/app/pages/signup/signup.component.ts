import { CommonModule } from '@angular/common';
import { Component, OnInit, OnDestroy } from '@angular/core';
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
  signupForm: FormGroup;
  errorMessage: string | null = null;
  private subscription: Subscription | undefined;

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

  ngOnInit() {
    this.subscription = this.userService.getCurrentUser().subscribe({
      next: (user) => {
        if (user) {
          console.log('Logged in: ', user);
          this.router.navigate(['/dashboard']);
        }
      },
      error: (error) => {
        this.errorMessage = UserService.convertAuthErrorToMessage(error);
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

  async loginWithGoogle() {
    try {
      await this.userService.loginWithGoogle();
    } catch (error) {
      this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error logging in with Google: ', error);
    }
  }

  async signupWithEmail() {
    try {
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
