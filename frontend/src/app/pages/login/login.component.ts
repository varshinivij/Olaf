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
  templateUrl: './login.component.html',
  styleUrl: './login.component.scss',
})
export class LoginComponent implements OnInit, OnDestroy {
  loginForm: FormGroup;
  errorMessage: string | null = null;
  private subscription: Subscription | null = null;

  constructor(
    private formBuilder: FormBuilder,
    private router: Router,
    private userService: UserService
  ) {
    this.loginForm = this.formBuilder.group({
      email: ['', [Validators.required, Validators.email]],
      password: ['', [Validators.required, Validators.minLength(6)]],
    });
  }

  ngOnInit(): void {
    if (this.userService) {
      this.subscription = this.userService.getCurrentUser().subscribe({
        next: (user) => {
          if (user) {
            console.log('Logged in: ', user);
            this.router.navigate(['/home']);
          }
        },
        error: (error) => {
          console.error('Error retrieving user data: ', error);
        },
      });
    }
  }

  ngOnDestroy() {
    this.subscription?.unsubscribe();
  }

  navigateToSignup() {
    this.router.navigate(['/signup']);
  }

  async loginWithGoogle() {
    try {
      await this.userService.loginWithGoogle();
    } catch (error) {
      // this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error logging in with Google: ', error);
    }
  }

  async loginWithEmail() {
    try {
      await this.userService.loginWithEmail(
        this.loginForm.value.email,
        this.loginForm.value.password
      );
    } catch (error) {
      // this.errorMessage = UserService.convertAuthErrorToMessage(error);
      console.error('Error logging in with email: ', error);
    }
  }
}
