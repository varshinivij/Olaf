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
  imports: [ReactiveFormsModule],
  templateUrl: './login.component.html',
  styleUrl: './login.component.scss',
  // placing UserService as a component-level provider creates a new instance of the
  // service rather than maintain a singleton (not what we want).
  // providers: [UserService],
})
export class LoginComponent implements OnInit, OnDestroy {
  loginForm: FormGroup;
  private subscription: Subscription | undefined;

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

  // this is my solution to redirecting reactively as user logs in. feel
  // free to modify if there's a better way. not sure about subscribing
  // directly inside event handler functions, might cause memory leaks.
  // placing it here will also automatically redirect to home if a
  // logged in user decides to visit /login through the URL.
  ngOnInit() {
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

  // prevent memory leaks
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
      console.error('Error logging in with email: ', error);
    }
  }
}
