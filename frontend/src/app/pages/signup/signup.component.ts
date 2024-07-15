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
  templateUrl: './signup.component.html',
  styleUrl: './signup.component.scss',
  // placing UserService as a component-level provider creates a new instance of the
  // service rather than maintain a singleton (not what we want).
  // providers: [UserService],
})
export class SignupComponent implements OnInit, OnDestroy {
  signupForm: FormGroup;
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
          this.router.navigate(['/home']);
        }
      },
      error: (error) => {
        console.error('Error retrieving Firestore: ', error);
      },
    });
  }

  ngOnDestroy() {
    this.subscription?.unsubscribe();
  }

  async loginWithGoogle() {
    try {
      await this.userService.loginWithGoogle();
    } catch (error) {
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
      console.error('Error signing up with email: ', error);
    }
  }
}
