import { Component } from '@angular/core';
import { UserService } from '../../services/user.service';
import { Router } from '@angular/router';
import {
  ReactiveFormsModule,
  FormBuilder,
  FormGroup,
  Validators,
} from '@angular/forms';

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
export class LoginComponent {
  loginForm: FormGroup | any = null;

  constructor(
    private userService: UserService,
    private router: Router,
    private formBuilder: FormBuilder
  ) {
    this.loginForm = this.formBuilder.group({
      email: ['', [Validators.required, Validators.email]],
      password: ['', [Validators.required, Validators.minLength(6)]],
    });
  }

  // eventually check if user exists - if not move to signup
  async loginWithGoogle() {
    try {
      await this.userService.loginWithGoogle();
      console.log(
        `Logged in with Google: ${JSON.stringify(this.userService.currentUser())}`
      );
      this.router.navigate(['/home']);
    } catch (error) {
      console.error('Error logging in with Google: ', error);
    }
  }

  navigateToSignup() {
    this.router.navigate(['/signup']);
  }

  async loginWithEmail() {
    try {
      await this.userService.loginWithEmail(
        this.loginForm.value.email,
        this.loginForm.value.password
      );
      console.log(
        `Logged in with email: ${JSON.stringify(this.userService.currentUser())}`
      );
      this.router.navigate(['/home']);
    } catch (error) {
      console.error('Error logging in with email: ', error);
    }
  }
}
