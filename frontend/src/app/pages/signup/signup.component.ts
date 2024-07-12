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
  templateUrl: './signup.component.html',
  styleUrl: './signup.component.scss',
  // placing UserService as a component-level provider creates a new instance of the
  // service rather than maintain a singleton (not what we want).
  // providers: [UserService],
})
export class SignupComponent {
  signupForm: FormGroup | any = null;

  constructor(
    private userService: UserService,
    private router: Router,
    private formBuilder: FormBuilder
  ) {
    this.signupForm = this.formBuilder.group({
      email: ['', [Validators.required, Validators.email]],
      password: ['', [Validators.required, Validators.minLength(6)]],
    });
  }

  // TODO - eventually make a new user record (done. -Cody)
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

  async signupWithEmail() {
    try {
      await this.userService.signupWithEmail(
        this.signupForm.value.email,
        this.signupForm.value.password
      );
      console.log(
        `Signed up with email: ${JSON.stringify(this.userService.currentUser())}`
      );
      this.router.navigate(['/home']);
    } catch (error) {
      console.error('Error signing up with email: ', error);
    }
  }
}
