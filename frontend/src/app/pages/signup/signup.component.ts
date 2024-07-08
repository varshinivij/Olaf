import { Component } from '@angular/core';
import { UserService } from '../../services/user.service';
import { Router } from '@angular/router';
import { ReactiveFormsModule, FormBuilder, FormGroup, Validators } from '@angular/forms';


@Component({
  selector: 'app-login',
  standalone: true,
  imports: [

  ],
  templateUrl: './signup.component.html',
  styleUrl: './signup.component.scss',
  providers: [UserService]
})
export class SignupComponent {
  signupForm: FormGroup|any = null;

  constructor(
    private userService: UserService,
    private router: Router,
    private formBuilder: FormBuilder
  
  ) {
    this.signupForm = this.formBuilder.group({
      email: ['', Validators.required, Validators.email],
      password: ['', Validators.required, Validators.minLength(6)]
    });
  }

  // TODO - eventually make a new user record
  async loginWithGoogle() {
    try {
      const result = await this.userService.loginWithGoogle();
      console.log('Logged in with Google: ', result);
      this.router.navigate(['/home']);
    } catch (error) {
      console.error('Error logging in with Google: ', error);
    }
  }

  signupWithEmail() {
    this.userService.signupWithEmail(this.signupForm.value.email, this.signupForm.value.password);
  }



}
