import { Component } from '@angular/core';
import {
    FormBuilder,
    FormGroup,
    ReactiveFormsModule,
    Validators,
} from '@angular/forms';
import { Router } from '@angular/router';
import { Subscription } from 'rxjs';

import { UserService } from '../../services/user.service';

import { HlmButtonDirective } from '@spartan-ng/ui-button-helm';
import { HlmInputDirective } from '@spartan-ng/ui-input-helm';
import { HlmH1Directive } from '@spartan-ng/ui-typography-helm';

@Component({
    selector: 'app-onboarding',
    standalone: true,
    imports: [ReactiveFormsModule, HlmButtonDirective, HlmInputDirective, HlmH1Directive],
    templateUrl: './onboarding.component.html',
    styleUrl: './onboarding.component.scss',
})
export class OnboardingComponent {
    onboardForm: FormGroup;
    errorMessage: string | null = null;
    private subscription: Subscription | null = null;

    constructor(
        private formBuilder: FormBuilder,
        private router: Router,
        private userService: UserService
    ) {
        this.onboardForm = this.formBuilder.group({
            name: ['', [Validators.required]],
        });
    }

    ngOnInit(): void {
        this.subscription = this.userService
            .getCurrentUser()
            .subscribe((user) => {
                if (user === null || user.name !== null) {
                    this.router.navigate(['/dashboard']);
                }
            });
    }

    ngOnDestroy() {
        this.subscription?.unsubscribe();
    }

    navigateToDashboard() {
        this.router.navigate(['/dashboard']);
    }

    async submitOnboardForm() {
        try {
            await this.userService.updateAccount({
                name: this.onboardForm.value.name,
            });
        } catch (error) {
            this.errorMessage =
                'Account could not be updated. Log in and try again later.';
            console.error('Error updating account: ', error);
        }
    }
}
