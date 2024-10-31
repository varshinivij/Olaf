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
import { HlmIconComponent } from '@spartan-ng/ui-icon-helm';

import { provideIcons } from '@ng-icons/core';
import { lucideLoaderCircle } from '@ng-icons/lucide';

@Component({
    selector: 'app-onboarding',
    standalone: true,
    imports: [
        ReactiveFormsModule,
        HlmButtonDirective,
        HlmIconComponent,
        HlmInputDirective,
        HlmH1Directive,
    ],
    providers: [provideIcons({ lucideLoaderCircle })],
    templateUrl: './onboarding.component.html',
    styleUrl: './onboarding.component.scss',
})
export class OnboardingComponent {
    loading = false;
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
        this.loading = true;
        try {
            if (this.onboardForm.value.name.trim() === '') {
                this.errorMessage = 'Please provide a valid name.';
            } else {
                await this.userService.updateAccount({
                    name: this.onboardForm.value.name,
                });
            }
        } catch (error) {
            this.errorMessage =
                'Account could not be updated. Log in and try again later.';
            console.error('Error updating account: ', error);
        }
        this.loading = false;
    }
}
