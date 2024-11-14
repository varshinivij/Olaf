import { AuthGuard, redirectUnauthorizedTo } from '@angular/fire/auth-guard';
import { Routes } from '@angular/router';

//--- Components ---
import { LoginComponent } from './pages/login/login.component';
import { SignupComponent } from './pages/signup/signup.component';
import { OnboardingComponent } from './pages/onboarding/onboarding.component';
import { DashboardComponent } from './pages/dashboard/dashboard.component';
import { WorkspaceComponent } from './pages/workspace/workspace.component';




const redirectUnauthorizedToLogin = () => redirectUnauthorizedTo(['/login']);

export const routes: Routes = [
  {
    path: '', // Default route
    component: LoginComponent,
    pathMatch: 'full',
  },
  { path: 'login', component: LoginComponent, pathMatch: 'full' },
  { path: 'signup', component: SignupComponent },
  {
    path: 'onboarding',
    component: OnboardingComponent,
    canActivate: [AuthGuard],
    data: { authGuardPipe: redirectUnauthorizedToLogin },
  },
  {
    path: 'dashboard',
    component: DashboardComponent,
    canActivate: [AuthGuard],
    data: { authGuardPipe: redirectUnauthorizedToLogin },
  },
  {
    path: 'workspace',
    component: WorkspaceComponent,
    canActivate: [AuthGuard],
    data: { authGuardPipe: redirectUnauthorizedToLogin },
  },
  { path: '**', redirectTo: 'dashboard' }, // wildcard route (page not found)
];
