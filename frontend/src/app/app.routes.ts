import { Routes } from '@angular/router';

import { authGuard } from './auth.guard';
//--- Components ---
import { LoginComponent } from './pages/login/login.component';
import { HomeComponent } from './pages/home/home.component';
import { SignupComponent } from './pages/signup/signup.component';

export const routes: Routes = [
  { path: '', redirectTo: '/login', pathMatch: 'full' }, // Default route
  { path: 'home', component: HomeComponent, canActivate: [authGuard] }, // route guard: must be logged in
  { path: 'login', component: LoginComponent, pathMatch: 'full' },
  { path: 'signup', component: SignupComponent },
];
