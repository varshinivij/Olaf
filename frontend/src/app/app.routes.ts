import { Routes } from '@angular/router';

import { authGuard } from './auth.guard';
//--- Components ---
import { LoginComponent } from './pages/login/login.component';
import { HomeComponent } from './pages/home/home.component';
import { SignupComponent } from './pages/signup/signup.component';
import { ProfileComponent } from './pages/profile/profile.component';
import { FileStorageComponent } from './components/file-storage/file-storage.component';
import { TestE2bComponent } from './components/test-e2b/test-e2b.component';
import { WorkspaceComponent } from './pages/workspace/workspace.component';

export const routes: Routes = [
  { path: '', redirectTo: '/teste2b', pathMatch: 'full' }, // Default route
  { path: 'home', component: HomeComponent },
  { path: 'login', component: LoginComponent, pathMatch: 'full' },
  { path: 'signup', component: SignupComponent },
  { path: 'profile', component: ProfileComponent },
  { path: 'teststorage', component: FileStorageComponent },
  { path: 'teste2b', component: TestE2bComponent },
  { path: 'workspace', component: WorkspaceComponent},
];
