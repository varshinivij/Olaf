import { inject } from '@angular/core';
import { CanActivateFn, Router } from '@angular/router';
import { firstValueFrom } from 'rxjs';

import { UserService } from './services/user.service';

// using Angular's new functional route guards rather than class-based.
// does not redirect automatically if user logs out in an unauthorized page
// (i.e. if user logs out in home page, you need to redirect to login yourself)

// referenced from https://dev.to/this-is-angular/demystifying-angular-route-guards-a-beginners-guide-to-secure-navigation-597b

export const authGuard: CanActivateFn = async (route, state) => {
  const router = inject(Router);
  const userService = inject(UserService);
  const user = await firstValueFrom(userService.getCurrentUser());

  if (user) {
    return true;
  } else {
    console.log('User is not logged in, access denied');
    router.navigate(['/login']);
    return false;
  }
};
