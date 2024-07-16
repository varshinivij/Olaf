import { Component, OnInit, OnDestroy } from '@angular/core';
import { UserService } from '../services/user.service';
import { User } from '../models/user';
import { User as FirebaseUser } from '@angular/fire/auth';
import { Subscription } from 'rxjs';
import { CommonModule, DatePipe } from '@angular/common';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatProgressBarModule } from '@angular/material/progress-bar';

@Component({
  standalone: true,
  selector: 'app-profile',
  templateUrl: './profile.component.html',
  styleUrls: ['./profile.component.scss'],
  imports: [CommonModule, MatIconModule, MatButtonModule, MatProgressBarModule],
  providers: [DatePipe],
})
export class ProfileComponent implements OnInit, OnDestroy {
  user: User | null = null;
  private userSubscription: Subscription | null = null;
  selectedFile: File | null = null;
  uploadProgress: number | null = null;

  constructor(private userService: UserService) {}

  ngOnInit() {
    this.userSubscription = this.userService.user$.subscribe({
      next: async (firebaseUser) => {
        if (firebaseUser) {
          try {
            this.user = await this.userService.readUserInfo(firebaseUser);
          } catch (error) {
            console.error('Error reading user info:', error);
          }
        }
      },
      error: (err) => {
        console.error('Error fetching user data:', err);
      },
    });
  }

  onFileSelected(event: any) {
    this.selectedFile = event.target.files[0];
  }

  async onUpload() {
    if (this.selectedFile && this.user) {
      this.uploadProgress = 0;
      const downloadUrl = await this.userService.uploadProfilePicture(
        this.user.id,
        this.selectedFile,
        (progress) => (this.uploadProgress = progress)
      );
      const firebaseUser = mapToFirebaseUser(this.user);
      await this.userService.updateUserInfo(firebaseUser, {
        profilePictureUrl: downloadUrl,
      });
      this.user.profilePictureUrl = downloadUrl;
      this.selectedFile = null;
      this.uploadProgress = null;
    }
  }

  ngOnDestroy() {
    if (this.userSubscription) {
      this.userSubscription.unsubscribe();
    }
  }
}

function mapToFirebaseUser(user: User): FirebaseUser {
  return {
    uid: user.id,
    email: user.email,
    displayName: user.name,
    emailVerified: false,
    isAnonymous: false,
    metadata: {
      creationTime: '',
      lastSignInTime: '',
    },
    providerData: [],
  } as unknown as FirebaseUser;
}
