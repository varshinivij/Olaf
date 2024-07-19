import { Component, OnInit, OnDestroy } from '@angular/core';
import { UserService } from '../services/user.service';
import { User } from '../models/user';
import { User as FirebaseUser } from '@angular/fire/auth';
import { Subscription } from 'rxjs';
import { CommonModule, DatePipe } from '@angular/common';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatProgressBarModule } from '@angular/material/progress-bar';
import { Timestamp } from '@angular/fire/firestore';

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
    this.userSubscription = this.userService.getCurrentUser().subscribe({
      next: (user) => {
        if (user) {
          this.user = user;
          this.user.createdAt =
            user.createdAt instanceof Timestamp
              ? user.createdAt.toDate()
              : user.createdAt;
          this.user.updatedAt =
            user.updatedAt instanceof Timestamp
              ? user.updatedAt.toDate()
              : user.updatedAt;
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
      try {
        const downloadUrl = await this.userService.uploadProfilePicture(
          this.selectedFile,
          (progress) => (this.uploadProgress = progress)
        );
        this.selectedFile = null;
        this.uploadProgress = null;

        this.user.profilePictureUrl = downloadUrl;

        await this.userService.updateUserInfo(
          {
            uid: this.user.id,
            email: this.user.email,
            displayName: this.user.name,
          } as unknown as FirebaseUser,
          {
            profilePictureUrl: downloadUrl,
          }
        );
      } catch (error) {
        console.error('Error uploading profile picture:', error);
      }
    } else {
      console.error('User ID is not defined or file is not selected.');
    }
  }

  ngOnDestroy() {
    if (this.userSubscription) {
      this.userSubscription.unsubscribe();
    }
  }
}
