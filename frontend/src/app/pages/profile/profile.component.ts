import { Component, OnInit, OnDestroy } from '@angular/core';
import { UserService } from '../../services/user.service';
import { User } from '../../models/user';
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
    this.userSubscription = this.userService.getCurrentUser().subscribe({
      next: (user) => {
        if (user) {
          this.user = user;
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
