import { ComponentFixture, TestBed } from '@angular/core/testing';
import { of, throwError } from 'rxjs';
import { UserService } from '../services/user.service';
import { User } from '../models/user';
import { ProfileComponent } from './profile.component';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatProgressBarModule } from '@angular/material/progress-bar';
import { CommonModule, DatePipe } from '@angular/common';

describe('ProfileComponent', () => {
  let component: ProfileComponent;
  let fixture: ComponentFixture<ProfileComponent>;
  let userServiceSpy: jasmine.SpyObj<UserService>;
  let mockUser: User;

  beforeEach(async () => {
    mockUser = {
      id: '1',
      email: 'test@example.com',
      name: 'Test User',
      profilePictureUrl: '',
      createdAt: new Date(),
      updatedAt: new Date(),
    };

    userServiceSpy = jasmine.createSpyObj('UserService', [
      'getCurrentUser',
      'uploadProfilePicture',
      'updateUserInfo',
    ]);
    userServiceSpy.getCurrentUser.and.returnValue(of(mockUser));
    userServiceSpy.uploadProfilePicture.and.returnValue(
      Promise.resolve('https://example.com/profile.jpg')
    );

    await TestBed.configureTestingModule({
      imports: [
        CommonModule,
        MatIconModule,
        MatButtonModule,
        MatProgressBarModule,
      ],
      declarations: [ProfileComponent],
      providers: [{ provide: UserService, useValue: userServiceSpy }, DatePipe],
    }).compileComponents();

    fixture = TestBed.createComponent(ProfileComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });

  it('should fetch the current user on init', () => {
    expect(userServiceSpy.getCurrentUser).toHaveBeenCalled();
    expect(component.user).toEqual(mockUser);
  });

  it('should handle file selection', () => {
    const file = new File([''], 'test.jpg', { type: 'image/jpeg' });
    const event = { target: { files: [file] } };

    component.onFileSelected(event);

    expect(component.selectedFile).toBe(file);
  });

  it('should upload the profile picture and update user info', async () => {
    const file = new File([''], 'test.jpg', { type: 'image/jpeg' });
    component.selectedFile = file;
    component.user = mockUser;

    await component.onUpload();

    expect(userServiceSpy.uploadProfilePicture).toHaveBeenCalledWith(
      file,
      jasmine.any(Function)
    );
    expect(userServiceSpy.updateUserInfo).toHaveBeenCalledWith(
      jasmine.objectContaining({
        uid: mockUser.id,
        email: mockUser.email,
        displayName: mockUser.name,
      }),
      jasmine.objectContaining({
        profilePictureUrl: 'https://example.com/profile.jpg',
      })
    );
    expect(component.user.profilePictureUrl).toBe(
      'https://example.com/profile.jpg'
    );
  });

  it('should handle upload errors', async () => {
    const file = new File([''], 'test.jpg', { type: 'image/jpeg' });
    component.selectedFile = file;
    component.user = mockUser;
    userServiceSpy.uploadProfilePicture.and.returnValue(
      Promise.reject('Upload error')
    );

    await component.onUpload();

    expect(component.uploadProgress).toBeNull();
  });

  it('should unsubscribe on destroy', () => {
    spyOn(component.userSubscription, 'unsubscribe');

    component.ngOnDestroy();

    expect(component.userSubscription.unsubscribe).toHaveBeenCalled();
  });
});
