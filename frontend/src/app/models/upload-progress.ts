export interface UploadProgress {
  file: File;
  progress: number;
  downloadURL: string | null;
  status: 'pending' | 'uploading' | 'completed' | 'error';
}
