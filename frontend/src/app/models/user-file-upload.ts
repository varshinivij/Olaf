export interface UserFileUpload {
  file: File;
  progress: number;
  downloadURL: string | null;
  status: 'pending' | 'uploading' | 'completed' | 'error';
}
