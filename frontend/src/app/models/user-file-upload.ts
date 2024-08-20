export interface UserFileUpload {
  file: File;
  progress: number;
  uploadPath: string;
  downloadURL: string | null;
  status: 'pending' | 'uploading' | 'completed' | 'error';
}
