export interface UserUploadTask {
  id: string;
  file: File | null;
  name: string;
  type: 'file' | 'folder'
  progress: number;
  uploadPath: string;
  downloadURL: string | null;
  status: 'pending' | 'uploading' | 'completed' | 'error';
  onCompleted?: (uploadRef: UserUploadTask) => void;
  onError?: (uploadRef: UserUploadTask) => void;
}
