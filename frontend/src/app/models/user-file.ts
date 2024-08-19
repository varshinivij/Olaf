export interface UserFile {
  id: string;
  name: string;
  path: string;
  size: number;
  type: string;
  extension: string;
  isFolder: boolean;
  uploadedOn: Date;
  storageLink: string;
}
