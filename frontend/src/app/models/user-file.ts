export interface UserFile {
  id: string;
  name: string;
  path: string;
  size: number;
  type: UserFileType;
  isFolder: boolean;
  uploadedOn: Date;
  storageLink: string;
}

// these types are defined within the Python Cloud Functions. Finding a better way
// to consolidate it would be nice.
export type UserFileType = 'code' | 'dataset' | 'text' | 'model' | 'folder' | 'unknown';

// allowing more sort keys would require creating more composite indexes in Firestore
export type UserFileSortKey = 'name' | 'type' | 'uploadedOn'
