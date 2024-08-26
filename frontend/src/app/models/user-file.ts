import { ExtensionType } from "./extension-type";

export interface UserFile {
  id: string;
  name: string;
  path: string;
  size: number;
  type: ExtensionType;
  extension: string;
  isFolder: boolean;
  uploadedOn: Date;
  storageLink: string;
}
