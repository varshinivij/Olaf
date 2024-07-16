export interface User {
  id: string;
  email: string;
  name: string;
  role: string;
  status: string;
  createdAt: Date;
  updatedAt: Date;
  organization?: string;
  profilePictureUrl?: string;
}
