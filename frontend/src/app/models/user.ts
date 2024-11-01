export interface User {
  id: string;
  email: string;
  name: string | null;
  role: string;
  status: string;
  createdAt: Date;
  updatedAt: Date;
  organization: string | null;
  profilePictureUrl: string | null;
}
