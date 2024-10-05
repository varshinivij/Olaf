import { ChatMessage } from './chat-message';

export interface Session {
  id: string;
  userId: string;
  projectId: string;
  name: string;
  context: string;
  history: ChatMessage[];
  sandboxId: string | null;
}
