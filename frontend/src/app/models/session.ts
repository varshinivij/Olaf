import { ChatMessage } from "./chat-message";
export interface Session {
    userId: string | null;
    history: ChatMessage[];
    files: File[];
    context: string;
    id: string;
}

export const createNewSession = (): Session => ({
  userId: '',
  history: [
    {
      type: 'text',
      role: 'assistant',
      content: 'Hello, how can I help you today?',
    },
  ],
  files: [],
  context: '',
  id: '',
});
