import { ChatMessage } from "./chat-message";
export interface Session {
    name?: string;
    userId: string | null;
    history: ChatMessage[];
    files: File[];
    context: string;
    id: string;
    sandboxId: string | null;
}

export const createNewSession = (): Session => ({
  name: '<untitled session>',
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
  sandboxId: null,
});
