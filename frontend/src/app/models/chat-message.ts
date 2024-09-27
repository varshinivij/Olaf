export interface ChatMessage {
  type: 'text' | 'code' | 'plan' | 'error' | 'image' | 'result';
  role: 'assistant' | 'user';
  content: string;
  isLive?: boolean;
}
