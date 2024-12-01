export interface ChatMessage {
  type: 'text' | 'code' | 'plan' | 'error' | 'image' | 'result' | "hidden";
  role: 'assistant' | 'user';
  content: string
}
