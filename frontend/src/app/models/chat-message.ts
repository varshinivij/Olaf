export interface ChatMessage {
  type: 'text' | 'code' | 'executedCode' | 'plan' | 'error' | 'image' | 'result' | "hidden";
  role: 'assistant' | 'user';
  content: string
}
