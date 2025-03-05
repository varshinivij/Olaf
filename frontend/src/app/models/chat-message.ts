export interface ChatMessage {
  type: 'text' | 'code' | 'executedCode' | 'plan' | 'error' | 'image' | 'result' | "hidden" | "terminal" | "terminalResult";
  role: 'assistant' | 'user';
  content: string
}
