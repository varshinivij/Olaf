export interface ChatMessage {
    type: 'text' | 'code' | 'plan' | 'error' | 'image'; 
    role: 'assistant' | 'user';
    content: string;
}
