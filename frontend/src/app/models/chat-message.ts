export interface ChatMessage {
    type: 'text' | 'code';
    role: 'assistant' | 'user';
    content: string;
}
