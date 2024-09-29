import { Pipe, PipeTransform } from '@angular/core';
import { ChatMessage } from '../models/chat-message';

@Pipe({
  name: 'codemessage',
  standalone: true,
})
export class CodeMessagePipe implements PipeTransform {
  transform(messages: ChatMessage[]): ChatMessage[] {
    return messages.filter(
      (message) =>
        message.type === 'code' ||
        message.type === 'result' ||
        message.type === 'error',
    );
  }
}
