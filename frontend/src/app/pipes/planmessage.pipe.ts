import { Pipe, PipeTransform } from '@angular/core';
import { ChatMessage } from '../models/chat-message';

@Pipe({
  name: 'planmessage',
  standalone: true,
})
export class PlanMessagePipe implements PipeTransform {
  transform(messages: ChatMessage[]): ChatMessage[] {
    return messages.filter((message) => message.type === 'plan');
  }
}
