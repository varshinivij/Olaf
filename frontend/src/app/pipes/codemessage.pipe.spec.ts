import { CodeMessagePipe } from './codemessage.pipe';
import { ChatMessage } from '../models/chat-message';

describe('CodemessagePipe', () => {
  it('create an instance', () => {
    const pipe = new CodeMessagePipe();
    expect(pipe).toBeTruthy();
  });
});
