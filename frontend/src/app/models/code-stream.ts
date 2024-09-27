export interface CodeStream {
    isOpen: boolean; // Indicates if a code block has started
    buffer: string;  // Temporary buffer to hold incomplete code
  }