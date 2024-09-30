/**
 * Special HTML event handler to vertically resize a textarea element
 * dynamically to match inputted text. Attach to the (input) handler.
 * Remember to set a max height attribute on the relevant textarea.
 *
 * When text is deleted, this function will set height to 0 then
 * recalculate scroll height, so set a min height attribute to prevent
 * shrinking too much when you delete text.
 *
 * @param $event The `$event` object from the textarea's (input) handler.
 */
export function adjustTextareaHeight($event: Event) {
  const textArea = $event.target as HTMLTextAreaElement;
  textArea.style.height = '0px';
  textArea.style.height = textArea.scrollHeight + 'px';
}
