/**
 * Special HTML function to vertically resize a textarea element
 * dynamically to match inputted text. Pass in a reference to the element.
 * Remember to set a max height attribute on the relevant textarea.
 *
 * When text is deleted, this function will set height to 0 then
 * recalculate scroll height, so set a min height attribute to prevent
 * shrinking too much when you delete text.
 *
 * @param element The reference to the textarea element.
 */
export function adjustTextareaHeight(element: HTMLTextAreaElement) {
  element.style.height = '0px';
  element.style.height = element.scrollHeight + 'px';
}
