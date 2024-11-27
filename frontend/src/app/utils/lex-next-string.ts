/**
 * [CURRENTLY UNUSED] Returns the lexicographically next string from a given string.
 * Used to Firebase query whether a string starts with a certain substring.
 * Reference: https://stackoverflow.com/q/46573804
 *
 * @param string The string to reference.
 * @returns A string of equal length lexicographically adjacent to the given string.
 */
export function lexNextString(string: string): string {
  const length = string.length;
  const stringFront = string.slice(0, length - 1);
  return stringFront + String.fromCharCode(string.charCodeAt(length - 1) + 1);
}
