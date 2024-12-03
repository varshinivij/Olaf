/**
 * Returns a string with all words in title case (for Angular templates,
 * use the titlecase pipe).
 *
 * Reference: https://stackoverflow.com/a/64910248
 *
 * @param string The string to convert.
 * @returns The string in title case.
 */
export function titleCase(string: string) {
  return string.replace(/(^|\s)\S/g, function (c) {
    return c.toUpperCase();
  });
}
