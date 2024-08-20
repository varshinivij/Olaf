// https://stackoverflow.com/questions/196972/convert-string-to-title-case-with-javascript/64910248#64910248

export function titleCase(str: string) {
  return str.replace(/(^|\s)\S/g, function (t) {
    return t.toUpperCase();
  });
}
