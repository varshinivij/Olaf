/**
 * Returns whether the contents of two arrays are equal (shallow equality).
 *
 * Reference: https://stackoverflow.com/a/16436975
 *
 * @param arr1 The first array.
 * @param arr2 The second array.
 * @returns Whether the contents of the arrays are equal.
 */
export function arraysEqual(arr1: any[] | null, arr2: any[] | null) {
  if (arr1 === arr2) {
    return true;
  }
  if (arr1 == null || arr2 == null) {
    return false;
  }
  if (arr1.length !== arr2.length) {
    return false;
  }

  for (let i = 0; i < arr1.length; ++i) {
    if (arr1[i] !== arr2[i]) {
      return false;
    }
  }
  return true;
}
