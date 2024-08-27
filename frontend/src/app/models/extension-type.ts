export const ExtensionTypeMap = {
  folder: 'folder',

  '.py': 'code',
  '.js': 'code',
  '.ts': 'code',

  '.txt': 'text',
  '.md': 'text',

  '.jpg': 'image',
  '.png': 'image',

  '.pdf': 'document',
  '.doc': 'document',
  '.docx': 'document',

  '.csv': 'dataset',
  '.json': 'dataset',

  '.h5': 'model',
} as const;

export type ExtensionType =
  (typeof ExtensionTypeMap)[keyof typeof ExtensionTypeMap] | 'unknown';

const typeImageMap: Partial<{[T in ExtensionType]: string}> = {
  code: 'assets/file-type-code.svg',
  dataset: 'assets/file-type-dataset.svg',
  folder: 'assets/file-type-folder.svg',
  model: 'assets/file-type-model.svg',
  text: 'assets/file-type-text.svg',
  unknown: 'assets/file-type-unknown.svg',
};

/**
 * Returns an image url for a given ExtensionType. If none are found
 * defaults to the image for unknown. A better way to do this later
 * may be to create a separate component that renders an img tag so we
 * can mix SVGs/mat icons/whatnot.
 *
 * @param type - The extension type.
 */
export function getImageUrlFromType(type: ExtensionType) {
  return typeImageMap[type] || typeImageMap['unknown']
}
