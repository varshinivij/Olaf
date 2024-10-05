export const ExtensionTypeMap = {
  folder: 'folder',

  '.py': 'code',
  '.js': 'code',
  '.ts': 'code',
  '.jsx': 'code',
  '.tsx': 'code',
  '.html': 'code',
  '.css': 'code',
  '.scss': 'code',
  '.java': 'code',
  '.cpp': 'code',
  '.c': 'code',
  '.php': 'code',
  '.rb': 'code',
  '.swift': 'code',
  '.go': 'code',
  '.rust': 'code',
  '.lua': 'code',
  '.pl': 'code',
  '.sh': 'code',
  '.bat': 'code',
  '.ps1': 'code',
  '.yaml': 'code',

  '.jpg': 'dataset',
  '.jpeg': 'dataset',
  '.png': 'dataset',
  '.gif': 'dataset',
  '.bmp': 'dataset',
  '.svg': 'dataset',
  '.tiff': 'dataset',
  '.raw': 'dataset',
  '.ico': 'dataset',

  '.csv': 'dataset',
  '.tsv': 'dataset',
  '.json': 'dataset',
  '.xml': 'dataset',

  '.h5': 'model',
  '.pth': 'model',
  '.onnx': 'model',
  '.pb': 'model',
  '.tflite': 'model',
  '.mlmodel': 'model',
  '.caffemodel': 'model',
  '.pt': 'model',
  '.ckpt': 'model',
  '.hdf5': 'model',
  '.joblib': 'model',
  '.sav': 'model',
  '.mdl': 'model',
  '.mod': 'model',
  '.rds': 'model',
  '.rdata': 'model',
  '.sas7bdat': 'model',
  '.feather': 'model',
  '.arff': 'model',
  '.mat': 'model',
  '.npy': 'model',
  '.npz': 'model',
  '.pkl': 'model',
  '.pkl.gz': 'model',
  '.pkl.bz2': 'model',
  '.pkl.xz': 'model',
  '.pkl.zst': 'model',
  '.pkl.sz': 'model',
  '.pkl.lz4': 'model',
  '.pkl.lzma': 'model',
  '.pkl.brotli': 'model',
  '.pkl.zlib': 'model',
  '.pkl.snappy': 'model',
  '.pkl.lzo': 'model',
  '.pkl.z': 'model',
  '.pkl.zpaq': 'model',

  '.pdf': 'document',
  '.doc': 'document',
  '.docx': 'document',
  '.xls': 'document',
  '.xlsx': 'document',
  '.ppt': 'document',
  '.pptx': 'document',
  '.odt': 'document',
  '.ods': 'document',
  '.odp': 'document',
  '.txt': 'document',
  '.md': 'document',
} as const;

export type ExtensionType =
  | (typeof ExtensionTypeMap)[keyof typeof ExtensionTypeMap]
  | 'unknown';

const typeImageMap: Partial<{ [T in ExtensionType]: string }> = {
  code: 'assets/file-type-code.svg',
  dataset: 'assets/file-type-dataset.svg',
  folder: 'assets/file-type-folder.svg',
  model: 'assets/file-type-model.svg',
  document: 'assets/file-type-text.svg',
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
  return typeImageMap[type] || typeImageMap['unknown'];
}

const TypeLucideIconMap: Partial<{ [T in ExtensionType]: string }> = {
  code: 'lucideFileCode',
  dataset: 'lucideFileChartColumn',
  folder: 'lucideFolder',
  model: 'lucideFileArchive',
  document: 'lucideFileText',
  unknown: 'lucideFileQuestion',
};

/**
 * Returns a Lucide Icon string for a given ExtensionType. If none are found
 * defaults to icon for unknown. This will eventually replace the raw svg
 * method when the file storage system also adopts it.
 *
 * @param type - The extension type.
 */
export function getLucideIconFromType(type: ExtensionType): string {
  return (TypeLucideIconMap[type] || TypeLucideIconMap['unknown'])!;
}
