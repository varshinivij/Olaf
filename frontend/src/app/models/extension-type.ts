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

const TypeLucideIconMap: Partial<{ [T in ExtensionType]: string }> = {
  code: 'lucideFileCode',
  dataset: 'lucideFileChartColumn',
  folder: 'lucideFolderOpen',
  model: 'lucideFileArchive',
  document: 'lucideFileText',
  unknown: 'lucideFileQuestion',
};

/**
 * Returns a Lucide Icon string for a given ExtensionType. If none are found
 * defaults to icon for unknown.
 *
 * @param type - The extension type.
 */
export function getLucideIconFromType(type: ExtensionType): string {
  return (TypeLucideIconMap[type] || TypeLucideIconMap['unknown'])!;
}
