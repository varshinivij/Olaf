export const ExtensionTypeMap: Record<string, string> = {
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
};

export type ExtensionType = typeof ExtensionTypeMap[keyof typeof ExtensionTypeMap]
