export type ProjectLanguage = 'Python'
export type ProjectModel = 'GPT-4o'

export interface Project {
  id: string;
  name: string;
  language: ProjectLanguage;
  model: ProjectModel;
  updatedAt: Date;
}
