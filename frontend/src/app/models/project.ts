export type ProjectLanguage = 'Python';
export type ProjectModel = 'GPT-4o';
export type ProjectAgent = 'BasicAgent' | 'L3-Reasoning' | 'Undefined';

export interface Project {
  id: string;
  name: string;
  language: ProjectLanguage;
  model: ProjectModel;
  updatedAt: Date;
  agent: ProjectAgent;
}
