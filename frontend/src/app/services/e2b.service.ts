import { Injectable } from '@angular/core';
import { CodeInterpreter } from '@e2b/code-interpreter';

@Injectable({
  providedIn: 'root'
})

/**
 * This logic will be transfored to a google cloud function to protect API keys
 * Instead a client will be given a sandbox id
 * They will then reach out to our api provide their auth and sandbox id
 * If sandbox id is valid we will make the connection
 * Otherwise, spin up a new sandbox and give them the updated id
 */
export class E2bService {
  private sandbox: any;

  async createSandbox() {
    this.sandbox = await CodeInterpreter.create();
  }

  async execCode(code: string): Promise<string> {
    if (!this.sandbox) {
      throw new Error('Sandbox not created. Call createSandbox() first.');
    }
    const execution = await this.sandbox.notebook.execCell(code);
    return execution.text;
  }

  async closeSandbox() {
    if (this.sandbox) {
      await this.sandbox.close();
      this.sandbox = null;
    }
  }
}