import { Injectable } from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class ThemeService {
  private themeKey = 'theme-preference';

  constructor() {
    this.loadTheme(); // Load saved theme on app startup
  }

  // Apply and persist the selected theme
  setTheme(theme: 'Light' | 'Dark'): void {
    if (theme === 'Dark') {
      document.body.classList.add('dark-mode');
    } else {
      document.body.classList.remove('dark-mode');
    }

    // Save to localStorage so it persists
    localStorage.setItem(this.themeKey, theme);
  }

  // Load and apply the saved theme from localStorage on startup
  loadTheme(): void {
    const savedTheme = (localStorage.getItem(this.themeKey) as 'Light' | 'Dark') || 'Light';
    this.setTheme(savedTheme);
  }
}