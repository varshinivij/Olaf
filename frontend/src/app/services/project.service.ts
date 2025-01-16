import { Injectable } from '@angular/core';
import {
  Firestore,
  collection,
  collectionData,
  deleteDoc,
  doc,
  getDocs,
  setDoc,
} from '@angular/fire/firestore';
import {
  BehaviorSubject,
  Observable,
  combineLatest,
  firstValueFrom,
  map,
  of,
  shareReplay,
  switchMap,
} from 'rxjs';

import { UserService } from './user.service';

import { Project, ProjectLanguage, ProjectModel, ProjectAgent } from '../models/project';
import { User } from '../models/user';

@Injectable({
  providedIn: 'root',
})
export class ProjectService {
  // behaviorsubjects can be modified to change initial settings
  private projects$: Observable<Project[]>;
  private searchFilter$ = new BehaviorSubject<string>('');
  private pageSize$ = new BehaviorSubject<number>(10);
  private pageNumber$ = new BehaviorSubject<number>(1);
  private totalPages$ = new BehaviorSubject<number>(1);

  constructor(
    private firestore: Firestore,
    private userService: UserService
  ) {
    this.projects$ = this.userService.getCurrentUser().pipe(
      switchMap((user: User | null) => {
        if (user) {
          const collectionRef = collection(
            this.firestore,
            'users',
            user.id,
            'projects'
          );
          return (collectionData(collectionRef) as Observable<Project[]>).pipe(
            map((projects: Project[]) =>
              projects.map((project: any) => {
                return {
                  ...project,
                  updatedAt: project.updatedAt.toDate(),
                } as Project;
              })
            )
          );
        } else {
          return of([]);
        }
      })
    );

    this.projects$ = combineLatest([
      this.projects$,
      this.searchFilter$,
      this.pageSize$,
      this.pageNumber$,
    ]).pipe(
      map(([projects, search, pageSize, pageNumber]) => {
        const filter = projects
          .filter((p) => p.name.toLowerCase().startsWith(search.toLowerCase()))
          // currently sorting by date
          .sort((a, b) => {
            return a.updatedAt.valueOf() - b.updatedAt.valueOf();
          });

        this.totalPages$.next(Math.max(1, Math.ceil(filter.length / pageSize)));

        const pageIndexStart = pageSize * (pageNumber - 1);
        const pageIndexLast = pageIndexStart + pageSize;
        return filter.slice(pageIndexStart, pageIndexLast);
      }),
      shareReplay(1)
    );
  }

  /**
   * Retrieves the user's projects as an Observable array. The search, path,
   * type, page, and sort filters are automatically applied.
   *
   * @returns An Observable array of the user's projects with filters applied.
   */
  getProjects(): Observable<Project[]> {
    return this.projects$;
  }

  /**
   * Retrieves the current search filter as an Observable.
   *
   * @returns An Observable of the search filter.
   */
  getSearchFilter(): Observable<string> {
    return this.searchFilter$.asObservable();
  }

  /**
   * Sets the search filter.
   *
   * @param search The search filter to apply.
   */
  setSearchFilter(search: string): void {
    this.searchFilter$.next(search);
    this.pageNumber$.next(1);
  }

  /**
   * Retrieves the current size of each page as an Observable.
   *
   * @returns An Observable of the current page size.
   */
  getPageSize(): Observable<number> {
    return this.pageSize$.asObservable();
  }

  /**
   * Sets the current size of each page.
   */
  setPageSize(size: number): void {
    this.pageSize$.next(size);
  }

  /**
   * Retrieves the current page number as an Observable.
   *
   * @returns An Observable of the current page number.
   */
  getPageNumber(): Observable<number> {
    return this.pageNumber$.asObservable();
  }

  /**
   * Sets the current page number.
   */
  setPageNumber(page: number): void {
    if (page <= 0 || page > this.totalPages$.getValue()) {
      return;
    }

    this.pageNumber$.next(page);
  }

  /**
   * Retrieves the total pages available as an Observable.
   *
   * @returns An Observable of the total pages.
   */
  getTotalPages(): Observable<number> {
    return this.totalPages$.asObservable();
  }

  /**
   * Creates and writes a project to Firestore under the path
   * users/{userId}/projects/{projectId} for the current user.
   *
   * @param name The name of the project. Cannot be blank.
   * @param language The programming language of the project.
   * @param name The LLM model of the project.
   */
  async createProject(
    name: string,
    language: ProjectLanguage,
    model: ProjectModel,
    agent: ProjectAgent,
  ): Promise<void> {
    if (!name.trim()) {
      return;
    }

    const user = await firstValueFrom(this.userService.getCurrentUser());
    const collectionRef = collection(
      this.firestore,
      'users',
      user!.id,
      'projects'
    );
    // automatically generates a new project id, no double write
    const documentRef = doc(collectionRef);

    const project: Project = {
      id: documentRef.id,
      name: name,
      language: language,
      model: model,
      agent: agent,
      updatedAt: new Date(),
    };

    await setDoc(documentRef, project);
  }

  /**
   * Updates an existing project in Firestore under the path
   * users/{userId}/projects/{projectId} for the current user.
   *
   * @param project The project to update.
   * @param data The data to update with. Cannot update id or blank names.
   */
  async updateProject(project: Project, data: Partial<Project>): Promise<void> {
    delete data.id; // do not update project id
    if (!data.name?.trim()) {
      delete data.name;
    }
    if (!Object.keys(data).length) {
      return; // return immediately if nothing to update
    }

    const user = await firstValueFrom(this.userService.getCurrentUser());
    const documentRef = doc(
      this.firestore,
      'users',
      user!.id,
      'projects',
      project.id
    );

    await setDoc(
      documentRef,
      { ...data, updatedAt: new Date() },
      { merge: true }
    );
  }

  /**
   * Deletes an existing project in Firestore under the path
   * users/{userId}/projects/{projectId} for the current user.
   *
   * @param project The project to delete.
   */
  async deleteProject(project: Project): Promise<void> {
    const user = await firstValueFrom(this.userService.getCurrentUser());
    const documentRef = doc(
      this.firestore,
      'users',
      user!.id,
      'projects',
      project.id
    );

    await deleteDoc(documentRef);
  }

  /**
   * Deletes all existing project in Firestore under the path
   * users/{userId}/projects/ for the current user.
   *
   * @param project The project to delete.
   */
  async deleteAllProjects(): Promise<void> {
    const user = await firstValueFrom(this.userService.getCurrentUser());
    const collectionRef = collection(
      this.firestore,
      'users',
      user!.id,
      'projects'
    );

    let promises: Promise<void>[] = [];
    const querySnapshot = await getDocs(collectionRef);
    querySnapshot.forEach((doc) => {
      promises.push(deleteDoc(doc.ref));
    });

    await Promise.all(promises);
  }
}
