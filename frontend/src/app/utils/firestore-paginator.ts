// Adapted from https://github.com/ibnYusrat/angularfire-paginator/
// This was originally a utility class that could be used directly in
// templates, but I am modifying it to assist with pagination while
// outsides services handle complex queries.

// Next and previous pagination doesn't seem to work.

// This code comes with ABSOLUTELY NO WARRANTY.

import {
  Firestore,
  QueryDocumentSnapshot,
  QueryConstraint,
  collection,
  collectionData,
  endBefore,
  limit,
  limitToLast,
  query,
  startAt,
} from '@angular/fire/firestore';
import { BehaviorSubject, Observable, filter, switchMap, tap } from 'rxjs';

export type PaginatorActions =
  | 'current'
  | 'first'
  | 'prev'
  | 'next'
  | 'last'
  | 'reset';

export class FirestorePaginator<T> {
  private items$: Observable<T[]>;
  private paging$: BehaviorSubject<PaginatorActions | null>;
  private firstItemId?: string; // the id of the first record of this query on the first page, used as marker for enabling previous
  private prevAnchor?: QueryDocumentSnapshot<any>; // anchor for querying previous page
  private nextAnchor?: QueryDocumentSnapshot<any>; // anchor for querying next page

  private firstEnabled = new BehaviorSubject(false);
  private lastEnabled = new BehaviorSubject(false);
  private nextEnabled = new BehaviorSubject(false);
  private previousEnabled = new BehaviorSubject(false);

  /**
   * AngularFire-based Firestore Paginator. Utility class to assist in paginating
   * Firestore results from a given collection. Returns items$ as an Observable
   * through this.getItems() and automatically updates through its public methods.
   *
   * Expects the documents in the collection to have an "id" field representing a unique
   * identifier in the collection.
   *
   * @param fs AngularFire Firestore service instance
   * @param path Firebase collection data path
   * @param pageSize Elements per page
   * @param queries Array of custom QueryConstraints (to allow complex queries)
   * @param stalled Sets stalled flag. Useful when you have to wait for other data to load before displaying paginator. The paginator will query after resume() is called.
   * @param debug Turns console logs on or off.
   */

  constructor(
    fs: Firestore,
    private readonly path: string,
    private pageSize: number,
    private queries: QueryConstraint[] | null = null,
    private stalled: boolean = false,
    private debug: boolean = false
  ) {
    this.paging$ = new BehaviorSubject('first' as PaginatorActions | null);

    this.items$ = this.paging$.pipe(
      tap((action) => {
        this.trace('start ------------------------ ');
        // disable all navigation buttons during query
        this.firstEnabled.next(false);
        this.lastEnabled.next(false);
        this.nextEnabled.next(false);
        this.previousEnabled.next(false);

        this.trace('page size ', this.pageSize, this.path, action);
      }),
      filter(() => {
        return this.stalled === false;
      }),
      switchMap(
        (pagingAction) => {
          this.trace('switchMap ------------------------ ');

          let queryConstraints = this.applyQueries([]);

          switch (pagingAction) {
            case 'current':
              queryConstraints = this.queryCurrent(queryConstraints);
              break;
            case 'next':
              queryConstraints = this.queryNext(queryConstraints);
              break;
            case 'prev':
              queryConstraints = this.queryPrev(queryConstraints);
              break;
            case 'last':
              queryConstraints = this.queryLast(queryConstraints);
              break;
            case 'first':
            case 'reset':
            default:
              queryConstraints = this.queryFirst(queryConstraints);
              break;
          }

          let finalQuery = query(collection(fs, path), ...queryConstraints);
          this.trace('query building done ------------------------ ');

          return (collectionData(finalQuery) as Observable<T[]>).pipe(
            tap((items: any) => {
              this.trace(
                'snapshot tapper ------------------ ',
                pagingAction,
                items,
                items.length
              );
              if (items.length) {
                const ps: number = +this.pageSize;
                switch (pagingAction) {
                  case 'reset':
                  case 'first':
                    this.firstItemId = items[0].id;
                    break;
                  case 'prev':
                    // in case of previous: if there are less items than page-size, it means that we reached the
                    // first page but the paging took place from an item that would usually display itself on the first page.
                    // e.g. only first 2 items returned because previous was called from item 3 with a page size of 4.
                    // This  can happen when the user uses prev all the way from the end of the list to the start.
                    // Forcefully refresh the page to first.
                    this.nextEnabled.next(true);
                    this.lastEnabled.next(true);

                    if (items.length < ps + 1) {
                      this.paginate('reset');
                      this.trace(
                        'enablePrev:',
                        this.previousEnabled,
                        'items.length',
                        items.length,
                        'ps',
                        ps,
                        pagingAction
                      );
                    }
                    break;
                  case 'next':
                    if (items.length < this.pageSize + 1) {
                      // but only if we are not already on the first page
                      if (items[0].id !== this.firstItemId) {
                        this.firstEnabled.next(true);
                        this.previousEnabled.next(true);
                        this.paginate('last');
                      }
                    }
                    this.trace(
                      'enableNext:',
                      this.nextEnabled,
                      'items.length',
                      items.length,
                      'ps',
                      ps,
                      pagingAction
                    );
                    break;
                  case 'last':
                    break;
                }

                // Check if we have a next page by the number of items.
                // If we have pagination-size + 1 results, then the next page exists.
                const nextPageExists = items.length == ps + 1
                this.lastEnabled.next(nextPageExists);
                this.nextEnabled.next(nextPageExists);

                // enablePrev if we are not at the very first element
                // enableFirst is just for convenience. Actually we could do with enablePrev
                // Todo: what happens if the first element (firstItemId) changes somewhere in between?
                const notFirstItem = items[0].id !== this.firstItemId
                this.firstEnabled.next(notFirstItem);
                this.previousEnabled.next(notFirstItem);

                // remember the anchors for moving to previous and next pages
                this.prevAnchor = items[1]; // item[1] because we have to use endBefore for previous and have to query one item extra
                this.nextAnchor = items[items.length - 1]; // last item because we have to use startAt for next and queried one item extra only for this purpose
                // this.trace('Anchor for prev', JSON.stringify(this.currentPrevAnchor.data()), "Anchor for next", JSON.stringify(this.currentNextAnchor.data()));
              } else {
                // no items were found for the new page, reset to first. But only if we did not just move to first anyway.
                if (pagingAction !== 'first' && pagingAction !== 'reset') {
                  this.paginate('reset');
                  this.trace('reset');
                }
              }
            })
          );
        } // pipe snapshotChanges
      ) // switchMap combineLatest
    ); // pipe combineLatest
  } // constructor

  private queryFirst(q: QueryConstraint[]): QueryConstraint[] {
    this.trace('first');
    return [...q, limit(this.pageSize + 1)];
  }

  private queryPrev(q: QueryConstraint[]): QueryConstraint[] {
    this.trace('Prev');
    if (this.prevAnchor) {
      this.trace('endBefore', JSON.stringify(this.prevAnchor));
      return [...q, endBefore(this.prevAnchor), limitToLast(this.pageSize + 1)];
    } else {
      return [...q, limit(this.pageSize + 1)];
    }
  }

  private queryNext(q: QueryConstraint[]): QueryConstraint[] {
    this.trace('Next');
    if (this.nextAnchor) {
      this.trace('startAt', JSON.stringify(this.nextAnchor));
      return [...q, startAt(this.nextAnchor), limit(this.pageSize + 1)];
    } else {
      return [...q, limit(this.pageSize + 1)];
    }
  }

  private queryLast(q: QueryConstraint[]): QueryConstraint[] {
    this.trace('Last');
    return [...q, limitToLast(this.pageSize)];
  }

  private queryCurrent(q: QueryConstraint[]): QueryConstraint[] {
    this.trace('Current');
    if (this.prevAnchor) {
      this.trace('startAt', JSON.stringify(this.prevAnchor));
      return [...q, startAt(this.prevAnchor), limit(this.pageSize + 1)];
    } else {
      return [...q, limit(this.pageSize + 1)];
    }
  }

  private applyQueries(q: QueryConstraint[]): QueryConstraint[] {
    if (this.queries) {
      this.queries.forEach((c) => {
        q.push(c);
        this.trace('query', c);
      });
    }
    return q;
  }

  /**
   * Returns the currently paginated items as an Observable.
   */
  public getItems(): Observable<T[]> {
    return this.items$;
  }

  /**
   * Returns whether the paginator is able to go to the first page
   * as an Observable (no effect on paginator, mostly for UI purposes).
   */
  public getFirstEnabled(): Observable<boolean> {
    return this.firstEnabled.asObservable();
  }

  /**
   * Returns whether the paginator is able to go to the last page
   * as an Observable (no effect on paginator, mostly for UI purposes).
   */
  public getLastEnabled(): Observable<boolean> {
    return this.lastEnabled.asObservable();
  }

  /**
   * Returns whether the paginator is able to go to the next page
   * as an Observable (no effect on paginator, mostly for UI purposes).
   */
  public getNextEnabled(): Observable<boolean> {
    return this.nextEnabled.asObservable();
  }

  /**
   * Returns whether the paginator is able to go to the previous page
   * as an Observable (no effect on paginator, mostly for UI purposes).
   */
  public getPreviousEnabled(): Observable<boolean> {
    return this.previousEnabled.asObservable();
  }

  public paginate(action: PaginatorActions) {
    this.paging$.next(action);
  }

  /**
   * Sets custom queries.
   *
   * The array structure is used because sorts are ordered.
   *
   * @param sort
   */
  public setQueries(queries: QueryConstraint[] | null) {
    this.queries = queries;
    this.paginate('reset');
  }

  public setPageSize(pageSize: number) {
    this.pageSize = pageSize;
    this.paginate('current');
  }

  private trace(...items: any) {
    if (this.debug) {
      console.log(...items);
    }
  }

  /**
   * call to resume loading documents in paginator.
   * can also be used to reload the current page if joined data was updated by simply calling resume() again.
   */
  public resume() {
    this.setLoadingAndRefresh(false);
  }

  /**
   * call to stall loading documents in paginator.
   */
  public stall() {
    this.setLoadingAndRefresh(true);
  }

  public first() {
    this.paginate('first');
  }

  public last() {
    this.paginate('last');
  }

  public next() {
    this.paginate('next');
  }

  public previous() {
    this.paginate('prev');
  }

  public prev() {
    this.previous();
  }

  private setLoadingAndRefresh(loading: boolean) {
    if (this.stalled) {
      this.trace('loading set to ', loading);
      this.stalled = loading;
      if (!loading) {
        this.paginate('reset');
      }
    } else {
      this.stalled = loading;
      if (!loading) {
        this.paginate('current');
      }
    }
  }
}
