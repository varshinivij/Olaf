import { TestBed } from '@angular/core/testing';

import { E2bService } from './e2b.service';

describe('E2bService', () => {
  let service: E2bService;

  beforeEach(() => {
    TestBed.configureTestingModule({});
    service = TestBed.inject(E2bService);
  });

  it('should be created', () => {
    expect(service).toBeTruthy();
  });
});
