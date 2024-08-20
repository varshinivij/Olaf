import { ComponentFixture, TestBed } from '@angular/core/testing';

import { TestE2bComponent } from './test-e2b.component';

describe('TestE2bComponent', () => {
  let component: TestE2bComponent;
  let fixture: ComponentFixture<TestE2bComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      imports: [TestE2bComponent]
    })
    .compileComponents();

    fixture = TestBed.createComponent(TestE2bComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
