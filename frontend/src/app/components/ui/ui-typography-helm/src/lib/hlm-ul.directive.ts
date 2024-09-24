import { Directive, computed, input } from '@angular/core';
import { hlm } from '@spartan-ng/ui-core';
import type { ClassValue } from 'clsx';

export const hlmUl = 'py-2.5 pl-4 text-sm list-disc [&>li]:pt-0.5';

@Directive({
	selector: '[hlmUl]',
	standalone: true,
	host: {
		'[class]': '_computedClass()',
	},
})
export class HlmUlDirective {
	public readonly userClass = input<ClassValue>('', { alias: 'class' });
	protected _computedClass = computed(() => hlm(hlmUl, this.userClass()));
}
