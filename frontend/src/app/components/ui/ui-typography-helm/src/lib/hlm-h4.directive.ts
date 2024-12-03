import { Directive, computed, input } from '@angular/core';
import { hlm } from '@spartan-ng/ui-core';
import type { ClassValue } from 'clsx';

export const hlmH4 = 'scroll-m-20 text-md font-semibold tracking-tight pt-3';

@Directive({
	selector: '[hlmH4]',
	standalone: true,
	host: {
		'[class]': '_computedClass()',
	},
})
export class HlmH4Directive {
	public readonly userClass = input<ClassValue>('', { alias: 'class' });
	protected _computedClass = computed(() => hlm(hlmH4, this.userClass()));
}
