import { Directive, computed, input } from '@angular/core';
import { hlm } from '@spartan-ng/ui-core';
import type { ClassValue } from 'clsx';

export const hlmH2 =
	'scroll-m-20 border-border border-b pb-1.5 text-lg font-semibold tracking-tight transition-colors [&:not(:first-child)]:pt-1.5';

@Directive({
	selector: '[hlmH2]',
	standalone: true,
	host: {
		'[class]': '_computedClass()',
	},
})
export class HlmH2Directive {
	public readonly userClass = input<ClassValue>('', { alias: 'class' });
	protected _computedClass = computed(() => hlm(hlmH2, this.userClass()));
}
