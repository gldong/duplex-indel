/* This file is duplicated from https://github.com/lh3/lianti (21a15c8)
	and originally written by Heng Li. We thank Heng Li for allowing us 
	to use this script for the adaptation of indel calling in this repo. */

/* The MIT License

   Copyright (c) 2008, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

// This is a simplified version of ksort.h

#ifndef AC_KSORT_H
#define AC_KSORT_H

#include <stdlib.h>
#include <string.h>

typedef struct {
	void *left, *right;
	int depth;
} ks_isort_stack_t;

#define KSORT_SWAP(type_t, a, b) { register type_t t=(a); (a)=(b); (b)=t; }

#define KSORT_INIT(name, type_t, __sort_lt)								\
	static inline void __ks_insertsort_##name(type_t *s, type_t *t)		\
	{																	\
		type_t *i, *j, swap_tmp;										\
		for (i = s + 1; i < t; ++i)										\
			for (j = i; j > s && __sort_lt(*j, *(j-1)); --j) {			\
				swap_tmp = *j; *j = *(j-1); *(j-1) = swap_tmp;			\
			}															\
	}																	\
	void ks_combsort_##name(size_t n, type_t a[])						\
	{																	\
		const double shrink_factor = 1.2473309501039786540366528676643; \
		int do_swap;													\
		size_t gap = n;													\
		type_t tmp, *i, *j;												\
		do {															\
			if (gap > 2) {												\
				gap = (size_t)(gap / shrink_factor);					\
				if (gap == 9 || gap == 10) gap = 11;					\
			}															\
			do_swap = 0;												\
			for (i = a; i < a + n - gap; ++i) {							\
				j = i + gap;											\
				if (__sort_lt(*j, *i)) {								\
					tmp = *i; *i = *j; *j = tmp;						\
					do_swap = 1;										\
				}														\
			}															\
		} while (do_swap || gap > 2);									\
		if (gap != 1) __ks_insertsort_##name(a, a + n);					\
	}																	\
	void ks_introsort_##name(size_t n, type_t a[])						\
	{																	\
		int d;															\
		ks_isort_stack_t *top, *stack;									\
		type_t rp, swap_tmp;											\
		type_t *s, *t, *i, *j, *k;										\
																		\
		if (n < 1) return;												\
		else if (n == 2) {												\
			if (__sort_lt(a[1], a[0])) { swap_tmp = a[0]; a[0] = a[1]; a[1] = swap_tmp; } \
			return;														\
		}																\
		for (d = 2; 1ul<<d < n; ++d);									\
		stack = (ks_isort_stack_t*)malloc(sizeof(ks_isort_stack_t) * ((sizeof(size_t)*d)+2)); \
		top = stack; s = a; t = a + (n-1); d <<= 1;						\
		while (1) {														\
			if (s < t) {												\
				if (--d == 0) {											\
					ks_combsort_##name(t - s + 1, s);					\
					t = s;												\
					continue;											\
				}														\
				i = s; j = t; k = i + ((j-i)>>1) + 1;					\
				if (__sort_lt(*k, *i)) {								\
					if (__sort_lt(*k, *j)) k = j;						\
				} else k = __sort_lt(*j, *i)? i : j;					\
				rp = *k;												\
				if (k != t) { swap_tmp = *k; *k = *t; *t = swap_tmp; }	\
				for (;;) {												\
					do ++i; while (__sort_lt(*i, rp));					\
					do --j; while (i <= j && __sort_lt(rp, *j));		\
					if (j <= i) break;									\
					swap_tmp = *i; *i = *j; *j = swap_tmp;				\
				}														\
				swap_tmp = *i; *i = *t; *t = swap_tmp;					\
				if (i-s > t-i) {										\
					if (i-s > 16) { top->left = s; top->right = i-1; top->depth = d; ++top; } \
					s = t-i > 16? i+1 : t;								\
				} else {												\
					if (t-i > 16) { top->left = i+1; top->right = t; top->depth = d; ++top; } \
					t = i-s > 16? i-1 : s;								\
				}														\
			} else {													\
				if (top == stack) {										\
					free(stack);										\
					__ks_insertsort_##name(a, a+n);						\
					return;												\
				} else { --top; s = (type_t*)top->left; t = (type_t*)top->right; d = top->depth; } \
			}															\
		}																\
	} \
	type_t ks_ksmall_##name(size_t n, type_t arr[], size_t kk)			\
	{																	\
		type_t *low, *high, *k, *ll, *hh, *mid;							\
		low = arr; high = arr + n - 1; k = arr + kk;					\
		for (;;) {														\
			if (high <= low) return *k;									\
			if (high == low + 1) {										\
				if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \
				return *k;												\
			}															\
			mid = low + (high - low) / 2;								\
			if (__sort_lt(*high, *mid)) KSORT_SWAP(type_t, *mid, *high); \
			if (__sort_lt(*high, *low)) KSORT_SWAP(type_t, *low, *high); \
			if (__sort_lt(*low, *mid)) KSORT_SWAP(type_t, *mid, *low);	\
			KSORT_SWAP(type_t, *mid, *(low+1));							\
			ll = low + 1; hh = high;									\
			for (;;) {													\
				do ++ll; while (__sort_lt(*ll, *low));					\
				do --hh; while (__sort_lt(*low, *hh));					\
				if (hh < ll) break;										\
				KSORT_SWAP(type_t, *ll, *hh);							\
			}															\
			KSORT_SWAP(type_t, *low, *hh);								\
			if (hh <= k) low = ll;										\
			if (hh >= k) high = hh - 1;									\
		}																\
	}																	\

#define ks_introsort(name, n, a) ks_introsort_##name(n, a)
#define ks_ksmall(name, n, a, k) ks_ksmall_##name(n, a, k)

#define ks_lt_generic(a, b) ((a) < (b))
#define ks_lt_str(a, b) (strcmp((a), (b)) < 0)

typedef const char *ksstr_t;

#define KSORT_INIT_GENERIC(type_t) KSORT_INIT(type_t, type_t, ks_lt_generic)
#define KSORT_INIT_STR KSORT_INIT(str, ksstr_t, ks_lt_str)

#endif
