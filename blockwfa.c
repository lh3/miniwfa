#include <stdio.h>
#include "blockwfa.h"
#include "kalloc.h"

/*
 * Default setting
 */
void bwf_opt_init(bwf_opt_t *opt)
{
	opt->x  = 2;
	opt->o1 = 4, opt->e1 = 2;
	opt->o2 = 24, opt->e2 = 1;
}

#define WF_NEG_INF 0x40000000

// Extend a diagonal along exact matches. This is a bottleneck and could be made faster with padding.
static inline int32_t wf_extend1(int32_t tl, const char *ts, int32_t ql, const char *qs, int32_t k, int32_t d)
{
	const char *ts_, *qs_;
	uint64_t cmp = 0;
	int32_t max_k = (ql - d < tl? ql - d : tl) - 1;
	ts_ = ts + 1;
	qs_ = qs + d + 1;
	while (k + 7 < max_k) {
		uint64_t x = *(uint64_t*)(ts_ + k); // warning: unaligned memory access
		uint64_t y = *(uint64_t*)(qs_ + k);
		cmp = x ^ y;
		if (cmp == 0) k += 8;
		else break;
	}
	if (cmp)
		k += __builtin_ctzl(cmp) >> 3; // on x86, this is done via the BSR instruction: https://www.felixcloutier.com/x86/bsr
	else if (k + 7 >= max_k)
		while (k < max_k && *(ts_ + k) == *(qs_ + k)) // use this for generic CPUs. It is slightly faster than the unoptimized version
			++k;
	return k;
}

typedef struct {
	int32_t lo, hi;
	int32_t *mem, *H, *E1, *E2, *F1, *F2;
} bwf_slice_t;

typedef struct {
	int32_t s, top, n, max_pen;
	bwf_slice_t *a;
} bwf_stripe_t;

void wf_stripe_add(void *km, bwf_stripe_t *wf, int32_t lo, int32_t hi)
{
	int32_t i, n, m2 = wf->max_pen * 2;
	bwf_slice_t *f;
	++wf->s;
	++wf->top;
	if (wf->top == wf->n) wf->top = 0;
	f = &wf->a[wf->top];
	f->lo = lo, f->hi = hi;
	n = hi - lo + 1;
	kfree(km, f->mem);
	f->mem = Kmalloc(km, int32_t, 5 * (n + m2));
	f->H = f->mem + wf->max_pen;
	f->E1 = f->H  + n + m2;
	f->E2 = f->E1 + n + m2;
	f->F1 = f->E2 + n + m2;
	f->F2 = f->F1 + n + m2;
	for (i = -wf->max_pen + 1; i < 0; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	for (i = n; i < n + wf->max_pen; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	f->H -= lo, f->E1 -= lo, f->E2 -= lo, f->F1 -= lo, f->F2 -= lo; // such that f->H[lo] points to 0
}

static bwf_stripe_t *wf_stripe_init(void *km, int32_t max_pen)
{
	int32_t i;
	bwf_stripe_t *wf;
	wf = Kcalloc(km, bwf_stripe_t, 1);
	wf->max_pen = max_pen;
	wf->n = max_pen + 1;
	wf->a = Kcalloc(km, bwf_slice_t, wf->n);
	for (i = 0; i < wf->n; ++i) {
		bwf_slice_t *f;
		wf_stripe_add(km, wf, 0, 0);
		f = &wf->a[wf->top];
		f->H[0] = f->E1[0] = f->E2[0] = f->F1[0] = f->F2[0];
	}
	wf->s = 0;
	wf->a[wf->top].H[0] = -1;
	return wf;
}

static void wf_stripe_destroy(void *km, bwf_stripe_t *wf)
{
	int32_t i;
	for (i = 0; i < wf->n; ++i)
		kfree(km, wf->a[i].mem);
	kfree(km, wf);
}

static inline bwf_slice_t *wf_stripe_get(bwf_stripe_t *wf, int32_t x)
{
	int32_t y = wf->top - x;
	if (y < 0) y += wf->n;
	return &wf->a[y];
}

#define wf_max(a, b) ((a) >= (b)? (a) : (b))

static void wf_next(void *km, const bwf_opt_t *opt, bwf_stripe_t *wf, int32_t lo, int32_t hi)
{
	int32_t *pHx, *pHo, *pE1, *pE2, *pF1, *pF2;
	int32_t *H, *E1, *E2, *F1, *F2, d;
	const bwf_slice_t *fx, *fo, *fe, *ft;
	wf_stripe_add(km, wf, lo, hi);
	ft = &wf->a[wf->top];
	fx = wf_stripe_get(wf, opt->x);
	fo = wf_stripe_get(wf, opt->o1 + opt->e1);
	fe = wf_stripe_get(wf, opt->e1);
	pHx = fx->H, pHo = fo->H, pE1 = fe->E1, pE2 = fe->E2, pF1 = fe->F1, pF2 = fe->F2;
	H = ft->H, E1 = ft->E1, E2 = ft->E2, F1 = ft->F1, F2 = ft->F2;
	for (d = lo; d <= hi; ++d) {
		int32_t h, f, e;
		F1[d] = wf_max(pHo[d-1], pF1[d-1]) + 1;
		F2[d] = wf_max(pHo[d-1], pF2[d-1]) + 1;
		f = wf_max(F1[d], F2[d]);
		E1[d] = wf_max(pHo[d+1], pE1[d+1]);
		E2[d] = wf_max(pHo[d+1], pE2[d+1]);
		e = wf_max(E1[d], E2[d]);
		h = wf_max(e, f);
		H[d] = wf_max(pHx[d], h);
	}
}

int32_t bwf_wfa_score(void *km, const bwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs)
{
	int32_t s, lo = 0, hi = 0;
	int32_t max_pen = opt->x;
	bwf_stripe_t *wf;

	max_pen = max_pen > opt->o1 + opt->e1? max_pen : opt->o1 + opt->e1;
	max_pen = max_pen > opt->o2 + opt->e2? max_pen : opt->o2 + opt->e2;
	wf = wf_stripe_init(km, max_pen);

	while (1) {
		int32_t d, *H = wf->a[wf->top].H;
		for (d = lo; d <= hi; ++d) {
			int32_t k;
			k = wf_extend1(tl, ts, ql, qs, H[d], d);
			if (k == tl - 1 && d + k == ql - 1)
				break;
			H[d] = k;
		}
		if (d <= hi) break;
		if (lo > -tl) --lo;
		if (hi < ql)  ++hi;
		wf_next(km, opt, wf, lo, hi);
	}
	s = wf->s;
	wf_stripe_destroy(km, wf);
	return s;
}
