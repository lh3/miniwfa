#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "miniwfa.h"
#include "kalloc.h"

#if defined(__clang__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("clang loop vectorize(enable)")
#elif defined(__GNUC__)
  #define PRAGMA_LOOP_VECTORIZE _Pragma("GCC ivdep")
#else
  #define PRAGMA_LOOP_VECTORIZE _Pragma("ivdep")
#endif

#define WF_NEG_INF (-0x40000000)

/*
 * Default setting
 */
void mwf_opt_init(mwf_opt_t *opt)
{
	memset(opt, 0, sizeof(*opt));
	opt->x  = 2;
	opt->o1 = 4, opt->e1 = 2;
	opt->o2 = 24, opt->e2 = 1;
}

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

/*
 * Structs and simple functions for traceback
 */
typedef struct {
	int32_t lo, hi;
	uint8_t *x;
} wf_tb1_t;

typedef struct {
	int32_t m, n;
	wf_tb1_t *a;
} wf_tb_t;

static wf_tb1_t *wf_tb_add(void *km, wf_tb_t *tb, int32_t lo, int32_t hi)
{
	wf_tb1_t *p;
	if (tb->n == tb->m) {
		tb->m += (tb->m>>1) + 4;
		tb->a = Krealloc(km, wf_tb1_t, tb->a, tb->m);
	}
	p = &tb->a[tb->n++];
	p->lo = lo, p->hi = hi;
	p->x = Kcalloc(km, uint8_t, hi - lo + 1);
	return p;
}

typedef struct {
	int32_t m, n;
	uint32_t *cigar;
} wf_cigar_t;

static void wf_cigar_push1(void *km, wf_cigar_t *c, int32_t op, int32_t len)
{
	if (c->n && op == (c->cigar[c->n-1]&0xf)) {
		c->cigar[c->n-1] += len<<4;
	} else {
		if (c->n == c->m) {
			c->m = c->m + (c->m>>1) + 4;
			c->cigar = Krealloc(km, uint32_t, c->cigar, c->m);
		}
		c->cigar[c->n++] = len<<4 | op;
	}
}

int32_t mwf_cigar2score(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar)
{
	int32_t k, s;
	for (k = 0, s = 0; k < n_cigar; ++k) {
		int32_t op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 1 || op == 2) {
			int32_t s1 = opt->o1 + len * opt->e1, s2 = opt->o2 + len * opt->e2;
			s += s1 < s2? s1 : s2;
		} else if (op == 8) s += opt->x;
	}
	return s;
}

/*
 * Core algorithm
 */
typedef struct {
	int32_t lo, hi;
	int32_t *mem, *H, *E1, *E2, *F1, *F2;
} wf_slice_t;

typedef struct {
	int32_t s, top, n, max_pen;
	wf_slice_t *a;
} wf_stripe_t;

void wf_stripe_add(void *km, wf_stripe_t *wf, int32_t lo, int32_t hi)
{
	int32_t i, n, m1 = wf->max_pen + 1, m2 = m1 * 2;
	wf_slice_t *f;
	++wf->s;
	++wf->top;
	if (wf->top == wf->n) wf->top = 0;
	f = &wf->a[wf->top];
	f->lo = lo, f->hi = hi;
	n = hi - lo + 1;
	kfree(km, f->mem);
	f->mem = Kmalloc(km, int32_t, 5 * (n + m2));
	f->H = f->mem + m1;
	f->E1 = f->H  + n + m2;
	f->F1 = f->E1 + n + m2;
	f->E2 = f->F1 + n + m2;
	f->F2 = f->E2 + n + m2;
	for (i = -m1; i < 0; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	for (i = n; i < n + m1; ++i)
		f->H[i] = f->E1[i] = f->E2[i] = f->F1[i] = f->F2[i] = WF_NEG_INF;
	f->H -= lo, f->E1 -= lo, f->E2 -= lo, f->F1 -= lo, f->F2 -= lo; // such that f->H[lo] points to 0
}

static wf_stripe_t *wf_stripe_init(void *km, int32_t max_pen)
{
	int32_t i;
	wf_stripe_t *wf;
	wf = Kcalloc(km, wf_stripe_t, 1);
	wf->max_pen = max_pen;
	wf->n = max_pen + 1;
	wf->a = Kcalloc(km, wf_slice_t, wf->n);
	for (i = 0; i < wf->n; ++i) {
		wf_slice_t *f;
		wf_stripe_add(km, wf, 0, 0);
		f = &wf->a[wf->top];
		f->H[0] = f->E1[0] = f->E2[0] = f->F1[0] = f->F2[0] = WF_NEG_INF;
	}
	wf->s = 0;
	wf->a[wf->top].H[0] = -1;
	return wf;
}

static void wf_stripe_destroy(void *km, wf_stripe_t *wf)
{
	int32_t i;
	for (i = 0; i < wf->n; ++i)
		kfree(km, wf->a[i].mem);
	kfree(km, wf);
}

static inline wf_slice_t *wf_stripe_get(wf_stripe_t *wf, int32_t x)
{
	int32_t y = wf->top - x;
	if (y < 0) y += wf->n;
	return &wf->a[y];
}

#define wf_max(a, b) ((a) >= (b)? (a) : (b))

static void wf_next_basic(void *km, void *km_tb, const mwf_opt_t *opt, wf_stripe_t *wf, wf_tb_t *tb, int32_t lo, int32_t hi)
{
	int32_t *H, *E1, *E2, *F1, *F2, d;
	const int32_t *pHx, *pHo1, *pHo2, *pE1, *pE2, *pF1, *pF2;
	const wf_slice_t *fx, *fo1, *fo2, *fe1, *fe2;
	wf_slice_t *ft;
	wf_stripe_add(km, wf, lo, hi);
	ft  = &wf->a[wf->top];
	fx  = wf_stripe_get(wf, opt->x);
	fo1 = wf_stripe_get(wf, opt->o1 + opt->e1);
	fo2 = wf_stripe_get(wf, opt->o2 + opt->e2);
	fe1 = wf_stripe_get(wf, opt->e1);
	fe2 = wf_stripe_get(wf, opt->e2);
	pHx = fx->H, pHo1 = fo1->H, pHo2 = fo2->H, pE1 = fe1->E1, pE2 = fe2->E2, pF1 = fe1->F1, pF2 = fe2->F2;
	H = ft->H, E1 = ft->E1, E2 = ft->E2, F1 = ft->F1, F2 = ft->F2;
	if (!tb) {
		PRAGMA_LOOP_VECTORIZE
		for (d = lo; d <= hi; ++d) {
			int32_t h, f, e;
			E1[d] = wf_max(pHo1[d-1], pE1[d-1]);
			E2[d] = wf_max(pHo2[d-1], pE2[d-1]);
			e = wf_max(E1[d], E2[d]);
			F1[d] = wf_max(pHo1[d+1], pF1[d+1]) + 1;
			F2[d] = wf_max(pHo2[d+1], pF2[d+1]) + 1;
			f = wf_max(F1[d], F2[d]);
			h = wf_max(e, f);
			H[d] = wf_max(pHx[d] + 1, h);
			// if (H[d] >= -1) fprintf(stderr, "s=%d, d=%d, k=%d, (%d,%d)\n", wf->s, d, H[d], E1[d], F1[d]);
		}
	} else {
		uint8_t *ax;
		ax = wf_tb_add(km_tb, tb, lo, hi)->x - lo;
		PRAGMA_LOOP_VECTORIZE
		for (d = lo; d <= hi; ++d) {
			int32_t h, f, e;
			uint8_t x = 0, ze, zf, z;
			x |= pHo1[d-1] >= pE1[d-1]? 0 : 0x08;
			E1[d] = wf_max(pHo1[d-1], pE1[d-1]);
			x |= pHo2[d-1] >= pE2[d-1]? 0 : 0x20;
			E2[d] = wf_max(pHo2[d-1], pE2[d-1]);
			ze = E1[d] >= E2[d]? 1 : 3;
			e = wf_max(E1[d], E2[d]);
			x |= pHo1[d+1] >= pF1[d+1]? 0 : 0x10;
			F1[d] = wf_max(pHo1[d+1], pF1[d+1]) + 1;
			x |= pHo2[d+1] >= pF2[d+1]? 0 : 0x40;
			F2[d] = wf_max(pHo2[d+1], pF2[d+1]) + 1;
			zf = F1[d] >= F2[d]? 2 : 4;
			f = wf_max(F1[d], F2[d]);
			z = e >= f? ze : zf;
			h = wf_max(e, f);
			z = pHx[d] + 1 >= h? 0 : z;
			H[d] = wf_max(pHx[d] + 1, h);
			ax[d] = x | z;
		}
	}
}

static uint32_t *wf_traceback(void *km, const mwf_opt_t *opt, wf_tb_t *tb, int32_t t_end, const char *ts, int32_t q_end, const char *qs, int32_t *n_cigar)
{ // adapted from ksw_backtrack() from ksw2
	wf_cigar_t cigar = {0,0,0};
	int32_t i = q_end, k = t_end, s = tb->n - 1, last = 0;
	while (i >= 0 && k >= 0) {
		int32_t k0 = k, j, x, state, ext;
		if (last == 0) { // if the previous state is 0, check exact matches
			while (i >= 0 && k >= 0 && qs[i] == ts[k])
				--i, --k;
			if (k0 - k > 0)
				wf_cigar_push1(km, &cigar, 7, k0 - k);
			if (i < 0 || k < 0) break;
		}
		assert(s >= 0);
		j = i - k - tb->a[s].lo;
		assert(j <= tb->a[s].hi - tb->a[s].lo);
		x = tb->a[s].x[j];
		state = x & 0x7;
		ext = state > 0? x>>(state+2)&1 : 0; // whether an extension
		if (state == 0) {
			wf_cigar_push1(km, &cigar, 8, 1);
			--i, --k, s -= opt->x;
		} else if (state == 1) {
			wf_cigar_push1(km, &cigar, 1, 1);
			--i, s -= ext? opt->e1 : opt->o1 + opt->e1;
		} else if (state == 3) {
			wf_cigar_push1(km, &cigar, 1, 1);
			--i, s -= ext? opt->e2 : opt->o2 + opt->e2;
		} else if (state == 2) {
			wf_cigar_push1(km, &cigar, 2, 1);
			--k, s -= ext? opt->e1 : opt->o1 + opt->e1;
		} else if (state == 4) {
			wf_cigar_push1(km, &cigar, 2, 1);
			--k, s -= ext? opt->e2 : opt->o2 + opt->e2;
		} else abort();
		last = state > 0 && ext? state : 0;
	}
	if (i >= 0) wf_cigar_push1(km, &cigar, 1, i + 1);
	else if (k >= 0) wf_cigar_push1(km, &cigar, 2, k + 1);
	for (i = 0; i < cigar.n>>1; ++i) { // reverse to the input order
		uint32_t t = cigar.cigar[i];
		cigar.cigar[i] = cigar.cigar[cigar.n - i - 1];
		cigar.cigar[cigar.n - i - 1] = t;
	}
	*n_cigar = cigar.n;
	return cigar.cigar;
}

void mwf_wfa_basic(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r)
{
	int32_t lo = 0, hi = 0, is_tb = !!(opt->flag&BWF_F_CIGAR);
	int32_t max_pen = opt->x;
	wf_stripe_t *wf;
	wf_tb_t tb = {0,0,0};
	void *km_tb;

	km_tb = is_tb? km_init2(km, 8000000) : 0; // this is slightly smaller than the kalloc block size
	max_pen = max_pen > opt->o1 + opt->e1? max_pen : opt->o1 + opt->e1;
	max_pen = max_pen > opt->o2 + opt->e2? max_pen : opt->o2 + opt->e2;
	wf = wf_stripe_init(km, max_pen);

	while (1) {
		int32_t d, *H = wf->a[wf->top].H;
		for (d = lo; d <= hi; ++d) {
			int32_t k;
			if (H[d] < -1 || d + H[d] < -1) continue;
			k = wf_extend1(tl, ts, ql, qs, H[d], d);
			if (k == tl - 1 && d + k == ql - 1)
				break;
			H[d] = k;
		}
		if (d <= hi) break;
		if (lo > -tl) --lo;
		if (hi < ql)  ++hi;
		wf_next_basic(km, km_tb, opt, wf, is_tb? &tb : 0, lo, hi);
	}
	r->s = wf->s;
	if (km && (opt->flag&BWF_F_DEBUG)) {
		km_stat_t st;
		km_stat(km, &st);
		fprintf(stderr, "cap=%ld, avail=%ld, n_blks=%ld\n", st.capacity, st.available, st.n_blocks);
	}
	if (is_tb) {
		r->cigar = wf_traceback(km, opt, &tb, tl-1, ts, ql-1, qs, &r->n_cigar);
		if (opt->flag&BWF_F_DEBUG) fprintf(stderr, "s_ori=%d, s_recal=%d\n", r->s, mwf_cigar2score(opt, r->n_cigar, r->cigar));
		km_destroy(km_tb);
	}
	wf_stripe_destroy(km, wf);
}
