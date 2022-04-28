#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include "kalloc.h"

static unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

/*
 * from ksort.h
 */
#define RS_MIN_SIZE 64
#define RS_MAX_BITS 8

#define KRADIX_SORT_INIT(name, rstype_t, rskey, sizeof_key) \
	typedef struct { \
		rstype_t *b, *e; \
	} rsbucket_##name##_t; \
	void rs_insertsort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		rstype_t *i; \
		for (i = beg + 1; i < end; ++i) \
			if (rskey(*i) < rskey(*(i - 1))) { \
				rstype_t *j, tmp = *i; \
				for (j = i; j > beg && rskey(tmp) < rskey(*(j-1)); --j) \
					*j = *(j - 1); \
				*j = tmp; \
			} \
	} \
	void rs_sort_##name(rstype_t *beg, rstype_t *end, int n_bits, int s) \
	{ \
		rstype_t *i; \
		int size = 1<<n_bits, m = size - 1; \
		rsbucket_##name##_t *k, b[1<<RS_MAX_BITS], *be = b + size; \
		assert(n_bits <= RS_MAX_BITS); \
		for (k = b; k != be; ++k) k->b = k->e = beg; \
		for (i = beg; i != end; ++i) ++b[rskey(*i)>>s&m].e; \
		for (k = b + 1; k != be; ++k) \
			k->e += (k-1)->e - beg, k->b = (k-1)->e; \
		for (k = b; k != be;) { \
			if (k->b != k->e) { \
				rsbucket_##name##_t *l; \
				if ((l = b + (rskey(*k->b)>>s&m)) != k) { \
					rstype_t tmp = *k->b, swap; \
					do { \
						swap = tmp; tmp = *l->b; *l->b++ = swap; \
						l = b + (rskey(tmp)>>s&m); \
					} while (l != k); \
					*k->b++ = tmp; \
				} else ++k->b; \
			} else ++k; \
		} \
		for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k-1)->e; \
		if (s) { \
			s = s > n_bits? s - n_bits : 0; \
			for (k = b; k != be; ++k) \
				if (k->e - k->b > RS_MIN_SIZE) rs_sort_##name(k->b, k->e, n_bits, s); \
				else if (k->e - k->b > 1) rs_insertsort_##name(k->b, k->e); \
		} \
	} \
	void radix_sort_##name(rstype_t *beg, rstype_t *end) \
	{ \
		if (end - beg <= RS_MIN_SIZE) rs_insertsort_##name(beg, end); \
		else rs_sort_##name(beg, end, RS_MAX_BITS, (sizeof_key - 1) * RS_MAX_BITS); \
	}

#define sort_key_64(a) ((a))
KRADIX_SORT_INIT(mwf64, uint64_t, sort_key_64, 8) 

/*
 * from fastcmp.c
 */
static int32_t mg_lis_64(void *km, int32_t n, const uint64_t *a, int32_t *b)
{
	int32_t i, k, L = 0, *M, *P = b;
	KMALLOC(km, M, n+1);
	for (i = 0; i < n; ++i) {
		int32_t lo = 1, hi = L, newL;
		while (lo <= hi) {
			int32_t mid = (lo + hi + 1) >> 1;
			if (a[M[mid]] < a[i]) lo = mid + 1;
			else hi = mid - 1;
		}
		newL = lo, P[i] = M[newL - 1], M[newL] = i;
		if (newL > L) L = newL;
	}
	k = M[L];
	memcpy(M, P, n * sizeof(int32_t));
	for (i = L - 1; i >= 0; --i) b[i] = k, k = M[k];
	kfree(km, M);
	return L;
}

static int32_t mg_fc_kmer(int32_t len, const char *seq, int32_t rid, int32_t k, uint64_t *a)
{
	int32_t i, l, n;
	uint64_t x, mask = (1ULL<<k*2) - 1;
	for (i = l = 0, x = 0, n = 0; i < len; ++i) {
		int32_t c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) {
			x = (x << 2 | c) & mask;
			if (++l >= k) a[n++] = (x<<1|rid) << 32 | i;
		} else l = 0, x = 0;
	}
	return n;
}

int32_t mg_fastcmp(void *km, int32_t l1, const char *s1, int32_t l2, const char *s2, int32_t k, int32_t max_occ)
{
	int32_t i, n_a, n_b, m_b, i0, mlen, n_lis, *lis;
	uint64_t *a, *b;

	if (l1 < k || l2 < k) return 0;
	assert(k >= 2 && k <= 15);

	// collect k-mers
	KMALLOC(km, a, l1 + l2);
	n_a = mg_fc_kmer(l1, s1, 0, k, a);
	n_a += mg_fc_kmer(l2, s2, 1, k, &a[n_a]);
	radix_sort_mwf64(a, a + n_a);

	// collect k-mer matches
	n_b = m_b = 0, b = 0;
	for (i0 = 0, i = 1; i <= n_a; ++i) {
		if (i == n_a || a[i0]>>33 != a[i]>>33) {
			if (i - i0 >= 2) {
				int32_t j, s, t;
				for (j = i0; j < i && (a[j]>>32&1) == 0; ++j) {}
				if (j > i0 && j < i && j - i0 <= max_occ && i - j <= max_occ) {
					for (s = i0; s < j; ++s)
						for (t = j; t < i; ++t) {
							if (n_b == m_b) KEXPAND(km, b, m_b);
							b[n_b++] = a[s]<<32 | (uint32_t)a[t];
						}
				}
			}
			i0 = i;
		}
	}
	kfree(km, a);

	// compute mlen
	radix_sort_mwf64(b, b + n_b);
	for (i = 0; i < n_b; ++i)
		b[i] = b[i]>>32 | b[i]<<32;
	KMALLOC(km, lis, n_b);
	n_lis = mg_lis_64(km, n_b, b, lis);
	for (i = 1, mlen = k; i < n_lis; ++i) {
		int32_t ll2 = (int32_t)(b[lis[i]]>>32) - (int32_t)(b[lis[i-1]]>>32);
		int32_t ll1 = (int32_t)b[lis[i]] - (int32_t)b[lis[i-1]];
		mlen += ll1 > k && ll2 > k? k : ll1 < ll2? ll1 : ll2;
	}
	kfree(km, lis);
	kfree(km, b);
	return mlen;
}
