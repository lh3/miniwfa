#include <assert.h>
#include <stdio.h>
#include "miniwfa.h"

// recompute score from CIGAR
int32_t mwf_cigar2score(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t *tl, int32_t *ql)
{
	int32_t k, s, x = 0, y = 0;
	for (k = 0, s = 0; k < n_cigar; ++k) {
		int32_t op = cigar[k]&0xf, len = cigar[k]>>4;
		if (op == 1 || op == 2) {
			int32_t s1 = opt->o1 + len * opt->e1, s2 = opt->o2 + len * opt->e2;
			s += s1 < s2? s1 : s2;
		} else if (op == 8) s += len * opt->x;
		if (op == 0 || op == 7 || op == 8) x += len, y += len;
		else if (op == 1) y += len;
		else if (op == 2) x += len;
	}
	if (tl) *tl = x;
	if (ql) *ql = y;
	return s;
}

void mwf_assert_cigar(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t tl0, int32_t ql0, int32_t s0)
{
	int32_t s, tl, ql;
	s = mwf_cigar2score(opt, n_cigar, cigar, &tl, &ql);
	assert(tl == tl0);
	assert(ql == ql0);
	if (s > s0) fprintf(stderr, "[%s] s0=%d, s=%d\n", __func__, s0, s);
}
