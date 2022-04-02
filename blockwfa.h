#ifndef BLOCKWFA_H
#define BLOCKWFA_H

#include <stdint.h>

#define BWF_F_CIGAR 0x1

typedef struct {
	int32_t x, o1, e1, o2, e2;
	int32_t flag;
} bwf_opt_t;

void bwf_opt_init(bwf_opt_t *opt);
int32_t bwf_wfa_score(void *km, const bwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs);

#endif
