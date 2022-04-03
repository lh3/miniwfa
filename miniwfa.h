#ifndef MINIWFA_H
#define MINIWFA_H

#include <stdint.h>

#define BWF_F_CIGAR 0x1
#define BWF_F_KMDBG 0x100

typedef struct {
	int32_t x, o1, e1, o2, e2;
	int32_t flag;
} mwf_opt_t;

void mwf_opt_init(mwf_opt_t *opt);
int32_t mwf_wfa_score(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs);

#endif
