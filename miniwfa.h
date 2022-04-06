#ifndef MINIWFA_H
#define MINIWFA_H

#include <stdint.h>

#define MWF_F_CIGAR    0x1
#define MWF_F_DEBUG    0x100

typedef struct {
	int32_t flag;
	int32_t x, o1, e1, o2, e2;
	int32_t step;
} mwf_opt_t;

typedef struct {
	int32_t s, n_cigar;
	uint32_t *cigar;
} mwf_rst_t;

#ifdef __cplusplus
extern "C" {
#endif

void mwf_opt_init(mwf_opt_t *opt);
void mwf_wfa(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r);

// These functions are in "mwf-dbg.c". For debugging only.
int32_t mwf_cigar2score(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t *tl, int32_t *ql);
void mwf_assert_cigar(const mwf_opt_t *opt, int32_t n_cigar, const uint32_t *cigar, int32_t tl0, int32_t ql0, int32_t s0);

#ifdef __cplusplus
}
#endif

#endif
