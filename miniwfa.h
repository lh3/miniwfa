#ifndef MINIWFA_H
#define MINIWFA_H

#include <stdint.h>

#define MWF_F_CIGAR    0x1
#define MWF_F_DEBUG    0x100

typedef struct {
	int32_t x, o1, e1, o2, e2;
	int32_t flag;
} mwf_opt_t;

typedef struct {
	int32_t s, n_cigar;
	uint32_t *cigar;
} mwf_rst_t;

#ifdef __cplusplus
extern "C" {
#endif

void mwf_opt_init(mwf_opt_t *opt);
void mwf_wfa_basic(void *km, const mwf_opt_t *opt, int32_t tl, const char *ts, int32_t ql, const char *qs, mwf_rst_t *r);

#ifdef __cplusplus
}
#endif

#endif
