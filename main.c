#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "kalloc.h"
#include "blockwfa.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	bwf_opt_t opt;
	int c, use_kalloc = 1;
	void *km = 0;

	bwf_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "K", 0)) >= 0) {
		if (o.opt == 'K') use_kalloc = !use_kalloc;
		else if (1) {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: wfa-test [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);

	km = use_kalloc? km_init() : 0;
	while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
		int32_t s;
		s = bwf_wfa_score(km, &opt, ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s);
		printf("%s\t%ld\t0\t%ld\t+\t%s\t%ld\t0\t%ld\t%d", ks1->name.s, ks1->seq.l, ks1->seq.l, ks2->name.s, ks2->seq.l, ks2->seq.l, s);
		putchar('\n');
	}
	if (use_kalloc) km_destroy(km);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
