#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "lib/wfa_lm.hpp"
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, mem_mode = 2;

	while ((c = ketopt(&o, argc, argv, 1, "cm:", 0)) >= 0) {
		if (c == 'm') mem_mode = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-wfalm [-m 1|2|3] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Notes: -m1 for recursive; -m2 for low-mem; -m3 for full\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);

	auto wfa = wfalm::make_convex_wfaligner(4, 2, 3, 1, 15);

	while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
		std::pair<std::vector<wfalm::CIGAROp>, int32_t> rst;
		if (mem_mode == 3) rst = wfa.wavefront_align(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		else if (mem_mode == 2) rst = wfa.wavefront_align_low_mem(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		else if (mem_mode == 1) rst = wfa.wavefront_align_recursive(ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		printf("%s\t%ld\t0\t%ld\t+\t%s\t%ld\t0\t%ld\t%d", ks1->name.s, ks1->seq.l, ks1->seq.l, ks2->name.s, ks2->seq.l, ks2->seq.l, rst.second);
		putchar('\n');
	}

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
