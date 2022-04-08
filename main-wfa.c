#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include "ketopt.h"
#include "wavefront/wavefront_align.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, cigar = 0, mem_mode = 2;

	while ((c = ketopt(&o, argc, argv, 1, "cm:", 0)) >= 0) {
		if (c == 'c') cigar = 1;
		else if (c == 'm') mem_mode = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: wfa-test <in1.fa> <in2.fa>\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);

	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
	attributes.distance_metric = gap_affine_2p;
	attributes.affine2p_penalties.mismatch = 4;       // X > 0
	attributes.affine2p_penalties.gap_opening1 = 4;   // O1 >= 0
	attributes.affine2p_penalties.gap_extension1 = 2; // E1 > 0
	attributes.affine2p_penalties.gap_opening2 = 24;  // O2 >= 0
	attributes.affine2p_penalties.gap_extension2 = 1; // E2 > 0
	attributes.alignment_scope = cigar? compute_alignment : compute_score;
	attributes.memory_mode = mem_mode <= 1? wavefront_memory_low : mem_mode == 2? wavefront_memory_med : wavefront_memory_high;
	attributes.heuristic.strategy = wf_heuristic_none;

	while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
		wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
		wavefront_align(wf_aligner, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		printf("%s\t%ld\t0\t%ld\t+\t%s\t%ld\t0\t%ld\t%d", ks1->name.s, ks1->seq.l, ks1->seq.l, ks2->name.s, ks2->seq.l, ks2->seq.l, wf_aligner->align_status.score);
		putchar('\n');
	}

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}