#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "ketopt.h"
#define USE_BIWFA
#include "wavefront/wavefront_align.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

int main(int argc, char *argv[])
{
	gzFile fp1, fp2;
	kseq_t *ks1, *ks2;
	ketopt_t o = KETOPT_INIT;
	int c, cigar = 0, mem_mode = 0, affine = 0;
	double t;

	while ((c = ketopt(&o, argc, argv, 1, "acm:", 0)) >= 0) {
		if (c == 'c') cigar = 1;
		else if (c == 'm') mem_mode = atoi(o.arg);
		else if (c == 'a') affine = 1;
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-wfa [-a] [-c] [-m 0|1|2|3] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Note: -m0 for linear-space\n");
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);

	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
	attributes.heuristic.strategy = wf_heuristic_none;
	attributes.alignment_scope = cigar? compute_alignment : compute_score;
	if (affine) {
		attributes.distance_metric = gap_affine;
		attributes.affine_penalties.mismatch = 4;       // X > 0
		attributes.affine_penalties.gap_opening = 4;   // O1 >= 0
		attributes.affine_penalties.gap_extension = 2; // E1 > 0
	} else {
		attributes.distance_metric = gap_affine_2p;
		attributes.affine2p_penalties.mismatch = 4;       // X > 0
		attributes.affine2p_penalties.gap_opening1 = 4;   // O1 >= 0
		attributes.affine2p_penalties.gap_extension1 = 2; // E1 > 0
		attributes.affine2p_penalties.gap_opening2 = 15;  // O2 >= 0
		attributes.affine2p_penalties.gap_extension2 = 1; // E2 > 0
	}
#ifdef USE_BIWFA
	if (cigar == 0) mem_mode = 3; // otherwise BiWFA segfaults
	if (mem_mode == 0) {
		if (cigar == 0) fprintf(stderr, "WARNING: apply -c in the linear-memory mode\n");
		cigar = 1;
		attributes.bidirectional_alignment = true;
	} else {
		attributes.memory_mode = mem_mode <= 1? wavefront_memory_low : mem_mode == 2? wavefront_memory_med : wavefront_memory_high;
	}
#else
	attributes.memory_mode = mem_mode <= 1? wavefront_memory_low : mem_mode == 2? wavefront_memory_med : wavefront_memory_high;
#endif

	t = cputime();
	wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
	while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
		wavefront_align(wf_aligner, ks2->seq.s, ks2->seq.l, ks1->seq.s, ks1->seq.l);
		printf("%s\t%ld\t0\t%ld\t+\t%s\t%ld\t0\t%ld\t%d", ks1->name.s, ks1->seq.l, ks1->seq.l, ks2->name.s, ks2->seq.l, ks2->seq.l, -wf_aligner->cigar.score);
		if (cigar) {
			putchar('\t');
			cigar_print(stdout, &wf_aligner->cigar, 1);
		}
		putchar('\n');
		fprintf(stderr, "T\t%s\t%s\t%.3f\n", ks1->name.s, ks2->name.s, cputime() - t);
		t = cputime();
	}
	wavefront_aligner_delete(wf_aligner);

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
