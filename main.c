#include <assert.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/time.h>
#include "ketopt.h"
#include "kalloc.h"
#include "miniwfa.h"
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
	mwf_opt_t opt;
	int c, use_kalloc = 1, approx = 0;
	double t;
	void *km = 0;

	mwf_opt_init(&opt);
	while ((c = ketopt(&o, argc, argv, 1, "cKdep:al:n:r", 0)) >= 0) {
		if (o.opt == 'K') use_kalloc = !use_kalloc;
		else if (o.opt == 'c') opt.flag |= MWF_F_CIGAR;
		else if (o.opt == 'd') opt.flag |= MWF_F_DEBUG;
		else if (o.opt == 'p') opt.flag |= MWF_F_CIGAR, opt.step = atoi(o.arg);
		else if (o.opt == 'a') opt.o2 = opt.o1, opt.e2 = opt.e1;
		else if (o.opt == 'e') opt.x = 1, opt.o1 = opt.o2 = 0, opt.e1 = opt.e2 = 1;
		else if (o.opt == 'l') opt.max_lag = atoi(o.arg);
		else if (o.opt == 'n') opt.max_width = atoi(o.arg);
		else if (o.opt == 'r') approx = 1;
		else if (1) {
			fprintf(stderr, "ERROR: unknown option\n");
			return 1;
		}
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: test-mwf [options] <in1.fa> <in2.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -c       generate CIGAR\n");
		fprintf(stderr, "  -p INT   step size (force -c; 0 to disable) [%d]\n", opt.step);
		fprintf(stderr, "  -a       mimic affine gap\n");
		fprintf(stderr, "  -e       mimic edit distance\n");
		fprintf(stderr, "  -K       disable the kalloc allocator\n");
		//fprintf(stderr, "  -n INT   start to apply heuristics if the size of WF is over INT [%d]\n", opt.max_width);
		//fprintf(stderr, "  -l INT   max lag [%d]\n", opt.max_lag);
		return 1;
	}

	fp1 = gzopen(argv[o.ind+0], "r");
	fp2 = gzopen(argv[o.ind+1], "r");
	assert(fp1 && fp2);
	ks1 = kseq_init(fp1);
	ks2 = kseq_init(fp2);

	t = cputime();
	while (kseq_read(ks1) >= 0 && kseq_read(ks2) >= 0) {
		mwf_rst_t rst;
		km = use_kalloc? km_init() : 0;
		if (approx) {
			rst.s = mg_fastcmp(km, ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, 13, 1);
			rst.cigar = 0;
		} else {
			mwf_wfa(km, &opt, ks1->seq.l, ks1->seq.s, ks2->seq.l, ks2->seq.s, &rst);
		}
		if (opt.flag & MWF_F_CIGAR) mwf_assert_cigar(&opt, rst.n_cigar, rst.cigar, ks1->seq.l, ks2->seq.l, rst.s);
		printf("%s\t%ld\t0\t%ld\t+\t%s\t%ld\t0\t%ld\t%d", ks1->name.s, ks1->seq.l, ks1->seq.l, ks2->name.s, ks2->seq.l, ks2->seq.l, rst.s);
		if (opt.flag & MWF_F_CIGAR) {
			int32_t i;
			putchar('\t');
			for (i = 0; i < rst.n_cigar; ++i)
				printf("%d%c", rst.cigar[i]>>4, "MIDNSHP=XBid"[rst.cigar[i]&0xf]);
		}
		putchar('\n');
		fflush(stdout);
		kfree(km, rst.cigar);
		if (use_kalloc) km_destroy(km);
		fprintf(stderr, "T\t%s\t%s\t%.3f\n", ks1->name.s, ks2->name.s, cputime() - t);
		t = cputime();
	}

	kseq_destroy(ks1);
	kseq_destroy(ks2);
	gzclose(fp1);
	gzclose(fp2);
	return 0;
}
