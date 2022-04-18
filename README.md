## Getting Started
```sh
# install and test
git clone https://github.com/lh3/miniwfa
cd miniwfa && make test-mwf
./test-mwf -c test/t3-?.fa

# use as a library
cp miniwfa.{c,h} kalloc.{c,h} your_src/
```

## Introduction

Miniwfa is a reimplementation of the WaveFront Alignment algorithm
([WFA][wfa-pub]) with 2-piece gap penalty. It was faster than WFA2-lib when
miniwfa was first developed. Now the linear-space [BiWFA][biwfa] algorithm
consistently uses less memory, though the relative performance between miniwfa
and BiWFA varies with input sequences. In addition, BiWFA is only directly
applicable to global alignment.

## Algorithm

When reporting alignment score only, miniwfa behaves largely the same as the
original WFA. Miniwfa differs mainly in traceback. In the high memory mode,
miniwfa packs traceback information into 7 bits per entry:
```txt
extD2<<6 | extI2<<5 | extD1<<4 | extI1<<3 | fromState
```
where `extD1`, for example, indicates whether we are extending the `D1`
deletion state, and `fromState` keeps the max state (5 possible values).
As such, miniwfa uses 1 byte for each traceback entry. The standard WFA uses 20
bytes per entry.

Miniwfa has a low-memory mode inspired by but different from wfalm. In the
wfalm terminology, miniwfa stores a "stripe" of traceback information every *p*
cycles. It finds the positions on these stripes that the optimal alignment path
goes through. Miniwfa then applies the WFA algorithm for a second time. In this
round, it resets the bandwidth to 1 whenever it reaches a stripe. Miniwfa
approximately uses (20*qs*<sup>2</sup>/*p*+*ps*) bytes of memory. Here, *s* is
the optimal alignment penalty, *p* is the distance between stripes and
*q*=max(*x*,*o*<sub>1</sub>+*e*<sub>1</sub>,*o*<sub>2</sub>+*e*<sub>2</sub>)
is the maximal penalty between adjacent entries. The time complexity is
*O*(*n*(*s*+*p*)) where *n* is the length of the longer sequence.

## Vectorization

Miniwfa relies on [compiler vectorization][auto-vec] to parallelize three key
loops. Old compilers (e.g.  gcc-4.8.5) that do not support auto vectorization
will lead to lower performance. Some compilers (e.g. gcc-6.4.0) cannot
vectorize one of the loops.  Clang-13.1.6 and gcc-10.3.0 are known to work well
with miniwfa. Also note that gcc does not vectorize loops with -O2 and
vectorization is influenced by SIMD.  If you use miniwfa as a library, it is
recommended to enable `-msse4 -O3` when compiling with gcc.

## Evaluation

We only use three pairs of sequences for evaluation. The first pair consists of
the two haplotypes of NA19240 around the C4A/C4B gene. They are 100-150kb in
length with a penalty of 26,917 under penalty *x*=4,
*o*<sub>1</sub>=4, *e*<sub>1</sub>=2, *o*<sub>2</sub>=15 and
*e*<sub>2</sub>=1.  The second pair consists of GRCh38 and CHM13 around MHC.
They are about 5Mb in length with a penalty of 229,868. The third pair consists
of the two MHC haplotypes in HG002 with a penalty of 267,637. These sequences
can be found [via Zenodo][seq-zenodo].

We checked out BiWFA on 2022-04-16 and compiled the code with
gcc-10.3.0 (no LTO) on a CentOS 7 server equipped with two Xeon 6230 CPUs.
The table below shows the timing and peak memory for miniwfa and BiWFA in its
linear mode. The table used to include WFA2-lib and wfalm, which were removed
because BiWFA is a clear winner now.

|Method             |Command line    |t<sub>MHC</sub> (s)|M<sub>MHC</sub> (GB)|t<sub>HG002</sub>|M<sub>HG002</sub>|t<sub>C4</sub>|M<sub>C4</sub>|
|:------------------|:---------------|------------------:|-------------------:|----------------:|----------------:|-------------:|-------------:|
|miniwfa high-mem   |test-mwf -c     |385   |50.6   |533   |68.1  |3.6   |0.73 |
|miniwfa low-mem    |test-mwf -cp5000|554   |4.1    |746   |5.3   |5.4   |0.22 |
|biwfa linear       |test-wfa -cm0   |417   |0.4    |2603  |0.4   |26.4  |0.05 |

## Historical notes on WFA and related algorithms

I got the following notes wrong a few times. Please let me know if you found
the narrative is still inaccurate.

Given a pair of strings, let *n* be the length of the longer string and *d* is
their edit distance. Ukkonen first found the *O*(*nd*) algorithm to compute
edit distances in 1983 and [published it in 1985][U85a]. Myers independently
published the same algorithm [in 1986][myers86]. He additionally introduced a
few variations including a linear-space algorithm based on the [Hirschberg's
algorithm][lin-space]. Sometimes the *O*(*nd*) algorithm is also attributed to
[Landau and Vishkin (1989)][lv89]. I don't know if the two authors were aware
of Ukkonen's work at the time of submission in 1986, but they did cite
Ukkonen (1985) in their published paper.

The search for a similar algorithm for linear or affine gap penalties took
three decades. To the best of my knowledge, [Xin et al (2017)][leap] first
found a version of this algorithm but apparently they have never published it
in a peer-reviewed journal. [Marco-Sola et al (2021)][wfa-pub] were probably
unaware of Xin et al. They independently published the WFA algorithm as we know
today. They also gave a highly efficient implementation, beating all global
alignment algorithms by a large margin.

A major concern with the original WFA is its large memory consumption. [Eizenga
and Paten (2022)][EP22] implemented wfalm to reduce its peak memory. This repo
was inspired by this work. Then [Marco-Sola et al (2022)][biwfa] developed
BiWFA that finds the alignment in linear space. This is so far the best
performing algorithm.

It is worth noting that also in 1985, Ukkonen published [another
algorithm][U85b] with expected *O*(*nt*) time complexity where *t* is a given
threshold. This algorithm guarantees to find the edit distance *d* if
*d*&le;*t*. It is somewhat similar to banded alignment but distinct from his
*O*(*nd*) algorithm published in the same year.

[biwfa-pub]: https://www.biorxiv.org/content/10.1101/2022.04.14.488380v1
[wfa-pub]: https://pubmed.ncbi.nlm.nih.gov/32915952/
[biwfa]: https://github.com/smarco/BiWFA-paper
[wfa2]: https://github.com/smarco/WFA2-lib
[wfalm]: https://github.com/jeizenga/wfalm
[seq-zenodo]: https://zenodo.org/record/6056061
[auto-vec]: https://en.wikipedia.org/wiki/Automatic_vectorization

[myers86]: https://link.springer.com/article/10.1007/BF01840446
[U85a]: https://www.sciencedirect.com/science/article/pii/S0019995885800462
[U85b]: https://www.sciencedirect.com/science/article/abs/pii/0196677485900239
[edlib]: https://github.com/Martinsos/edlib
[lin-space]: https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm
[leap]: https://www.biorxiv.org/content/10.1101/133157v3
[myers-bit]: https://dl.acm.org/doi/10.1145/316542.316550
[EP22]: https://www.biorxiv.org/content/10.1101/2022.01.12.476087v1
[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
