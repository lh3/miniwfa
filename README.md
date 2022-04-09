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
([WFA][wfa-pub]) with 2-piece affine gap penalty. When reporting base alignment
for megabase-long sequences, miniwfa is sometimes a few times faster and tends
to use less memory in comparison to [WFA2-lib][wfa] and [wfalm][wfalm] in their
low-memory mode.

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

We only use two pairs of sequences for evaluation. The first pair consists of
the two haplotypes of NA19240 around the C4A/C4B gene. They are 100-150kb in
length with a penalty of 27k under the minimap2 penalty (*x*=4,
*o*<sub>1</sub>=4, *e*<sub>1</sub>=2, *o*<sub>2</sub>=24 and
*e*<sub>2</sub>=1).  The second pair consists of GRCh38 and CHM13 around MHC.
They are about 5Mb in length with a penalty of 232k. These sequences can be
found [via Zenodo][seq-zenodo].

We checked out WFA2-lib and wfalm on 2022-04-07 and compiled the code with
gcc-10.3.0 (no LTO) on a CentOS 7 server equipped with two Xeon 6230 CPUs.

|Method             |Path|Command line    |t<sub>MHC</sub> (s)|M<sub>MHC</sub> (GB)|t<sub>C4</sub> (s)|M<sub>C4</sub> (MB)|
|:------------------|:---|:---------------|------------------:|-------------------:|-----------------:|------------------:|
|miniwfa score-only |N   |test-mwf        |414   |0.5    |3.1   |33   |
|miniwfa high-mem   |Y   |test-mwf -c     |425   |51.6   |3.7   |736  |
|miniwfa low-mem    |Y   |test-mwf -cp5000|646   |6.2    |6.1   |266  |
|wfa2-lib score-only|N   |test-wfa        |349   |0.6    |2.3   |55   |
|wfa2-lib high-mem  |Y   |test-wfa -cm3   |      |~1000  |8.9   |14332|
|wfa2-lib med-mem   |Y   |test-wfa -cm2   |2396  |34.4   |17.6  |1173 |
|wfa2-lib low-mem   |Y   |test-wfa -cm1   |3310  |23.6   |18.5  |888  |
|wfalm high-mem     |Y   |test-wfalm -m3  |      |~1000  |15.2  |13883
|wfalm low-mem      |Y   |test-wfalm -m2  |2476  |38.1   |25.7  |1241 |
|wfalm recursive    |Y   |test-wfalm -m1  |6846  |3.1    |56.1  |300  |

When only calculating the alignment score, WFA2-lib is the fastest, probably
due to its better engineering. When reporting the alignment path, miniwfa is
the fastest. The recursive algorithm in wfalm uses the least memory but it is
an order of magnitude slower. At present, WFA2-lib and wfalm use 20 bytes per
traceback entry. I expect them to use much less memory and spend less time on
memory allocation when they adopt the 1-byte-per-entry representation.

[wfa-pub]: https://pubmed.ncbi.nlm.nih.gov/32915952/
[wfa]: https://github.com/smarco/WFA2-lib
[wfalm]: https://github.com/jeizenga/wfalm
[seq-zenodo]: https://zenodo.org/record/6056061
[auto-vec]: https://en.wikipedia.org/wiki/Automatic_vectorization
