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
for megabase-long sequences, miniwfa is sometimes a few times faster and uses a
fraction of memory in comparison to [WFA2-lib][wfa] and [wfalm][wfalm] in their
low-memory mode.

## Algorithm

When reporting alignment score only, miniwfa behaves largely the same as the
original WFA. WFA2-lib is faster in this case, probably due to its better
engineering.  Miniwfa differs mainly in traceback. In the high memory mode,
miniwfa packs traceback information into 7 bits:
```txt
extD2<<6 | extI2<<5 | extD1<<4 | extI1<<3 | fromState
```
where `extD1`, for example, indicates whether we are extending the `D1`
deletion state, and `fromState` keeps the max state (5 distinct values).
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
is the maximal penalty between adjacent entries.

## Evaluation

We only use two pairs of sequences for evaluation. The first pair consists of
the two haplotypes of NA19240 around the C4A/C4B gene. They are 100-150kb in
length with a penalty of 27k. The second pair consists of GRCh38 and
CHM13 around MHC. They are about 5Mb in length with an edit distance of 232kb.
These sequences can be found [via Zenodo][seq-zenodo]. We compiled the code
with gcc-10.3.0 (no LTO) on a CentOS 7 server equipped with two Xeon 6130 CPUs.

|Method |Path|CMD             |t<sub>MHC</sub> (s)|M<sub>MHC</sub> (GB)|t<sub>C4</sub> (s)|M<sub>C4</sub> (MB)|Comment|
|:------|:---|:---------------|------------------:|-------------------:|-----------------:|------------------:|:------|
|miniwfa|N   |test-mwf        |440   |0.5    |4.8   |33   |Score only|
|miniwfa|Y   |test-mwf -c     |558   |51.6   |6.7   |736  |High-mem|
|miniwfa|Y   |test-mwf -cp5000|787   |6.2    |11.2  |266  |Low-mem|
|wfa2-lib|N  |test-wfa        |333   |0.6    |2.9   |55   |Score only|
|wfa2-lib|Y  |test-wfa -cm1   |4915  |22.4   |20.6  |888  |Low-mem|
|wfa2-lib|Y  |test-wfa -cm2   |2460  |34.4   |21.5  |1173 |Med-mem|
|wfa2-lib|Y  |test-wfa -cm3   |      |~1000  |8.6   |14332|High-mem|
|wfalm   |Y  |test-wfalm -m1  |      |       |73.9  |300  |Recursive|
|wfalm   |Y  |test-wfalm -m2  |2734  |38.1   |30.5  |1241 |Low-mem|
|wfalm   |Y  |test-wfalm -m3  |      |~1000  |15.4  |13883|High-mem|

[wfa-pub]: https://pubmed.ncbi.nlm.nih.gov/32915952/
[wfa]: https://github.com/smarco/WFA2-lib
[wfalm]: https://github.com/jeizenga/wfalm
[seq-zenodo]: https://zenodo.org/record/6056061
