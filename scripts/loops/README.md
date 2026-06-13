## Loop Analysis

### Loop Annotation

We call loops with `cooltools dots` using the ICE-balanced merged contact matrices as input with default parameters @ 25,10,5Kb resolutions.
For all subsequent analyses we only use loops with a `qvalue < 0.1` with the `donut` kernel.
Columns are:

### Loop Comparison

We evaluate wheter a loop is replicated or not between two conditions using the `IDR2D` R package.
We rank loops by their `-log10(q.value)`, and we resolve ambiguous overlaps by the picking the most significant (`'ambiguity_resolution_method='value'`) and we allow +/- 5 bins offset when matching loops to compare (`max_gap = 5 * resolution`).
Columns are:

### Loop Valency

How many loops is an individual loop anchor a part of, described in [this paper](https://www.nature.com/articles/s41592-021-01248-7#Sec10).
Columns are:

### Loop Nesting

We also want to analyze how "nested" regions of the genome are within loops. Since loops can overlap, some genomic regions are "within" multiple loops.
If a region is nested within several loops, is it's genomic accessability more important/tighly regulated/constrained? 
Are heavily nested regions enriched for Differential HiC signals or DEGs or DEG-linked loci (e.g. enahcners) or functional elements e.g. TFBS?
So nesting can be simply computed per genomic bin as the number of loops overlapping that bin, and we can collapse contiguous sets of bins (which must be overlapped by the same set of loops) into segments.
We can further extend this by 1) computing summary stats (not just counts) across all overlapping loops per bin e.g. mean loop log10(qvalue) and 2) stratify loops by differential status i.e. nesting of only loops only detected in CTCF.BIALLELIC when comparing CTCF.BIALLELIC vs CTCF.WT vs nesting of commonly detected loops (i.e. loops w/ IDR < 0.1).
Columns are:

#### Rolling Nesting Correlation 

We can also calculate a rolling correlation of nesting value between two comparisons across genomic segments coverted by a loop in either condition. So a "genomic segment" is largest set of adjacent bins with >=1 overlapping loop in either condition. For each segment we compute the correlation of the nesting value across all bins and repeat separately per condition. So this data is the segment-wise correlation test results (corr value + pvalue) for all pairs of conditions, for loops stratified by differential status and across summary stats.
Columns are:

### Generate Loop results

Commands to generate loop results
```bash
# generate loop annotations 
mamba activate r
Rscript scripts/loops/make.loop.calling.cmds.R
mamba activate cooltools
parallel -j $(nproc) --bar --eta :::: ./results/loops/all.loop.cmds.txt
# Use IDR2D to define which loops are reproducible between conditions
mamba activate r
Rscript scripts/loops/run.IDR2D.loops.R
# calculate loop valency i.e. how many loops each loop anchor is a part of
mamba activate r
Rscript scripts/loops/calculate.loop.valency.R
# calculate loop nesting i.e. for each bin, how many loops overlap that bin
mamba activate r
Rscript scripts/loops/calculate.loop.nesting.level.R 
```
