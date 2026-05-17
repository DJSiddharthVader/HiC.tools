## Loop Analysis

### Loop Annotation

We call loops with `cooltools dots` using the ICE-balanced merged contact matrices as input with default parameters @ 25,10,5Kb resolutions.
For all subsequent analyses we only use loops with a `qvalue < 0.1` with the `donut` kernel.

### Loop Comparison

We evaluate wheter a loop is replicated or not between two conditions using the `IDR2D` R package.
We rank loops by their `-log10(q.value)`, and we resolve ambiguous overlaps by the picking the most significant (`'ambiguity_resolution_method='value'`) and we allow +/- 5 bins offset when matching loops to compare (`max_gap = 5 * resolution`).

### Loop Valency

How many loops is an individual loop anchor a part of, described in [this paper](https://www.nature.com/articles/s41592-021-01248-7#Sec10).
Columns are:
| Column Name           | Example Row1       | Example Row2       | Column Description |
| --------------------- | ------------------ | ------------------ | ------------- |
|  type               | cis                | cis                  | loops called on cis contacts | 
|  weight             | balanced           | balanced             | loops called with balanced matrix |
|  resolution         | 10000              | 10000                | resoluion loops are called at |
|  SampleID           | CTCF.iN.BIALLELIC  | CTCF.iN.BIALLELIC    | merged matrix samples were called in |
|  kernel             | donut              | donut                | cooltools kernel loops stats come from |
|  chr                | chr1               | chr1                 | self-explanatory | 
|  anchor.position    | 9720000            | 9860000              | genomic bin where loop anchor is|
|  valency            | 1                  | 6                    | how many unique loops this anchor anchors | 
|  length-min         | 140000             | 140000               |  summary stats over lengths of all loops at this anchor |
|  length-mean        | 140000             | 508333.3333333333    |  |
|  length-median      | 140000             | 485000               |  |
|  length-max         | 140000             | 950000               |  |
|  count-min          | 78                 | 20                   |  |
|  count-mean         | 78                 | 46.666666666666664   | summary stats over total contacts supporting each loop at this anchor   | 
|  count-median       | 78                 | 33.5                 |  |
|  count-max          | 78                 | 91                   |  |
|  enrichment-mean    | 28.133223724607973 | 13.065525417692793   | summary stats over enrichment of all loops at this anchor |
|  enrichment-min     | 28.133223724607973 | 4.391584583785407    |  |
|  enrichment-max     | 28.133223724607973 | 28.133223724607973   |  |
|  enrichment-median  | 28.133223724607973 | 8.688043000412927    |  |
|  log10.qval-mean    | 9.788118488133792  | 8.93324939463884     | summary stats over log(qvalues) of all loops at this anchor |   
|  log10.qval-min     | 9.788118488133792  | 3.430779068372195    |  |
|  log10.qval-max     | 9.788118488133792  | 21.09790044160207    |  |
|  log10.qval-median  | 9.788118488133792  | 6.791862325480757    |  |


### Loop Nesting

We also want to analyze how "nested" regions of the genome are within loops. Since loops can overlap fully (but not partially) some genomic regions are "within" multiple loops.
If a region (or functional element is nested within several loops, is it's genomic accessability more important? 
Are heavily nested regions enriched for Differential HiC signals or functional elements e.g. cCRES?


| Column Name           | Example Row1       | Example Row2       | Column Description |
| --------------------- | ------------------ | ------------------ | ------------- |
| method                | cooltools          | cooltools          | method to call loops |
| type                  | cis                | cis                | loops called on cis contacts |
| weight                | balanced           | balanced           | loops called with balanced matrix |
| resolution            | 10000              | 10000              | resoluion loops are called at |
| SampleID              | All.iN.WT          | All.iN.WT          | merged matrix samples were called in |
| kernel                | donut              | donut              | cooltools kernel loops stats come from |
| chr                   | chr1               | chr18              | self-explanatory |
| nest.start            | 115080000          | 57570000           | contiguous region that is nested within the exact same set of loops     |
| nest.end              | 115600000          | 57620000           |  |
| nesting.lvl           | 8                  | 3                  | how many unique loops this region is within |
| metric.min.length     | 890000             | 1320000            | summary stats over total contacts supporting each loop at this anchor   | 
| metric.mean.length    | 1412500            | 1510000            |  |
| metric.max.length     | 1760000            | 1780000            |  |
| metric.sum.length     | 11300000           | 4530000            |  |
| metric.min.count      | 22                 | 24                 |  summary stats over lengths of all loops at this anchor |
| metric.mean.count     | 30.25              | 30.333333333333332 |  |
| metric.max.count      | 55                 | 42                 |  |
| metric.sum.count      | 242                | 91                 |  |
| metric.min.enrichment | 5.7467889552097535 | 6.309215723321976  | summary stats over enrichment of all loops at this anchor |
| metric.mean.enrichment| 9.54862719604052   | 7.956159202393153  |  |
| metric.max.enrichment | 23.092135519693127 | 10.102903568115485 |  |
| metric.sum.enrichment | 76.38901756832416  | 23.868477607179457 |  |
| metric.min.log10_qval | 2.7506571069237604 | 2.5241860990303    | summary stats over log(qvalues) of all loops at this anchor | 
| metric.mean.log10_qval| 3.4835348226901517 | 4.648892171887046  |  |
| metric.max.log10_qval | 4.596627111104732  | 7.206963447063814  |  |
| metric.sum.log10_qval | 27.868278581521214 | 13.946676515661137 |  |


### Generate Loop results

Commands to generate loop results
```bash
# generate loop annotations 
mamba activate r
Rscript scripts/loops/
# Use IDR2D to define which loops are reproducible between conditions
mamba activate r
Rscript scripts/loops/run.IDR2D.loops.R
# calculate loop valency i.e. how many loops each loop anchor is a part of
mamba activate r
Rscript scripts/loops/calculate.loop.valency.R
# calculate loop nesting i.e. for each bin, how many loops overlap that bin
mamba activate r
Rscript scripts/loops/calculate.loop.nesting.level.R 
mamba activate cooltools
parallel -j $(nproc) --bar --eta :::: ./results/loops/all.loop.nesting.bedtools.cmds.txt
```
