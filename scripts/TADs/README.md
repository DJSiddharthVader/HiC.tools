## Generating TAD Results

### TADs

We generate the TAD annotations we use 3 different programs:

1. [HiTAD](https://xiaotaowang.github.io/TADLib/domaincaller.html)
1. [cooltools insulation](https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html)
2. [TADCompare consensusTAD()](https://pubmed.ncbi.nlm.nih.gov/32211023/)

For `HiTAD` we only generate single-level TADs (not hierarchical) using default parameters. This produces a set of bin-pairs, each pair marking the start and end of the predicted TAD.

Both of these tools also output the Diamond Insulation (DI) score calcualted per genomic bin.
This can be used to score regions to compare the relative insulation of regions. 

The TADs called from all methods are stored as folows
Columns are: 
| Column Name | Example Row 1 | Example Row 2 | Column Description |
| ----------- | ------------- | ------------- | ------------------ |
| resolution       |  100000    | 100000    | resolution TADs were called at | 
| method           |  hiTAD     | hiTAD     | which TAD calling method was used| 
| TAD.params       |  NA        | NA        | method-specific params for TAD calling |
| Sample.Group     |  All.iN.WT | All.iN.WT | Samples used as input to call TADs| 
| chr              |  chr1      | chr1      | chromosome TAD is on |
| start            |  3800000   | 4500000   | start of TAD in bp |
| end              |  4500000   | 5500000   | end of TAd in bp |
| TAD.length       |  700000    | 1000000   | TAD length in bp |
| TAD.bins         |  7         | 10        | tad length in bins |
| TAD.start.score  |  1.748     | 2.929     | method specific TAD score for the start bin |
| TAD.end.score    |  3.175     | 1.835     | method specific TAD score for the end bin |
| TAD.inner.min    |  -17.04    | -25.13    | summary stats computed over all bin scores within the TAD |
| TAD.inner.q25    |  -7.916    | -2.44925  |  |
| TAD.inner.mean   |  -3.921014 | -2.77422  |  |
| TAD.inner.median |  -0.4689   | -0.88975  |  |
| TAD.inner.q75    |  1.3594    | 0.917325  |  |
| TAD.inner.max    |  3.175     | 2.929     |  |
| TAD.inner.total  |  -27.4471  | -27.7422  |  |

### TAD Boundaries

These are generaed from just pivoting the TAD results so each row is a single TAD boundary, all boundaries are 1 bin long.
Note that only the "start" boundaries are actually called by the tools, the "end" is just the bin before the next called boundary.
If you want to limit the analysis to only "start" boundaries, just filter on the column "boundary.side"
Columns are:
| Column Name    | Example Row 1                         | Example Row 2                         | Column Description | 
|------------    | -------------                         | -------------                         | ------------------ | 
| resolution     | 10000                                 | 10000                                 | resolution |
| normalization  | balanced                              | balanced                              | input matrix normalization |
| TAD.method     | hiTAD                                 | hiTAD                                 | tool used to call TADs |
| TAD.params     | NA                                    | NA                                    | tool-specific params used to call TADs |
| TAD.metric     | ADI                                   | ADI                                   | tool-specific metric used to calculated boundary score |
| Sample.Group   | All.iN.WT                             | All.iN.WT                             | biological condition of input data |
| SampleID       | All.iN.WT.Merged.Merged               | All.iN.WT.Merged.Merged               | specific Sample boundaries were called in |
| isMerged       | Merged                                | Merged                                | whether merged matrices were used as input to call TADs |
| FeatureID      | 10000#hiTAD#NA#ADIchr1#900000#1030000 | 10000#hiTAD#NA#ADIchr1#900000#1030000 | uniqueID for identifying individual TADs that Boundaries belong to|
| boundary.side  | start                                 | end                                   | whether boundary is the start or end of a TAD |
| chr            | chr1                                  | chr1                                  | chromosome |
| start          | 900000                                | 1030000                               | boundary start coord |
| end            | 910000                                | 1040000                               | start + resolution |
| boundary.score | -5.734                                | 0.4358                                | tool-specific score for the boundary bin |
| TAD.length     | 130000                                | 130000                                | size of TAD in bp this boundary is part of |
| TAD.bins       | 13                                    | 13                                    | size of TAD in bins this boundary is part of |

### TAD Measure of Concordance (MoC)

For each pair of conditions we can compare how similar a set of TAD annotations are by calculating the MoC. We limit comparisons to be between the same resolution, chr and TAD calling method + params, so we only compare biological conditions. Note that numerator/denominator are arbitrary, the MoC is not directionally specific.
Columns are:
| Column Name              | Example Row 1      | Example Row 2     | Column Description | 
|------------              | -------------      | -------------     | ------------------ | 
| resolution               | 10000              | 10000             | resolution |
| TAD.method               | hiTAD              | hiTAD             | tool used to call TADs |
| TAD.params               | NA                 | NA                | tool-specific params used to call TADs |
| TAD.metric               | ADI                | ADI               | tool-specific metric used to calculated boundary score |
| isMerged                 | Merged             | Merged            | whether merged matrices were used as input to call TADs |
| Sample.Group.Numerator   | CTCF.iN.BIALLELIC  | CTCF.iN.BIALLELIC | one of the biological conditions who's TADs are being compared |
| n.TADs.Numerator         | 1311               | 1311              | total number of TADs from numerator condition |
| Sample.Group.Denominator | All.iN.WT          | CTCF.iN.DEL       | |
| n.TADs.Denominator       | 1414               | 1276              | |
| chr                      | chr1               | chr1              | chromosome |
| n.Overlaps               | 3401               | 2559              | total number of pairs of overlapping TADs used to calculate MoC |
| MoC                      | 0.8409605634496815 | 0.669946690981735 | calculate MoC value |

### TADCompare 

We also compare TADs called using the r package TADCompare. It calcuatles a bin-wise spectral decomposition-based "boundary score" and comapre these scores around each TAD boundary used as input. The TAD Boundaries come from the called TADs above. The results include 
Gap.Score is a directionally specific z-score, the sign will line up with the "Enriched.Condition" column i.e. -ve => enriched in denominator and vice versa. 
Columns are:
| Column Name              | Example Row 1        | Example Row 2         | Column Description | 
|------------              | -------------        | -------------         | ------------------ | 
| resolution               | 10000                | 10000                 | resolution |
| normalization            | balanced             | balanced              | input matrix normalization |
| TADCompare.params        | 3#15#0.2             | 3#15#0.2              | TADCompare specific params |
| TAD.method               | hiTAD                | hiTAD                 | tool used to call TADs |
| TAD.params               | NA                   | NA                    | tool-specific params used to call TADs |
| Sample.Group.Numerator   | CTCF.iN.BIALLELIC    | CTCF.iN.BIALLELIC     | biological condition of numerator in comparison
| Sample.Group.Denominator | All.iN.WT            | All.iN.WT             | |
| chr                      | chr1                 | chr1                  | chromosome |
| start                    | 1940000              | 2900000               | start of TAD boundary bin|
| end                      | 1950000              | 2910000               | start + resolution |
| isBoundary               | FALSE                | FALSE                 | whether the bin was TAD Boundary called by TAD.method in either condition |
| DifferenceType           | Shifted              | Shifted               | Type of difference detected by TADCompare |
| Enriched.Condition       | All.iN.WT            | All.iN.WT             | which condition is enriched for the TAD boundary difference |
| TAD.Score.Numerator      | 0.8797242872092477   | 0.6460212128404218    | bin-specific boundary score for the boundary being compared |
| TAD.Score.Denominator    | 1.9732218728763105   | 1.817865422191164     | |
| Gap.Score                | -3.4442393791143675  | -3.726171792796507    | z-score of the difference between TAD.Score.* |
| p.value                  | 5.726684378441336e-4 | 1.9440992804322307e-4 | p-value calculated from Gap.Score |
| p.adj.gw                 | 0.06789254167713184  | 0.03075384322974094   | BH-adjusted p-value for all compared boundaries across the genome per condition |

### Generate TAD Results

Commands to generate TAD + TAD Comparison data
```bash
# generate commands to run shell tools for TAD calling
mamba activtte r
Rscript ./scripts/TADs/make.TAD.calling.cmds.R -t $(nproc)
# now run those generated commands with GNU parallel
mamba activtte TADs
parallel -j 1 --eta --bar :::: ./results/TADs/all.TAD.calling.cmds.txt
# Generate Consensus TAD results from set of individual matrices with spectralTAD
mamba activtte r
Rscript ./scripts/TADs/run.ConsensusTADs.R
# Coallate TAD results into single, structured output files
mamba activtte r
Rscript ./scripts/TADs/coallate.all.TAD.results.R
# Compute TAD MoCs for all sets of TADs
mamba activtte r
Rscript ./scripst/TADs/calculate.TAD.MoCs.R
# Run TADCompare to generated differential TAD results
# requires 120Gb for the largest matrix comparison (i.e. chr1 @5Kb)
mamba activtte r
Rscript ./scripts/TADs/run.TADCompare.R -t $(nproc)
```

Ultimately these commands will generate the followig output files
```bash
./results/TADs/
├── results_TADs/                # Nested directory structure of individually generated results
├── all.ConsensusTAD.TADs.tsv    # combined file with all TADs called by ConsensusTAD
├── all.hiTAD.TADs.tsv           # combined file with all TADs called by hiTAD
├── all.cooltools.TADs.tsv       # combined file with all TADs called by cooltools
├── all.all.TADs.tsv             # combined file with all TADs called by all methods
├── all.all.TAD.MoCs.tsv         # Computed Measure of Concordance of TADs for all pairs of conditions
├── results_TADCompare/
├── all.TADCompare.results.tsv
└── all.TADCompare.n.results.tsv
```

