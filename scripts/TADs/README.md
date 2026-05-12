# Generating TAD Results

## Commands to Run

Commands to generate TAD + TAD Comparison data
```bash
# generate commands to run shell tools for TAD calling
Rscript ./scripts/TADs/make.TADCaller.cmds.R -t 8 # use 8 threads
# now run those generated commands with GNU parallel
parallel -j 1 --eta --bar ::::: ./results/TADs/all.TAD.calling.cmds.txt
# Generate Consensus TAD results from set of individual matrices with spectralTAD
Rscript ./scripts/TADs/run.ConsensusTADs.R
# Coallate TAD results into single, structured output files
Rscript ./scripts/TADs/coallate.TAD.results.R
# Compute TAD MoCs for all sets of TADs
Rscript ./scripst/TADs/calculate.TAD.MoCs.R
# Run TADCompare to generated differential TAD results
# requires 120Gb for the largest matrix comparison (i.e. chr1 @5Kb)
Rscript ./scripts/TADs/run.TADCompare.R -t 8 
```

## Output File Descriptions

The TAD results files are as follows:
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

### TAD Annotation

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

### TAD MoC Data

### TADCompare Data

For 2 different TAD annotation sets (start/end pairs) of the same region (e.g. chr16) we can calculate the similarity of the 2 sets by computing the Measure of Concordance as defined in [this paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1596-9#Sec21). 

We also compare TADs called from merged matrices using `TADCompare`, since each merged matrix represents a biological condition (e.g. 16p.iN.WT) it is an expliciit differential analysis 

