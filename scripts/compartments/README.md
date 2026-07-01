## Compartment Analysis

### Compartment Annotation

We call compartments using `cooltools eigs-cis`, and phase the calls using `genecov` track, compuated by `cooltools genome genecov`, using all default params @ 100,50,25,10,5Kb.
```bash
./compartments/
├── all.compartment.calling.cmds.txt
├── all.compartment.bins.tsv 
├── all.compartments.tsv
└── all.compartment.switches.tsv
```

The binwise PC1 scores are the most granular data, can compute other summaries from this. Columns are:
| Column Name             | Example Row 1                 | Example Row 2                 | Column Description |
|------------             | -------------                 | -------------                 | ------------------ |
| compartment.params      | cooltools#genecov#balanced#E1 | cooltools#genecov#balanced#E1 | hyper-params for how comparments were called |
| Comp.method             | cooltools                     | cooltools                     | tools used |
| track.type              | genecov                       | genecov                       | track type used to orient PC1 scores |
| normalization           | balanced                      | balanced                      | wheter input matric was balanced or not |
| resolution              | 100000                        | 100000                        |                    |
| isMerged                | Merged                        | Merged                        | is input matrix merged |
| SampleID                | All.iN.WT.Merged.Merged       | All.iN.WT.Merged.Merged       |  SampleId of matrix compartments were called on |
| Sample.Group            | All.iN.WT                     | All.iN.WT                     |                    |
| chr                     | chr1                          | chr1                          |                    |
| start                   | 600000                        | 900000                        |                    |
| end                     | 700000                        | 1000000                       |                    |
| score.source            | E1                            | E1                            |                    |
| score                   | -0.06375528953866347          | 1.4529687030279492            |                    |
| score.change            | 1.5167239925666127            | -0.20742163683109704          |                    |
| n.compartment.lvls      | 20                            | 20                            | how many levels abs(PC1 score) was binned into for strenght.lvl |
| strength.lvl            | 1                             | 12                            | qunatile of abs(PC1) strenght of the bin |
| strength.lvl.change     | 1                             | -2                            | stenghth.lvl difference between this bin and the next adjacent bin (not next row, adjacent bin) |

| compartment             | B                             | A                             | annotated compartment for the current bin |
| compartment.change      | B->A                          | A->A                          | compartment of next adjacent bin |
| does.compartment.switch | TRUE                          | FALSE                         | is the next adjacent bin assigned a different compartment? |
| compartment.switch.type | switch                        | weaker                        | is the strength of the next bin either i) a switch or or ii) same compartment, but with a strenght change |


Also the switches and the difference in PC1 betweeen each switch and its next adjacent bin (which must be in a differential compartment)

### Compartment Comparisons

We can compare compartments in multiple ways. The results files are:
```bash
./compartments/
├── all.compartment.MoCs.tsv
├── all.compartment.segment.test.results.tsv
├── all.compartment.switch.comparisons.tsv
└── all.saddle.data.cmds.txt
```

#### Measure of Concordance 

First compare whether comaprtments overlap by computing MoC (same as with TADs) per compartment between conditions
Columns are:
| Column Name        | Example Row 1                 | Example Row 2                 | Column Description |
|------------        | -------------                 | -------------                 | ------------------ |
| Sample.Group.P1    | NIPBL.iN.DEL                  | NIPBL.iN.DEL                  | Conditions being compared
| Sample.Group.P2    | NIPBL.iN.WT                   | NIPBL.iN.WT                   |                    |
| resolution         | 50000                         | 100000                        |                    |
| compartment.params | cooltools#genecov#balanced#E1 | cooltools#genecov#balanced#E1 |                    |
| compartment        | B                             | B                             | which type of compartments are being compared |
| chr                | chr1                          | chr1                          |                    |
| n.Overlaps         | 136                           | 80                            | number of overlapping pairs of compartments |
| MoC                | 0.3061653019480192            | 0.3844047347551683            | computed measure of conrcordance |
| n.regions.P1       | 138                           | 83                            | number of compartments in condition P1 |
| n.regions.P2       | 164                           | 90                            | number of compartments detected in condition P2 |


### Compare PC1 scores binwise

1. For each compartment in each denominator condition, compare PC1 scores for the corresponding bins in the numerator condition 
   - compute summary stats over the binwise differenceds in PC1
   - compute differences test results (e.g. KS-test) between PC1 scores per compartment
segment-wise test results:
| Column Name              | Example Row 1                 | Example Row 2                 | Column Description |
|------------              | -------------                 | -------------                 | ------------------ |
| Sample.Group.Numerator   | NIPBL.iN.DEL                  | NIPBL.iN.DEL                  | Conditions being compared |
| Sample.Group.Denominator | NIPBL.iN.WT                   | NIPBL.iN.WT                   |                    |
| resolution               | 100000                        | 100000                        |                    |
| compartment.params       | cooltools#genecov#balanced#E1 | cooltools#genecov#balanced#E1 |                    |
| chr                      | chr1                          | chr1                          |                    |
| start                    | 900000                        | 900000                        | start of compartment in WT |
| end                      | 1400000                       | 1400000                       | end of compartment in WT |
| SegmentID                | chr1#900000#1400000           | chr1#900000#1400000           |                    |
| compartment.Denominator  | A                             | A                             | type of compartment in WT |
| n.bins                   | 5                             | 5                             | lenght of WT compartment in bins | 
| feature                  | n.bins.diff.compartment       | cosine.dist                   | type of statistic computed to compare binwise PC1 scores between conditions |
| value                    | 0                             | 0.9988142119718464            | computed value of the listed statistic |

#### Switch Differences

We can also compare the strength an location of switches bewteen conditions. Here I map every switch in WT to its nearest switch in DELs, so we can plot the distance and "strength difference" difference between conditions.
Columns are:
| Column Name                    | Example Row 1                 | Example Row 2                 | Column Description |
|------------                    | -------------                 | -------------                 | ------------------ |
| Sample.Group.Numerator         | NIPBL.iN.DEL                  | NIPBL.iN.DEL                  | Conditions being compared |
| Sample.Group.Denominator       | NIPBL.iN.WT                   | NIPBL.iN.WT                   |
| resolution                     | 100000                        | 100000                        |
| compartment.params             | cooltools#genecov#balanced#E1 | cooltools#genecov#balanced#E1 |
| compartment.change.Numerator   | B->A                          | B->A                          | compartment change of numerator switch |
| compartment.change.Denominator | A->B                          | A->B                          | Compartmetn change of denominator switch
| chr                            | chr1                          | chr1                          | |
| start                          | 17500000                      | 17500000                      | |
| end                            | 17600000                      | 17600000                      | |
| distance                       | 3000000                       | 3000000                       | distance between the two swtiches |
| feature                        | score                         | score.change                  | swithc feature being compared between switches |
| diff                           | -0.5302659900969029           | 0.6823725849841974            | difference in feature values of switches (x-y) |
| FC                             | -0.13015713529876818          | -0.28082016544743144          | fold-change of feautre values of swiches (x/y) |

#### Saddle Data

Compute saddle plot data to ensure that A and B comparments preferentiall contacts other A/B compartments respectively (as expected 

### Generate Compartment Results

Generate Compartment annotations
```bash
Rscritp scripts/compartments/make.compartment.calling.cmds.R
parallel -j 1 --eta --bar ::::: ./results/compartments/all.compartment.calling.cmds.txt
# combine individual results files into a single tidy file
Rscript scripts/compartments/coallate.all.compartment.results.R
```
Compute comparison results of PC1 scores between conditions
```bash
# compute differential compartment statistics between conditions
Rscript scripts/compartments/compute.compartment.comparison.results.R
```
Generate saddle statistics data
```bash
Rscritp scripts/compartments/make.saddle.data.cmds.R
parallel -j 1 --eta --bar ::::: ./results/compartments/all.saddle.data.cmds.txt
```

