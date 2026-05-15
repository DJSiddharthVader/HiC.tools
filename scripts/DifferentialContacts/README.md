## Differential Interacting Regiosn (DIRs) 

### Differential Contact Analysis

We use the [multiHiCCompare](https://bioconductor.org/packages/release/bioc/vignettes/multiHiCcompare/inst/doc/multiHiCcompare.html) R package to look for individual difference in piexl contacts between conditions. (i.e. number of contacts between a specific pair of bins).
Before testing pixels, we first filter out any pixels that are not detected in at least 80% of input matrices (i.e. 0 paired-end reads supporting that specific pixel). So if we are comparing 6 WTs vs 6 DELs then we only test pixels that are detected in at least 10 samples total.
We also preform two sets of comparisons 1) comapring all replicates between a pair of conditions and 2) comapring a pair of merged matrices, 1 merged matrix per condition. Since this is an explicitly differential analysis, directly modeling the variance within a condition (i.e. beteween technical + biological replicates) is valuable and leads to more stable results. 
Using a GLM it tests whether the contacts for a specifi bin-pair (i.e. contact matrix pixel) are differential between 2 conditions e.g. DEL vs WT, acounting for distance-expectation and sparsity/low-coverge pixels. 
We use the `fastlo()` implementation of the cyclic LOESS normalizaiton and we use the function `hic_exactTest()` since we are only every comparing a single binary condition (e.g. DEL vs WT).
We calulate results separately for every resolution (@ 100,50,25,10Kb) + comparison + chromosome, and pool results genome-wide (i.e. by resolution by comparison) to apply BH adjustment to the raw p-values.

### Generating results

Run the script, it will generate results files, 1 per combination of comparison + chr + resolution + merged vs individual + pixel filtering thresholds.
```bash
Rscript ./scripts/DifferentialContacts/run.multiHiCCompare.R
```
Now we can combine all these results into a single tidy file `./results/DifferentialContacts/all.multiHiCCompare.results.tsv`
```bash
Rscript ./scripts/DifferentialContacts/coallate.multiHiCCompare.results.R
```
The combined results have the following columns:

| Column Name            | Example Row 1               | Example Row 2             | Column Description | 
| -----------            | -------------               | -------------             | ------------------ | 
| zero.p                 | 0.8                         | 0.8                       | filtering by the proportion of zero pixels across samples for an interaction | 
| A.min                  | 5                           | 5                         | filtering by a minimum average pixel value across samples |
| merged                 | Individual                  | Individual                | whether input is 2 merged matrices or 2 sets of individual replicate matrices |
| resolution             | 25000                       | 25000                     | resolution DIRs are called at |
| SampleID.Numerator     | CTCF.iN.BIALLELIC           | CTCF.iN.BIALLELIC         | Comparison info, +ve log2(FC) means increased contacts in numerator, relative to denominator |
| SampleID.Denominator   | All.iN.WT                   | All.iN.WT                 | ~ |
| Edit.Numerator         | CTCF                        | CTCF                      | ~ |
| Celltype.Numerator     | iN                          | iN                        | ~ |
| Genotype.Numerator     | BIALLELIC                   | BIALLELIC                 | ~ |
| Edit.Denominator       | All                         | All                       | ~ |
| Celltype.Denominator   | iN                          | iN                        | ~ |
| Genotype.Denominator   | WT                          | WT                        | ~ |
| bin.pair.idx           | chr10#104325000#104425000   | chr10#131975000#132325000 | unique index for each pixel in the genome (chr#anchor.left#anchor.right) |
| chr                    | chr10                       | chr10                     | chromosome |
| region1                | 104325000                   | 131975000                 | left anchor position in bp |
| region2                | 104425000                   | 132325000                 | right anchor position in bp |
| distance.bins          | 4                           | 14                        |  distance between anchors in number of bins |
| logFC                  | 1.381440110604435           | 1.5453056614655554        | log2(fold-change) between pixel contacts  in numerator vs denominator | 
| logCPM                 | 5.354427957484957           | 6.005682413414727         | log(contacts-per-million) for the pixel across all samples from both conditions |
| p.value                | 2.2418922644089334e-20      | 3.731730920115992e-19     | differential pixel raw p-value |
| p.adj                  | 3.376065560973413e-16       | 7.825439739483235e-15     | BH adjusted p-value, adjust pixels from each chr separately |                     
| p.adj.gw               | 4.34455180973675e-15        | 4.167963124443381e-14     | BH adjusted p-value, adjust 
| log.p.value            | 19.649385261543422          | 18.42808967899892         | log10() of each specified p-value above |                     
| log.p.adj              | 15.471589128251205          | 14.106491248533576        | ~ |                     
| log.p.adj.gw           | 14.362055019352907          | 13.380076132060164        | ~ |                     

