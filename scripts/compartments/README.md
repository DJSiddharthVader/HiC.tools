## Compartment Analysis

### Compartment Annotation

We call compartments using `cooltools eigs-cis`, and phase the calls using `genecov` track, compuated by `cooltools genome genecov`, using all default params @ 100,50,25,10,5Kb.

### Generate Compartment Results

Generate Compartment annotations
```bash
Rscritp scripts/compartments/make.compartment.calling.cmds.R
parallel -j 1 --eta --bar ::::: ./results/compartments/all.compartment.calling.cmds.txt
```
Generate saddle statistics data
```bash
Rscritp scripts/compartments/make.saddle.data.cmds.R
parallel -j 1 --eta --bar ::::: ./results/compartments/all.saddle.data.cmds.txt
```

