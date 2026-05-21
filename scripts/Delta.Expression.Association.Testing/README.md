## HiF ~ Differential Gene Association Testing 

### Generate Association data

### Generate Testing Results

### Generate Association Testing Results

```bash
# Parse + clean up + combine gene-linked functional loci annotations that define "indirect assciations"
Rscript ./scripts/Delta.Expression.Association.Testing/coallate.all.gene.associated.functioal.loci.R
# Generate clean association tables linking HiFs called in individual conditions to DEG results directly and indirectly
Rscript ./scripts/Delta.Expression.Association.Testing/link.HiFs.to.Genes.R
# Generate clean association tables linking differential HiFs called between conditions to DEG results directly and indirectly
Rscript ./scripts/Delta.Expression.Association.Testing/link.diff.HiFs.to.Genes.R
# Loop over generated association data to calculate statistical test results + stratifications
Rscript ./scripts/Delta.Expression.Association.Testing/test.HiF.DEG.functional.associations.R
```
