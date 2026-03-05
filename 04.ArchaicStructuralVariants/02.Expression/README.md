**Detailed descriptions are provided in the manuscript.**

# Requirements
Before running the pipeline, you should make sure these software/packages are installed:

- edgeR
- 


# RNA-seq reads mapping
@anguo

# Gene expression quantifing and normalization
@anguo

After getting rawly quantified gene expression levels, we subsequently normalized them into TMM via `edgeR`
```
library("edgeR")

stats <- read.table("input.gene_counts",header = T,sep = "\t", check.names = F)
y <- DGEList(counts = stats[, 2:dim(stats)[2]], genes = stats[, 1])
y <- calcNormFactors(y)
cpms <- cpm(y,log = F)
cpms_with_genes <- cbind(Gene = stats[, 1], cpms)
write.table(cpms_with_genes, file = "input.gene_counts.TMM", sep = "\t", row.names = FALSE, quote = FALSE)
```
