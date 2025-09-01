# Transcriptomics_Coffee-Rust_clusterprofiler
Description

# 1. Load the required libraries.
```r
library(DESeq2)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(org.At.tair.db)   # Replace with Coffea annotation DB if available
library(tidyverse)
```

# 2. Set working directory to where the featureCounts matrix is stored
```r
setwd("D:/lucianoboa/royatranscriptomics/analysis/featureCounts")
```

# 3. Import raw counts matrix
```r
countData <- read.table("counts_matrix_complete_royatranscriptomics.txt",
                        header = TRUE, row.names = 1, sep = "\t")

# Keep only count columns (adjust indices to your file)
countData <- countData[, 6:21]

# Rename columns to simpler sample IDs
colnames(countData) <- c("H10","H11","H12","H13","H14","H15","H16","H9",
                         "T1","T2","T3","T4","T5","T6","T7","T8")

# 4. Define experimental groups
group <- rep(NA, ncol(countData))
names(group) <- colnames(countData)

group[c("T1","T2","T3","T4","T5","T6","T7","T8")] <- "Group_1"
group[c("H9","H11","H13","H14","H15","H16")]      <- "Group_2"
group[c("H10","H12")]                             <- "Group_3"

group <- factor(group)
```

# 5. Create colData with sample metadata
```r
colData <- data.frame(
  row.names = colnames(countData),
  condition = group
)
```

# 6. Build DESeq2 object
```r
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)
```

# 7. Filter out very lowly expressed genes
```r
dds <- dds[rowSums(counts(dds)) > 10, ]
```

# 8. Run DESeq normalization and modeling
```r
dds <- DESeq(dds)
```

# 9. Extract contrasts between conditions
```r
res_G2_vs_G1 <- results(dds, contrast = c("condition", "Group_2", "Group_1"))
res_G3_vs_G1 <- results(dds, contrast = c("condition", "Group_3", "Group_1"))
res_G2_vs_G3 <- results(dds, contrast = c("condition", "Group_2", "Group_3"))
```

# 10. Apply filtering criteria (padj < 0.05 & |log2FC| > 1)
```r
deg_G2_vs_G1 <- subset(res_G2_vs_G1, padj < 0.05 & abs(log2FoldChange) > 1)
up_G2_vs_G1 <- subset(deg_G2_vs_G1, log2FoldChange > 1)
print(up_G2_vs_G1)
summary(up_G2_vs_G1)
down_G2_vs_G1 <- subset(deg_G2_vs_G1, log2FoldChange < -1)
print(down_G2_vs_G1)
summary(down_G2_vs_G1)

```r
# print(up_G2_vs_G1)
log2 fold change (MLE): condition Group_2 vs Group_1 
Wald test p-value: condition Group 2 vs Group 1 
DataFrame with 3422 rows and 6 columns
               baseMean log2FoldChange     lfcSE      stat      pvalue        padj
              <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
LOC113741252   412.4153        1.30811  0.371678   3.51947 4.32417e-04 2.53074e-03
LOC113741420    98.5419        1.37903  0.506075   2.72494 6.43128e-03 2.52903e-02
LOC113742192   116.4594        6.84668  0.942959   7.26085 3.84659e-13 1.52541e-11
LOC113742438  3631.4381        4.85972  0.762585   6.37270 1.85729e-10 4.71183e-09
LOC113688219    40.7452        1.36299  0.509389   2.67573 7.45657e-03 2.84774e-02
...                 ...            ...       ...       ...         ...         ...
CoarCp002       8.77996        2.55909  0.952475   2.68678  0.00721437   0.0277172
CoarCp011     197.52014        1.18503  0.435629   2.72029  0.00652257   0.0255670
CoarCp013      80.82201        1.52478  0.559059   2.72740  0.00638347   0.0251354
CoarCp020    3325.26408        1.24010  0.510120   2.43100  0.01505722   0.0486725
CoarCp022     106.04523        1.25783  0.480158   2.61961  0.00880300   0.0324645

# summary(up_G2_vs_G1)
out of 3422 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 3422, 100%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

# print(down_G2_vs_G1)
log2 fold change (MLE): condition Group_2 vs Group_1 
Wald test p-value: condition Group 2 vs Group 1 
DataFrame with 4666 rows and 6 columns
              baseMean log2FoldChange     lfcSE      stat      pvalue        padj
             <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
LOC113740232   282.978       -1.83803  0.299178  -6.14359 8.06746e-10 1.81509e-08
LOC113732076   683.863       -1.22732  0.308030  -3.98441 6.76492e-05 5.02582e-04
LOC113741337  4943.925       -2.32863  0.259013  -8.99039 2.46351e-19 2.41482e-17
LOC113741513  2248.809       -2.16801  0.270904  -8.00286 1.21567e-15 7.00191e-14
LOC113743322   303.338       -1.78522  0.295549  -6.04036 1.53770e-09 3.28783e-08
...                ...            ...       ...       ...         ...         ...
LOC113718036 578.70047       -1.72956  0.198105  -8.73048 2.53587e-18 2.10522e-16
LOC140032825   5.82434       -5.87491  1.652946  -3.55421 3.79122e-04 2.25582e-03
LOC113719623 356.76611       -1.62249  0.275877  -5.88120 4.07305e-09 8.04933e-08
LOC113717559 500.93778       -1.16120  0.301700  -3.84885 1.18675e-04 8.25936e-04
LOC113719218   7.64243       -5.61186  2.059496  -2.72487 6.43265e-03 2.52933e-02

# summary(down_G2_vs_G1)

out of 4666 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 4666, 100%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 1)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

```



```r
deg_G3_vs_G1 <- subset(res_G3_vs_G1, padj < 0.05 & abs(log2FoldChange) > 1)
up_G3_vs_G1 <- subset(deg_G3_vs_G1, log2FoldChange > 1)
down_G3_vs_G1 <- subset(deg_G3_vs_G1, log2FoldChange < -1)
```

```

deg_G2_vs_G3 <- subset(res_G2_vs_G3, padj < 0.05 & abs(log2FoldChange) > 1)
up_G2_vs_G3 <- subset(deg_G2_vs_G3, log2FoldChange > 1)
down_G2_vs_G3 <- subset(deg_G2_vs_G3, log2FoldChange < -1)


# 11. Export DEG tables following strict filtering protocol

## Group_2 vs Group_1

write.csv(as.data.frame(deg_G2_vs_G1), file = "DEG_G2_vs_G1_strict.csv")
write.csv(as.data.frame(up_G2_vs_G1), file = "Up_G2_vs_G1_strict.csv")
write.csv(as.data.frame(down_G2_vs_G1), file = "Down_G2_vs_G1_strict.csv")

## Group_3 vs Group_1
write.csv(as.data.frame(deg_G3_vs_G1), file = "DEG_G3_vs_G1_strict.csv")
write.csv(as.data.frame(up_G3_vs_G1), file = "Up_G3_vs_G1_strict.csv")
write.csv(as.data.frame(down_G3_vs_G1), file = "Down_G3_vs_G1_strict.csv")

## Group_2 vs Group_3
write.csv(as.data.frame(deg_G2_vs_G3), file = "DEG_G2_vs_G3_strict.csv")
write.csv(as.data.frame(up_G2_vs_G3), file = "Up_G2_vs_G3_strict.csv")
write.csv(as.data.frame(down_G2_vs_G3), file = "Down_G2_vs_G3_strict.csv")
