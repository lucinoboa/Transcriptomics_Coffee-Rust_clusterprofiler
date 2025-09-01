# Transcriptomics_Coffee-Rust_clusterprofiler
Description

## 1. Load the required libraries.
```r
library(DESeq2)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(org.At.tair.db)   # Replace with Coffea annotation DB if available
library(tidyverse)
```

## 2. Set working directory to where the featureCounts matrix is stored
```r
setwd("D:/lucianoboa/royatranscriptomics/analysis/featureCounts")
```

## 3. Import raw counts matrix
```r
countData <- read.table("counts_matrix_complete_royatranscriptomics.txt",
                        header = TRUE, row.names = 1, sep = "\t")

# Keep only count columns (adjust indices to your file)
countData <- countData[, 6:21]

# Rename columns to simpler sample IDs
colnames(countData) <- c("H10","H11","H12","H13","H14","H15","H16","H9",
                         "T1","T2","T3","T4","T5","T6","T7","T8")

## 4. Define experimental groups
group <- rep(NA, ncol(countData))
names(group) <- colnames(countData)

group[c("T1","T2","T3","T4","T5","T6","T7","T8")] <- "Group_1"
group[c("H9","H11","H13","H14","H15","H16")]      <- "Group_2"
group[c("H10","H12")]                             <- "Group_3"

group <- factor(group)
```

## 5. Create colData with sample metadata
```r
colData <- data.frame(
  row.names = colnames(countData),
  condition = group
)
```

## 6. Build DESeq2 object
```r
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData   = colData,
  design    = ~ condition
)
```

## 7. Filter out very lowly expressed genes
```r
dds <- dds[rowSums(counts(dds)) > 10, ]
```

## 8. Run DESeq normalization and modeling
```r
dds <- DESeq(dds)
```

## 9. Extract contrasts between conditions
```r
res_G2_vs_G1 <- results(dds, contrast = c("condition", "Group_2", "Group_1"))
res_G3_vs_G1 <- results(dds, contrast = c("condition", "Group_3", "Group_1"))
res_G2_vs_G3 <- results(dds, contrast = c("condition", "Group_2", "Group_3"))
```

## 10. Apply filtering criteria (padj < 0.05 & |log2FC| > 1)

### Group_2 vs Group_1
```r
deg_G2_vs_G1 <- subset(res_G2_vs_G1, padj < 0.05 & abs(log2FoldChange) > 1)
up_G2_vs_G1 <- subset(deg_G2_vs_G1, log2FoldChange > 1)
print(up_G2_vs_G1)
summary(up_G2_vs_G1)
down_G2_vs_G1 <- subset(deg_G2_vs_G1, log2FoldChange < -1)
print(down_G2_vs_G1)
summary(down_G2_vs_G1)
```

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

### Group_3 vs Group_1
```r
deg_G3_vs_G1 <- subset(res_G3_vs_G1, padj < 0.05 & abs(log2FoldChange) > 1)
up_G3_vs_G1 <- subset(deg_G3_vs_G1, log2FoldChange > 1)
print(up_G3_vs_G1)
summary(up_G3_vs_G1)
down_G3_vs_G1 <- subset(deg_G3_vs_G1, log2FoldChange < -1)
print(down_G3_vs_G1)
summary(down_G3_vs_G1)
```

```r
# print(up_G3_vs_G1)
log2 fold change (MLE): condition Group_3 vs Group_1 
Wald test p-value: condition Group 3 vs Group 1 
DataFrame with 1897 rows and 6 columns
              baseMean log2FoldChange     lfcSE      stat      pvalue        padj
             <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
LOC113741252   412.415        1.76366  0.543072   3.24756 1.16399e-03  0.01464447
LOC113742115   407.926        1.97722  0.694882   2.84540 4.43557e-03  0.03953960
LOC113732466  3101.236        1.56046  0.514725   3.03163 2.43234e-03  0.02538484
LOC113742438  3631.438        4.40771  1.115847   3.95010 7.81184e-05  0.00180478
LOC113691221  1254.375        1.19275  0.335174   3.55859 3.72852e-04  0.00603730
...                ...            ...       ...       ...         ...         ...
LOC113722912   892.204        1.60588  0.481418   3.33573 8.50753e-04 0.011426050
LOC113719025   387.244        1.74022  0.441461   3.94196 8.08172e-05 0.001840978
LOC113717824   235.603        3.76434  1.219002   3.08805 2.01476e-03 0.022119169
LOC113719365   327.516        1.57965  0.507100   3.11507 1.83903e-03 0.020701190
LOC113719217   462.872        1.00932  0.244271   4.13195 3.59701e-05 0.000966077

# summary(up_G3_vs_G1)

out of 1897 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1897, 100%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 9)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

# print(down_G3_vs_G1)
log2 fold change (MLE): condition Group_3 vs Group_1 
Wald test p-value: condition Group 3 vs Group 1 
DataFrame with 1986 rows and 6 columns
              baseMean log2FoldChange     lfcSE      stat      pvalue        padj
             <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
LOC113741337 4943.9253       -1.60467  0.379097  -4.23287 2.30730e-05 0.000671612
LOC113741513 2248.8095       -1.48210  0.396469  -3.73826 1.85300e-04 0.003512902
LOC113741638   62.8550       -3.73787  0.925385  -4.03926 5.36202e-05 0.001337966
LOC113734549 1342.5441       -1.21913  0.336343  -3.62466 2.89339e-04 0.004935160
LOC113690361   18.1007       -7.43403  2.258591  -3.29145 9.96733e-04 0.012864664
...                ...            ...       ...       ...         ...         ...
LOC113717553  420.7201       -1.94575  0.609438  -3.19270 1.40948e-03 1.69519e-02
LOC113719438 2139.8028       -1.11992  0.248078  -4.51439 6.34988e-06 2.34232e-04
LOC113691294  430.5862       -1.37843  0.438586  -3.14289 1.67291e-03 1.93112e-02
LOC113718036  578.7005       -1.95680  0.293418  -6.66898 2.57586e-11 6.55426e-09
LOC140034518   62.5567       -3.49675  1.194840  -2.92654 3.42751e-03 3.27131e-02

# summary(down_G3_vs_G1)

out of 1986 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 1986, 100%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 9)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

### Group_2 vs Group_3
```r
deg_G2_vs_G3 <- subset(res_G2_vs_G3, padj < 0.05 & abs(log2FoldChange) > 1)
up_G2_vs_G3 <- subset(deg_G2_vs_G3, log2FoldChange > 1)
print(up_G2_vs_G3)
summary(up_G2_vs_G3)
down_G2_vs_G3 <- subset(deg_G2_vs_G3, log2FoldChange < -1)
print(down_G2_vs_G3)
summary(down_G2_vs_G3)
```

```r
# print(up_G2_vs_G3)
log2 fold change (MLE): condition Group_2 vs Group_3 
Wald test p-value: condition Group_2 vs Group_3 
DataFrame with 1885 rows and 6 columns
              baseMean log2FoldChange     lfcSE      stat      pvalue       padj
             <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
LOC113741638   62.8550        3.28912  0.950235   3.46138 5.37423e-04 0.02135215
LOC113737686  892.1137        1.42073  0.484868   2.93015 3.38802e-03 0.04608584
LOC140006165   48.2494        4.24343  1.333040   3.18327 1.45621e-03 0.03094668
LOC113740114   63.4060        3.53997  0.836220   4.23330 2.30289e-05 0.00351784
LOC113715221  551.7163        1.20288  0.410861   2.92770 3.41477e-03 0.04628937
...                ...            ...       ...       ...         ...        ...
LOC140034370   58.2204        3.31742   1.14561   2.89577  0.00378234  0.0493319
LOC140034371  115.7935        3.25234   1.12360   2.89458  0.00379661  0.0494061
LOC140034372   53.3361        3.66315   1.13080   3.23942  0.00119773  0.0285191
LOC140034373   59.1059        3.72689   1.14405   3.25763  0.00112347  0.0278503
LOC140034374   62.6359        3.58789   1.15058   3.11832  0.00181888  0.0342338

# summary(up_G2_vs_G3)

out of 1885 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 1885, 100%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

# print(down_G2_vs_G3)
log2 fold change (MLE): condition Group_2 vs Group_3 
Wald test p-value: condition Group_2 vs Group_3 
DataFrame with 939 rows and 6 columns
               baseMean log2FoldChange     lfcSE      stat      pvalue      padj
              <numeric>      <numeric> <numeric> <numeric>   <numeric> <numeric>
LOC113740232  282.97783       -1.45196  0.450647  -3.22195  0.00127320 0.0291239
LOC113723694  537.12123       -1.22434  0.396816  -3.08540  0.00203279 0.0362575
LOC113695710    4.84374       -5.87803  1.862323  -3.15629  0.00159791 0.0322861
LOC113695794   11.39724       -7.27834  2.489723  -2.92335  0.00346283 0.0467836
LOC113703325 1154.66147       -1.69421  0.449510  -3.76902  0.00016389 0.0122169
...                 ...            ...       ...       ...         ...       ...
LOC113719613    358.772       -1.57262  0.506506  -3.10483 0.001903882 0.0350913
LOC113716647    582.794       -1.11540  0.311008  -3.58641 0.000335263 0.0169348
LOC113735054   1312.026       -1.15765  0.350956  -3.29856 0.000971811 0.0263914
LOC113694112   1559.010       -2.95066  0.816800  -3.61247 0.000303297 0.0161577
LOC113718708    442.632       -2.41845  0.718281  -3.36700 0.000759897 0.0246344

# summary(down_G2_vs_G3)

out of 939 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 939, 100%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```


## 11. Export DEG tables following strict filtering protocol
```r
### Group_2 vs Group_1
write.csv(as.data.frame(deg_G2_vs_G1), file = "DEG_G2_vs_G1_strict.csv")
write.csv(as.data.frame(up_G2_vs_G1), file = "Up_G2_vs_G1_strict.csv")
write.csv(as.data.frame(down_G2_vs_G1), file = "Down_G2_vs_G1_strict.csv")

### Group_3 vs Group_1
write.csv(as.data.frame(deg_G3_vs_G1), file = "DEG_G3_vs_G1_strict.csv")
write.csv(as.data.frame(up_G3_vs_G1), file = "Up_G3_vs_G1_strict.csv")
write.csv(as.data.frame(down_G3_vs_G1), file = "Down_G3_vs_G1_strict.csv")

### Group_2 vs Group_3
write.csv(as.data.frame(deg_G2_vs_G3), file = "DEG_G2_vs_G3_strict.csv")
write.csv(as.data.frame(up_G2_vs_G3), file = "Up_G2_vs_G3_strict.csv")
write.csv(as.data.frame(down_G2_vs_G3), file = "Down_G2_vs_G3_strict.csv")
```
