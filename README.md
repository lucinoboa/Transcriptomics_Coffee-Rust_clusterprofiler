# Transcriptomics_Coffee-Rust_enrichment

# Differential Gene Expression Analysis of Royal Coffee Transcriptome: Low vs High Severity
This project presents a comprehensive workflow for analyzing transcriptomic data from coffee plants infected with rust (Hemileia vastatrix) under different disease severity conditions (Low vs High Severity). The analysis includes:
- Loading raw count matrices and defining experimental groups.
- Differential expression analysis using DESeq2.
- Filtering significant differentially expressed genes (DEGs) with strict criteria (adjusted p-value < 0.05, |log2 fold change| > 1).
- Annotating DEGs with functional categories (COG) using a custom annotation file.
- Generating a summarized functional table comparing up- and down-regulated genes across categories.

## 0. Load required libraries.
```r
library(DESeq2)       # For differential expression analysis
library(tidyverse)    # For data manipulation, joining, summarizing
library(readr)        # Reading CSV/TSV files
library(dplyr)        # Filtering, summarizing
```

## 1. Set working directory
```r
setwd("D:/lucianoboa/royatranscriptomics/analysis/featureCounts")
```

## 2. Load count matrix
```r
countData <- read.table("counts_matrix_complete_royatranscriptomics.txt",
                        header = TRUE, row.names = 1, sep = "\t")
```

### Select only count columns (columns 6-21)
```r
countData <- countData[, 6:21]
```

### Rename columns for simplicity
```r
colnames(countData) <- c("H10","H11","H12","H13","H14","H15","H16","H9",
                         "T1","T2","T3","T4","T5","T6","T7","T8")
```

## 3. Define experimental groups (Low vs High severity + Control)
```r
group <- rep(NA, ncol(countData))
names(group) <- colnames(countData)

group[c("H9","H10","H11","H12","H16")] <- "Low_Severity"
group[c("H13","H14","H15")] <- "High_Severity"
group[c("T1","T2","T3","T4","T5","T6","T7","T8")] <- "Control"

# Convert to factor
group <- factor(group)
colData <- data.frame(row.names = colnames(countData), group = group)
```

## 4. DESeq2 analysis
```r
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ group)

# Filter out low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq2
dds <- DESeq(dds)
```

# 5. Extract results for Low vs High Severity
```r
res_Low_vs_High <- results(dds, contrast = c("group", "Low_Severity", "High_Severity"))
res_Low_vs_High <- res_Low_vs_High[order(res_Low_vs_High$padj), ]

# Extract significant DEGs (strict filter)
deg <- subset(res_Low_vs_High, padj < 0.05 & abs(log2FoldChange) > 1)

# Export all DEGs
write.csv(deg, "DEG_Low_vs_High_strict.csv")
```
[Check out the file: DEG_Low_vs_High_strict.csv](DEG_Low_vs_High_strict.csv)


# 6. Separate up- and down-regulated DEGs
```r
deg_df <- as.data.frame(deg)
deg_df$gene_id <- rownames(deg_df)

deg_up <- deg_df %>% filter(log2FoldChange > 1)
deg_down <- deg_df %>% filter(log2FoldChange < -1)

write.csv(deg_up, "Up_DEG_Low_vs_High_strict.csv", row.names = FALSE)
write.csv(deg_down, "Down_DEG_Low_vs_High_strict.csv", row.names = FALSE)
```

[Check out the file: Up_DEG_Low_vs_High_strict.csv](Up_DEG_Low_vs_High_strict.csv)
[Check out the file: Down_DEG_Low_vs_High_strict.csv](Down_DEG_Low_vs_High_strict.csv)


## 7. Load annotation file
```r
annotation <- read_delim("fullAnnotation.tsv.txt", delim = "\t", col_types = cols())

# Fix duplicated column names if needed
colnames(annotation) <- make.unique(colnames(annotation))
```

### Add human-readable COG names
```r
cog_dict <- c(
  "C" = "Energy production and conversion",
  "D" = "Cell cycle control, cell division, chromosome partitioning",
  "E" = "Amino acid transport and metabolism",
  "F" = "Nucleotide transport and metabolism",
  "G" = "Carbohydrate transport and metabolism",
  "H" = "Coenzyme transport and metabolism",
  "I" = "Lipid transport and metabolism",
  "J" = "Translation, ribosomal structure and biogenesis",
  "K" = "Transcription",
  "L" = "Replication, recombination and repair",
  "M" = "Cell wall/membrane/envelope biogenesis",
  "N" = "Cell motility",
  "O" = "Posttranslational modification, protein turnover, chaperones",
  "P" = "Inorganic ion transport and metabolism",
  "Q" = "Secondary metabolites biosynthesis, transport and catabolism",
  "R" = "General function prediction only",
  "S" = "Function unknown",
  "T" = "Signal transduction mechanisms",
  "U" = "Intracellular trafficking, secretion, vesicular transport",
  "V" = "Defense mechanisms",
  "W" = "Extracellular structures",
  "Y" = "Nuclear structure",
  "Z" = "Cytoskeleton"
)

# Add readable COG name
annotation <- annotation %>%
  mutate(COG_name = cog_dict[COG_category])
```

## 6. Annotate DEGs with COG
```r
deg_up <- as_tibble(deg_up)
deg_down <- as_tibble(deg_down)

deg_up_annot <- deg_up %>%
  left_join(dplyr::select(annotation, gene_id, COG_category, COG_name), by = "gene_id")
deg_down_annot <- deg_down %>%
  left_join(dplyr::select(annotation, gene_id, COG_category, COG_name), by = "gene_id")
```

# 7. Summarize functional categories
```r
# Universe
universe_summary <- annotation %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Universe = n_distinct(gene_id))

# Up DEGs
up_summary <- deg_up_annot %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Up = n_distinct(gene_id))

# Down DEGs
down_summary <- deg_down_annot %>%
  filter(!is.na(COG_name)) %>%
  group_by(COG_name) %>%
  summarise(Down = n_distinct(gene_id))
```

### Merge summaries
```r
summary_table <- universe_summary %>%
  full_join(up_summary, by = "COG_name") %>%
  full_join(down_summary, by = "COG_name") %>%
  replace(is.na(.), 0) %>%
  arrange(desc(Universe))

write_csv(summary_table, "Summary_by_COG_categories_Low_vs_High.csv")
print(n=21, summary_table)
```
