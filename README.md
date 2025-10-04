# Cyanobacterial *mcy* Gene Cluster Analysis

## Overview

This project provides a robust pipeline for analyzing population genetics and mutation patterns in cyanobacterial *mcy* gene clusters across different taxonomic orders. It integrates C# and R scripts to automate data extraction, alignment, feature matrix construction, statistical analysis, and dimensionality reduction (PCA), supporting downstream visualization and interpretation.

---

## Features

- **Automated extraction and alignment of gene sequences**
- **Feature matrix construction for mutation analysis**
- **Calculation of population genetics statistics (polymorphic sites, nucleotide diversity, Tajima's D, mutation rate per base)**
- **Dimensionality reduction and clustering (PCA)**
- **Comprehensive summary tables and publication-ready plots**

---

## Directory Structure
```
E:/Research/Research/
??? DataSet/                 # Raw sequence data
??? ReferenceGenes/          # Reference gene sequences
??? ExtractedGenes/          # Extracted gene FASTA files
??? MSA/                     # Multiple sequence alignments and stats
?   ??? feature_matrices/    # Mutation feature matrices (CSV)
?   ??? plots/               # Output plots and summary tables
?   ??? [Order]/             # Per-order alignments and stats (replace [Order] with actual order names)


```

---

## Installation

### Prerequisites

- **.NET 8 SDK** (for C# components)
- **R (? 4.0)** with the following packages:
  - `pegas`
  - `ape`
  - `readr`
  - `dplyr`
  - `ggplot2`
  - `stringr`
  - `Biostrings` (for FASTA filtering, if needed)

### R Package Installation
```
install.packages(c("pegas", "ape", "readr", "dplyr", "ggplot2", "stringr", "Biostrings"))
```

---

## Usage

### 1. Sequence Extraction and Alignment (C#)

- Extract gene sequences from raw data and reference sets.
- Align sequences using external tools (e.g., MUSCLE, MAFFT).
- Store aligned FASTA files in `MSA/[Order]/`.

### 2. Feature Matrix Construction (C#)

- Generate binary mutation matrices for each gene and order.
- Output CSV files to `MSA/feature_matrices/`.

### 3. Population Genetics Statistics (R)

- Run the provided R script to calculate:
  - Polymorphic sites
  - Nucleotide diversity
  - Tajima's D
  - Mutation rate per base

**Example R script:**
```R
library(pegas)
library(ape)

genes <- c("mcyA", "mcyB", "mcyE", "mcyH")
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
base_path <- "E:/Research/Research/MSA"

compute_stats <- function(order, gene) {
  fasta_file <- file.path(base_path, order, paste0(tolower(gene), "_aligned.fasta"))
  
  if (!file.exists(fasta_file)) return(NULL)
  
  gen <- tryCatch(
    read.dna(fasta_file, format = "fasta"),
    error = function(e) NULL
  )
  if (is.null(gen)) return(NULL)
  
  alignment_length <- tryCatch(ncol(gen), error = function(e) NA)
  poly_sites <- tryCatch(length(seg.sites(gen)), error = function(e) NA)
  pi <- tryCatch(nuc.div(gen), error = function(e) NA)
  tajima <- tryCatch(tajima.test(gen)$D, error = function(e) NA)
  
  mutation_rate <- if (!is.na(alignment_length) && alignment_length > 0 && !is.na(poly_sites)) {
    poly_sites / alignment_length
  } else {
    NA
  }
  
  stats <- data.frame(
    Gene = gene,
    Order = order,
    PolymorphicSites = poly_sites,
    NucleotideDiversity = pi,
    TajimasD = tajima,
    MutationRatePerBase = mutation_rate
  )
  
  out_file <- file.path(base_path, order, paste0(tolower(gene), "_stats.tsv"))
  write.table(stats, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(stats)
}

all_stats <- do.call(
  rbind,
  lapply(genes, function(gene) {
    do.call(
      rbind,
      lapply(orders, function(order) {
        compute_stats(order, gene)
      })
    )
  })
)

summary_file <- file.path(base_path, "plots", "gene_order_stats_table.tsv")

if (!dir.exists(file.path(base_path, "plots"))) {
  dir.create(file.path(base_path, "plots"), recursive = TRUE)
}

write.table(all_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Summary table written to:", summary_file, "\n")
```

### 4. Statistics Aggregation (C#)

- Use `GeneOrderStatsTableBuilder.cs` to aggregate per-gene/order stats into a summary table for reporting and visualization.

### 5. PCA and Clustering (R)

- Perform PCA on mutation feature matrices for each gene.
- Visualize PC1 vs PC2, colored by order, with confidence ellipses.
- Save principal component tables and plots for downstream analysis.

**Example PCA R code:**
```R
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
plot_dir <- "E:/Research/Research/MSA/plots"
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
gene <- "mcyA"

files <- file.path(matrix_dir, paste0(orders, "_", gene, "mutation_matrix.csv"))

all_mats <- lapply(seq_along(files), function(i) {
  f <- files[i]
  mat <- read_csv(f, show_col_types = FALSE)
  mat$Order <- orders[i]
  mat
})

combined <- bind_rows(all_mats)
rownames(combined) <- paste(combined$Order, combined$SeqHeader, sep = "")

mat <- combined %>%
  select(starts_with("Pos"))

mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]
mat <- mat[, apply(mat, 2, function(x) var(x, na.rm = TRUE) > 0), drop = FALSE]

if (ncol(mat) < 2) {
  stop("Not enough variable columns for PCA.")
}

pca <- prcomp(mat, center = TRUE, scale. = TRUE)
scores <- as.data.frame(pca$x)
scores$Order <- combined$Order
scores$SeqHeader <- combined$SeqHeader

p <- ggplot(scores, aes(x = PC1, y = PC2, color = Order, fill = Order)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE) +
  geom_point(size = 2) +
  labs(
    title = "PCA of mcyA Mutation Matrix (All Orders)",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal(base_size = 12)

plot_file <- file.path(plot_dir, "mcyA_PCA_all_orders.png")
ggsave(plot_file, p, width = 8, height = 6, dpi = 300)

pc_table <- scores %>%
  select(Order, SeqHeader, PC1, PC2, PC3, PC4, PC5)

pc_file <- file.path(matrix_dir, "mcyA_PCs1-5.csv")
write_csv(pc_table, pc_file)

```

---

## Output

- **Summary tables:** Per-gene/order statistics in TSV and CSV format.
- **Feature matrices:** Binary mutation matrices for each gene/order.
- **PCA plots:** Visualizations of genetic structure and clustering.
- **Principal component tables:** For further clustering or downstream analysis.

---

## Troubleshooting

- Ensure all sequences in each alignment file are the same length (including gaps).
- Check that all required R packages are installed and up to date.
- Review log and error messages for missing files or data inconsistencies.

---

## License

This project is provided for academic research purposes.  
Please cite appropriately if used in publications.

---

## Contact

For questions or contributions, please contact the project maintainer.

