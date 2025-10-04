# Install required packages if not already installed
if (!requireNamespace("pegas", quietly = TRUE)) install.packages("pegas")
if (!requireNamespace("ape", quietly = TRUE)) install.packages("ape")

library(pegas)
library(ape)

# 1. Define genes and orders
genes <- c("mcyA", "mcyB", "mcyE", "mcyH")
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
base_path <- "E:/Research/Research/MSA"

# 2. Function to compute stats for one gene/order
compute_stats <- function(order, gene) {
  fasta_file <- file.path(base_path, order, paste0(tolower(gene), "_aligned.fasta"))
  if (!file.exists(fasta_file)) return(NULL)
  
  # Read alignment
  gen <- tryCatch(read.dna(fasta_file, format="fasta"), error=function(e) NULL)
  if (is.null(gen)) return(NULL)
  
  # Alignment length
  alignment_length <- tryCatch(ncol(gen), error=function(e) NA)
  
  # Polymorphic sites
  poly_sites <- tryCatch(length(seg.sites(gen)), error=function(e) NA)
  
  # Nucleotide diversity
  pi <- tryCatch(nuc.div(gen), error=function(e) NA)
  
  # Tajima's D
  tajima <- tryCatch(tajima.test(gen)$D, error=function(e) NA)
  
  # Mutation rate per base
  mutation_rate <- if (!is.na(alignment_length) && alignment_length > 0 && !is.na(poly_sites)) poly_sites / alignment_length else NA
  
  stats <- data.frame(
    Gene = gene,
    Order = order,
    PolymorphicSites = poly_sites,
    NucleotideDiversity = pi,
    TajimasD = tajima,
    MutationRatePerBase = mutation_rate
  )
  
  # Write per-gene/order TSV
  out_file <- file.path(base_path, order, paste0(tolower(gene), "_stats.tsv"))
  write.table(stats, out_file, sep="\t", row.names=FALSE, quote=FALSE)
  
  return(stats)
}

# 3. Loop over all genes and orders, collect stats
all_stats <- do.call(rbind, lapply(genes, function(gene) {
  do.call(rbind, lapply(orders, function(order) {
    compute_stats(order, gene)
  }))
}))

# 4. Write combined summary table
summary_file <- file.path(base_path, "plots", "gene_order_stats_table.tsv")
if (!dir.exists(file.path(base_path, "plots"))) dir.create(file.path(base_path, "plots"), recursive=TRUE)
write.table(all_stats, summary_file, sep="\t", row.names=FALSE, quote=FALSE)

cat("Summary table written to:", summary_file, "\n")