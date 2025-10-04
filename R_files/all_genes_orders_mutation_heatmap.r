library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# 1. Set directory and list all CSVs
matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
csv_files <- list.files(matrix_dir, pattern="*_mutation_matrix.csv", full.names=TRUE)

# 2. Read and combine all CSVs
all_long <- do.call(rbind, lapply(csv_files, function(f) {
  # Extract gene and order from filename
  fname <- basename(f)
  m <- str_match(fname, "^([A-Za-z]+)_([A-Za-z]+)_mutation_matrix\\.csv$")
  if (is.na(m[1,1])) return(NULL)
  order <- m[1,2]
  gene <- m[1,3]
  mat <- read.csv(f, row.names=1, check.names=FALSE)
  mat$SeqHeader <- rownames(mat)
  long <- mat %>%
    pivot_longer(-SeqHeader, names_to="Position", values_to="Mutation")
  long$Gene <- gene
  long$Order <- order
  long
}))

# 3. Make factors for better ordering in plot
all_long$Gene <- factor(all_long$Gene, levels=c("mcyA", "mcyB", "mcyE", "mcyH"))
all_long$Order <- factor(all_long$Order, levels=c("Chroococcales", "Nostocales", "Oscillatoriales"))

# 3b. Make gene names expressions for italic facet labels
all_long$Gene <- recode(
  as.character(all_long$Gene),
  "mcyA" = "italic(mcyA)",
  "mcyB" = "italic(mcyB)",
  "mcyE" = "italic(mcyE)",
  "mcyH" = "italic(mcyH)"
)

# 4. Plot: Heatmap, faceted by gene, y=Order, x=Position, color=Mutation
p <- ggplot(all_long, aes(x=Position, y=Order, fill=factor(Mutation))) +
  geom_tile(color=NA) +
  scale_fill_manual(values=c("0"="white", "1"="red"), name="Mutation", labels=c("No", "Yes")) +
  facet_wrap(~ Gene, ncol=1, scales="free_x", labeller = label_parsed) + # Italic facet labels
  labs(title="Mutation Presence/Absence Across Genes and Orders",
       x="Alignment Position", y="Order") +
  theme_minimal(base_size=12) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size=14),
        legend.position="top")

# 5. Save the plot
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
heatmap_file <- file.path(plot_dir, "all_genes_orders_mutation_heatmap.png")
ggsave(heatmap_file, p, width=12, height=6, dpi=300)

cat("Heatmap saved to:", heatmap_file, "\n")  