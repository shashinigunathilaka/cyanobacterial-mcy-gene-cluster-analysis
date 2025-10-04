library(ggplot2)
library(dplyr)

# 1. Set output directory for plots
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 2. List of genes and orders
genes <- c("mcyA", "mcyB", "mcyE", "mcyH")
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
base_path <- "E:/Research/Research/MSA"

# 3. Read and combine all data
all_data <- do.call(rbind, lapply(genes, function(gene) {
do.call(rbind, lapply(orders, function(order) {
file <- file.path(base_path, order, paste0(tolower(gene), "_mutation_summary.tsv"))
if (file.exists(file)) {
df <- read.table(file, header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)
df <- df[!is.na(as.numeric(df$Mutations)), ]
df$Mutations <- as.numeric(df$Mutations)
df$Order <- order
df$Gene <- gene
df
}
}))
}))

# 4. Remove any NULLs (if a file was missing)
all_data <- all_data[!is.na(all_data$Gene), ]

# 4b. Make gene names expressions for italic facet labels
all_data$Gene <- recode(
all_data$Gene,
"mcyA" = "italic(mcyA)",
"mcyB" = "italic(mcyB)",
"mcyE" = "italic(mcyE)",
"mcyH" = "italic(mcyH)"
)

# 5. Faceted boxplot by gene, colored by order, with italic facet labels
p <- ggplot(all_data, aes(x=Order, y=Mutations, fill=Order)) +
geom_boxplot(outlier.shape=NA, alpha=0.7) +
facet_wrap(~ Gene, scales="free_y", labeller = label_parsed) + # Italic facet labels
labs(title="Mutations by Order and Gene (mcy genes)", y="Number of Mutations") +
theme_minimal() +
theme(axis.text.x = element_text(angle=45, hjust=1))

# 6. Save the plot
plot_file <- file.path(plot_dir, "mcy_mutations_by_order_and_gene.png")
ggsave(plot_file, p, width=10, height=6, dpi=300)

cat("Plot saved to:", plot_file, "\n")