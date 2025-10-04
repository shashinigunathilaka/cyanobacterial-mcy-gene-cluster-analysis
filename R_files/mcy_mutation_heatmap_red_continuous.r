library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Define genes, orders, and file paths
genes <- c("mcyA", "mcyB", "mcyE", "mcyH")
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
base_path <- "E:/Research/Research/MSA"

# 2. Read and combine all data
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

# 3. Calculate average mutations per gene/order
heatmap_data <- all_data %>%
group_by(Gene, Order) %>%
summarise(AvgMutations = mean(Mutations), .groups = "drop")

# 3b. Make gene names expressions for italic y-axis labels
heatmap_data$Gene <- recode(
heatmap_data$Gene,
"mcyA" = "italic(mcyA)",
"mcyB" = "italic(mcyB)",
"mcyE" = "italic(mcyE)",
"mcyH" = "italic(mcyH)"
)

# 4. Plot heatmap: continuous red gradient, no borders, no numbers, italic y-axis
p <- ggplot(heatmap_data, aes(x=Order, y=Gene, fill=AvgMutations)) +
geom_tile(color=NA) +
scale_fill_gradient(
low = "#fff5f0", # very light red
high = "#de2d26", # strong red
name = "Avg Mutations",
na.value = "grey90"
) +
scale_y_discrete(labels = parse(text = levels(factor(heatmap_data$Gene)))) +
labs(title="Average mcy Gene Mutations by Order",
x="Order", y="Gene") +
theme_minimal(base_size=14) +
theme(
axis.text.x = element_text(angle=45, hjust=1),
panel.grid = element_blank()
)

# 5. Save the heatmap
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
heatmap_file <- file.path(plot_dir, "mcy_mutation_heatmap_red_continuous.png")
ggsave(heatmap_file, p, width=6, height=4, dpi=300)

cat("Heatmap saved to:", heatmap_file, "\n")