library(reshape2)
library(ggplot2)
library(dplyr)

# Read the summary CSV
data <- read.csv("E:/Research/Research/MSA/feature_matrices/mutation_type_summary.csv")

# List of mutation type columns (adjust if your column names differ)
mutation_type_cols <- c("Missense", "Frameshift", "FrameshiftDeletion", "FrameshiftInsertion", "Complex", "InFrame")

# Ensure these columns are numeric
data[mutation_type_cols] <- lapply(data[mutation_type_cols], as.numeric)

# Calculate total mutations per row (gene/order)
data$total <- rowSums(data[mutation_type_cols])

# Calculate percentages for each mutation type
for (col in mutation_type_cols) {
    data[[paste0(col, "_pct")]] <- 100 * data[[col]] / data$total
}

# Melt for ggplot (long format)
df_long <- melt(
    data,
    id.vars = c("Gene", "Order"),
    measure.vars = paste0(mutation_type_cols, "_pct"),
    variable.name = "MutationType",
    value.name = "Percent"
)

# Clean up MutationType names
df_long$MutationType <- gsub("_pct", "", df_long$MutationType)

# Make mutation type names more readable
df_long$MutationType <- recode(
    df_long$MutationType,
    "Missense" = "Missense",
    "Frameshift" = "Frameshift",
    "FrameshiftDeletion" = "Frameshift Deletion",
    "FrameshiftInsertion" = "Frameshift Insertion",
    "Complex" = "Complex",
    "InFrame" = "In-frame"
)

# Remove bars with 0% (or NA)
df_long <- subset(df_long, Percent > 0 & !is.na(Percent))

# Make gene names expressions for italic facet labels
df_long$Gene <- recode(
  df_long$Gene,
  "mcyA" = "italic(mcyA)",
  "mcyB" = "italic(mcyB)",
  "mcyE" = "italic(mcyE)",
  "mcyH" = "italic(mcyH)"
  # Add more if needed
)

# Calculate mean and SD for error bars
summary_df <- df_long %>%
    group_by(Gene, MutationType) %>%
    summarise(
        Mean = mean(Percent),
        SD = sd(Percent),
        N = n(),
        .groups = "drop"
    )

# Create plots directory if it doesn't exist
if (!dir.exists("E:/Research/Research/MSA/plots")) dir.create("E:/Research/Research/MSA/plots", recursive = TRUE)

# Plot with only the top error bar, y-axis always including 100, and italic facet labels
p <- ggplot(summary_df, aes(x=MutationType, y=Mean, fill=MutationType)) +
    geom_bar(stat="identity", width=0.4) +
    geom_errorbar(aes(ymin=Mean, ymax=Mean+SD), width=0.15, color="black", linewidth=0.7) +
    facet_wrap(~ Gene, ncol=2, labeller = label_parsed) +  # Italic facet labels
    scale_fill_brewer(palette="Set2") +
    scale_y_continuous(
        breaks=seq(0, 100, by=20),
        limits=c(0, 100),
        expand = expansion(mult = c(0, 0.05))
    ) +
    theme_minimal(base_size = 14) +
    labs(
        title="Mutation Type Percentages for mcy Genes",
        y="Percentage (%)",
        x=NULL,
        fill="Mutation Type"
    ) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        strip.text = element_text(face="bold"),
        plot.title = element_text(hjust=0.5)
    )

# Save the plot to the specified folder
ggsave("E:/Research/Research/MSA/plots/mutation_type_percentages.png", plot=p, width=8, height=8, dpi=300)
print(p)