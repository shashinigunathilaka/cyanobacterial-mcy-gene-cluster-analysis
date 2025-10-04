***********************************************************************mcyA***********************************************************************
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# 1. List the three files for mcyA
matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
gene <- "mcyA"
files <- file.path(matrix_dir, paste0(orders, "_", gene, "_mutation_matrix.csv"))

# 2. Read and combine, adding Order info
all_mats <- lapply(seq_along(files), function(i) {
  f <- files[i]
  if (!file.exists(f)) stop("File not found: ", f)
  mat <- read_csv(f, show_col_types = FALSE)
  mat$Order <- orders[i]
  mat
})
combined <- bind_rows(all_mats)

# 3. Set rownames and remove SeqHeader for PCA
rownames(combined) <- paste(combined$Order, combined$SeqHeader, sep="_")
mat <- combined %>% select(starts_with("Pos"))

# 4. Keep only columns present in all orders (no NA)
mat <- mat[, colSums(is.na(mat)) == 0, drop=FALSE]

# 5. Remove zero-variance columns
mat <- mat[, apply(mat, 2, function(x) var(x, na.rm=TRUE) > 0), drop=FALSE]

# 6. Check if enough columns remain
if (ncol(mat) < 2) stop("Not enough variable columns for PCA.")

# 7. PCA
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
scores <- as.data.frame(pca$x)
scores$Order <- combined$Order
scores$SeqHeader <- combined$SeqHeader

# 8. Plot PC1 vs PC2, color by Order, add ellipses
p <- ggplot(scores, aes(x=PC1, y=PC2, color=Order, fill=Order)) +
  stat_ellipse(geom="polygon", alpha=0.2, level=0.95, show.legend=FALSE) +
  geom_point(size=2) +
  labs(title="PCA of mcyA Mutation Matrix (All Orders)",
       x="PC1", y="PC2") +
  theme_minimal(base_size=12)

# 9. Save the plot
plot_file <- file.path(plot_dir, "mcyA_PCA_all_orders.png")
ggsave(plot_file, p, width=8, height=6, dpi=300)
cat("PCA plot saved to:", plot_file, "\n")

# 10. Save PC1-PC5 table
pc_table <- scores %>%
  select(Order, SeqHeader, PC1, PC2, PC3, PC4, PC5)
pc_file <- file.path(matrix_dir, "mcyA_PCs1-5.csv")
write_csv(pc_table, pc_file)
cat("PC table saved to:", pc_file, "\n")    

***********************************************************************mcyB***********************************************************************

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# 1. List the three files for mcyB
matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
gene <- "mcyB"
files <- file.path(matrix_dir, paste0(orders, "_", gene, "_mutation_matrix.csv"))

# 2. Read and combine, adding Order info
all_mats <- lapply(seq_along(files), function(i) {
  f <- files[i]
  if (!file.exists(f)) stop("File not found: ", f)
  mat <- read_csv(f, show_col_types = FALSE)
  mat$Order <- orders[i]
  mat
})
combined <- bind_rows(all_mats)

# 3. Set rownames and remove SeqHeader for PCA
rownames(combined) <- paste(combined$Order, combined$SeqHeader, sep="_")
mat <- combined %>% select(starts_with("Pos"))

# 4. Keep only columns present in all orders (no NA)
mat <- mat[, colSums(is.na(mat)) == 0, drop=FALSE]

# 5. Remove zero-variance columns
mat <- mat[, apply(mat, 2, function(x) var(x, na.rm=TRUE) > 0), drop=FALSE]

# 6. Check if enough columns remain
if (ncol(mat) < 2) stop("Not enough variable columns for PCA.")

# 7. PCA
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
scores <- as.data.frame(pca$x)
scores$Order <- combined$Order
scores$SeqHeader <- combined$SeqHeader

# 8. Plot PC1 vs PC2, color by Order, add ellipses
p <- ggplot(scores, aes(x=PC1, y=PC2, color=Order, fill=Order)) +
  stat_ellipse(geom="polygon", alpha=0.2, level=0.95, show.legend=FALSE) +
  geom_point(size=2) +
  labs(title="PCA of mcyB Mutation Matrix (All Orders)",
       x="PC1", y="PC2") +
  theme_minimal(base_size=12)

# 9. Save the plot
plot_file <- file.path(plot_dir, "mcyB_PCA_all_orders.png")
ggsave(plot_file, p, width=8, height=6, dpi=300)
cat("PCA plot saved to:", plot_file, "\n")

# 10. Save PC1-PC5 table
pc_table <- scores %>%
  select(Order, SeqHeader, PC1, PC2, PC3, PC4, PC5)
pc_file <- file.path(matrix_dir, "mcyB_PCs1-5.csv")
write_csv(pc_table, pc_file)
cat("PC table saved to:", pc_file, "\n")

***********************************************************************mcyE***********************************************************************

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# 1. List the three files for mcyE
matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
gene <- "mcyE"
files <- file.path(matrix_dir, paste0(orders, "_", gene, "_mutation_matrix.csv"))

# 2. Read and combine, adding Order info
all_mats <- lapply(seq_along(files), function(i) {
  f <- files[i]
  if (!file.exists(f)) stop("File not found: ", f)
  mat <- read_csv(f, show_col_types = FALSE)
  mat$Order <- orders[i]
  mat
})
combined <- bind_rows(all_mats)

# 3. Set rownames and remove SeqHeader for PCA
rownames(combined) <- paste(combined$Order, combined$SeqHeader, sep="_")
mat <- combined %>% select(starts_with("Pos"))

# 4. Keep only columns present in all orders (no NA)
mat <- mat[, colSums(is.na(mat)) == 0, drop=FALSE]

# 5. Remove zero-variance columns
mat <- mat[, apply(mat, 2, function(x) var(x, na.rm=TRUE) > 0), drop=FALSE]

# 6. Check if enough columns remain
if (ncol(mat) < 2) stop("Not enough variable columns for PCA.")

# 7. PCA
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
scores <- as.data.frame(pca$x)
scores$Order <- combined$Order
scores$SeqHeader <- combined$SeqHeader

# 8. Plot PC1 vs PC2, color by Order, add ellipses
p <- ggplot(scores, aes(x=PC1, y=PC2, color=Order, fill=Order)) +
  stat_ellipse(geom="polygon", alpha=0.2, level=0.95, show.legend=FALSE) +
  geom_point(size=2) +
  labs(title="PCA of mcyE Mutation Matrix (All Orders)",
       x="PC1", y="PC2") +
  theme_minimal(base_size=12)

# 9. Save the plot
plot_file <- file.path(plot_dir, "mcyE_PCA_all_orders.png")
ggsave(plot_file, p, width=8, height=6, dpi=300)
cat("PCA plot saved to:", plot_file, "\n")

# 10. Save PC1-PC5 table
pc_table <- scores %>%
  select(Order, SeqHeader, PC1, PC2, PC3, PC4, PC5)
pc_file <- file.path(matrix_dir, "mcyE_PCs1-5.csv")
write_csv(pc_table, pc_file)
cat("PC table saved to:", pc_file, "\n")

***********************************************************************mcyH***********************************************************************

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# 1. List the three files for mcyH
matrix_dir <- "E:/Research/Research/MSA/feature_matrices"
plot_dir <- "E:/Research/Research/MSA/plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
orders <- c("Chroococcales", "Nostocales", "Oscillatoriales")
gene <- "mcyH"
files <- file.path(matrix_dir, paste0(orders, "_", gene, "_mutation_matrix.csv"))

# 2. Read and combine, adding Order info
all_mats <- lapply(seq_along(files), function(i) {
  f <- files[i]
  if (!file.exists(f)) stop("File not found: ", f)
  mat <- read_csv(f, show_col_types = FALSE)
  mat$Order <- orders[i]
  mat
})
combined <- bind_rows(all_mats)

# 3. Set rownames and remove SeqHeader for PCA
rownames(combined) <- paste(combined$Order, combined$SeqHeader, sep="_")
mat <- combined %>% select(starts_with("Pos"))

# 4. Keep only columns present in all orders (no NA)
mat <- mat[, colSums(is.na(mat)) == 0, drop=FALSE]

# 5. Remove zero-variance columns
mat <- mat[, apply(mat, 2, function(x) var(x, na.rm=TRUE) > 0), drop=FALSE]

# 6. Check if enough columns remain
if (ncol(mat) < 2) stop("Not enough variable columns for PCA.")

# 7. PCA
pca <- prcomp(mat, center=TRUE, scale.=TRUE)
scores <- as.data.frame(pca$x)
scores$Order <- combined$Order
scores$SeqHeader <- combined$SeqHeader

# 8. Plot PC1 vs PC2, color by Order, add ellipses
p <- ggplot(scores, aes(x=PC1, y=PC2, color=Order, fill=Order)) +
  stat_ellipse(geom="polygon", alpha=0.2, level=0.95, show.legend=FALSE) +
  geom_point(size=2) +
  labs(title="PCA of mcyH Mutation Matrix (All Orders)",
       x="PC1", y="PC2") +
  theme_minimal(base_size=12)

# 9. Save the plot
plot_file <- file.path(plot_dir, "mcyH_PCA_all_orders.png")
ggsave(plot_file, p, width=8, height=6, dpi=300)
cat("PCA plot saved to:", plot_file, "\n")

# 10. Save PC1-PC5 table
pc_table <- scores %>%
  select(Order, SeqHeader, PC1, PC2, PC3, PC4, PC5)
pc_file <- file.path(matrix_dir, "mcyH_PCs1-5.csv")
write_csv(pc_table, pc_file)
cat("PC table saved to:", pc_file, "\n")