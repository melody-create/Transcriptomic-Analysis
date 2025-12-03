# Load libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

setwd("C:/Users/")
getwd()

cts <- read.csv("cts.csv", row.names = 1)
coldata <- read.csv("coldata.csv", row.names = 1)

rownames(coldata)
colnames(cts)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition + age)

dds$condition <- relevel(dds$condition, ref = "WT")

dds <- DESeq(dds)
resultsNames(dds)

res_HT <- results(dds, contrast = c("condition", "HT", "WT"))
res_MUT <- results(dds, contrast = c("condition", "MUT", "WT"))

# find significant genes
sig_HT <- rownames(subset(res_HT, pvalue < 0.05 & abs(log2FoldChange) > 1))
sig_MUT <- rownames(subset(res_MUT, pvalue < 0.05 & abs(log2FoldChange) > 1))

sig_genes <- intersect(sig_HT, sig_MUT)
cat("Number of significant genes in BOTH HTvWT AND MUTvWT:", length(sig_genes), "\n")

vst_data <- vst(dds, blind = TRUE)
expr_matrix <- assay(vst_data)

# Subset to significant genes
expr_sig <- expr_matrix[sig_genes, ]

# Save that subset of significant genes
write.csv(expr_sig, "sig_gene_expression_log2.csv", quote = FALSE)

coldata <- read.csv("coldata.csv", row.names = 1)
coldata <- coldata[colnames(expr_sig), , drop = FALSE]

# Define colors
condition_colors <- c(WT = "blue", HT = "red", MUT = "orange")
age_col_fun <- c(three = "gray", five = "black")

top_anno <- HeatmapAnnotation(
  Condition = coldata$condition,
  Age = coldata$age,
  col = list(Condition = condition_colors,
             Age = age_col_fun),
  annotation_legend_param = list(
    Condition = list(title = "Condition"),
    Age = list(title = "Age")
  )
)

# define clustering
col_clust <- hclust(dist(t(expr_sig)), method = "complete")
row_clust <- hclust(dist(expr_sig), method = "complete")

col_fun_expr_gr <- colorRamp2(c(min(expr_sig), mean(expr_sig), max(expr_sig)),
                           c("white", "green", "darkgreen"))

ht_expr <- Heatmap(expr_sig,
                   name = "log2(Expression)",
                   top_annotation = top_anno,
                   col = col_fun_expr_pu,
                   cluster_rows = as.dendrogram(row_clust),
                   cluster_columns = as.dendrogram(col_clust),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   heatmap_legend_param = list(title = "log2(Expression)"))

draw(ht_expr)


################################################
#       AVERAGE
################################################

# create combined group
coldata <- coldata[colnames(expr_sig), , drop = FALSE]
group2 <- paste0(coldata$condition, "_age", coldata$age)   # e.g. "WT_age3"
groups2 <- unique(group2)

# choose aggregation function
agg_fun <- rowMeans 

expr_avg_cond_age <- do.call(cbind, lapply(groups2, function(g){
  cols <- which(group2 == g)
  if(length(cols) == 1) {
    expr_sig[, cols, drop = FALSE]
  } else {
    rowMeans(expr_sig[, cols, drop = FALSE])
  }
}))
colnames(expr_avg_cond_age) <- groups2
rownames(expr_avg_cond_age) <- rownames(expr_sig)

# save averaged matrix
write.csv(as.data.frame(expr_avg_cond_age),
          file = "expr_avg_by_condition_age.csv", row.names = TRUE, quote = FALSE)

mat_to_plot <- expr_avg_cond_age          

# Define colors
ann <- HeatmapAnnotation(
  Condition = cond,
  Age = age,
  col = list(
    Condition = c(WT = "blue", HT = "red", MUT = "orange"),
    Age = c(three = "gray", five = "black"))
  )

# clustering on averaged columns/rows
col_clust_avg <- hclust(dist(t(mat_to_plot)), method = "complete")
row_clust_avg <- hclust(dist(mat_to_plot), method = "complete")

# color functions
col_fun_avg <- colorRamp2(c(min(mat_to_plot), mean(mat_to_plot), max(mat_to_plot)),
                          c("white", "mediumorchid1", "purple4"))

# Heatmaps
ht_avg <- Heatmap(mat_to_plot,
                  name = "log2(Expression) (avg)",
                  top_annotation = ann,
                  col = col_fun_avg,
                  cluster_rows = as.dendrogram(row_clust_avg),
                  cluster_columns = as.dendrogram(col_clust_avg),
                  show_column_names = FALSE,
                  show_row_names = FALSE)

draw(ht_avg)
