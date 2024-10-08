setwd("INSERT PATH TO RESPOSITORY HERE")

# Installing Limma

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("limma", force = TRUE)


# CSV insertion

# 10 minute counts CSV organization
counts10 <- read.csv("(CSV) 10min Dex vs control_count matrix.csv", stringsAsFactors = FALSE, skip = 1)
gene_labels10 <- counts10$Gene

counts10_raw <- counts10[, grepl("crtl|Dex", colnames(counts10)) & !grepl("\\.norm|\\.rlog", colnames(counts10))]
counts10_raw <- cbind(Gene = gene_labels10, counts10_raw)

counts10_norm <- counts10[, grepl("crtl_\\w+\\.norm|Dex_\\w+\\.norm", colnames(counts10))]
counts10_norm <- cbind(Gene = gene_labels10, counts10_norm)

counts10_rlog <- counts10[, grepl("crtl_\\w+\\.rlog|Dex_\\w+\\.rlog", colnames(counts10))]
counts10_rlog <- cbind(Gene = gene_labels10, counts10_rlog)

# 30 minute counts CSV organization
counts30 <- read.csv("(CSV) 30min Dex vs control_count matrix.csv", stringsAsFactors = FALSE, skip = 1)
gene_labels30 <- counts30$Gene

counts30_raw <- counts30[, grepl("crtl_\\w+|Dex_\\w+", colnames(counts30)) & !grepl("\\.norm|\\.rlog", colnames(counts30))]
counts30_raw <- cbind(Gene = gene_labels30, counts30_raw)

counts30_norm <- counts30[, grepl("crtl_\\w+\\.norm|Dex_\\w+\\.norm", colnames(counts30))]
counts30_norm <- cbind(Gene = gene_labels30, counts30_norm)

counts30_rlog <- counts30[, grepl("crtl_\\w+\\.rlog|Dex_\\w+\\.rlog", colnames(counts30))]
counts30_rlog <- cbind(Gene = gene_labels30, counts30_rlog)

# 120 minute counts CSV organization
counts120 <- read.csv("(CSV) 120min Dex vs control_count matrix.csv", stringsAsFactors = FALSE, skip = 1)
gene_labels120 <- counts120$Gene

counts120_raw <- counts120[, grepl("crtl_\\w+|Dex_\\w+", colnames(counts120)) & !grepl("\\.norm|\\.rlog", colnames(counts120))]
counts120_raw <- cbind(Gene = gene_labels120, counts120_raw)

counts120_norm <- counts120[, grepl("crtl_\\w+\\.norm|Dex_\\w+\\.norm", colnames(counts120))]
counts120_norm <- cbind(Gene = gene_labels120, counts120_norm)

counts120_rlog <- counts120[, grepl("crtl_\\w+\\.rlog|Dex_\\w+\\.rlog", colnames(counts120))]
counts120_rlog <- cbind(Gene = gene_labels120, counts120_rlog)

# Running Limma

if (!requireNamespace("limma", quietly = TRUE)) {
  install.packages("limma")
}

library(limma)

# Limma voom function
func_limma_voom <- function(norm_counts, condition_data) {
  
  rownames(norm_counts) <- make.unique(norm_counts$Gene)  
  norm_counts <- norm_counts[, -1]  # Removing Gene Column
  
  design <- model.matrix(~ condition_data$condition)
  v <- voom(as.matrix(norm_counts), design, plot = TRUE)
  
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  res <- topTable(fit, coef = 2, number = Inf) 
  res$Gene <- rownames(res)  # Restoring Gene Column
  return(res)
}

# Condition data
condition_data_10 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts10_norm)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts10_norm)))))))
condition_data_30 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts30_norm)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts30_norm)))))))
condition_data_120 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts120_norm)))), 
                                                      rep("Dex", sum(grepl("Dex", colnames(counts120_norm)))))))

# Running limma voom
limma_results_10 <- func_limma_voom(counts10_norm, condition_data_10)
limma_results_30 <- func_limma_voom(counts30_norm, condition_data_30)
limma_results_120 <- func_limma_voom(counts120_norm, condition_data_120)

# Check results
head(limma_results_10)
head(limma_results_30)
head(limma_results_120)
