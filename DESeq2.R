if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

setwd("INSERT PATH TO RESPOSITORY HERE")

# Installing DESeq2

BiocManager::install("DESeq2")

library(DESeq2)

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

# Condition data
condition_data_10 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts10_norm)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts10_norm)))))))
condition_data_30 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts30_norm)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts30_norm)))))))
condition_data_120 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts120_norm)))), 
                                                      rep("Dex", sum(grepl("Dex", colnames(counts120_norm)))))))

# DESeq2 Function (with integer check)
func_DESeq2 <- function(count_data, condition_data) {
  # Exclude gene labels and convert to matrix
  count_data_matrix <- as.matrix(count_data[-1])  # Exclude first column (gene labels)
  
  # Debug: Print dimensions and column names
  cat("Dimensions of count_data_matrix:", dim(count_data_matrix), "\n")
  cat("Number of rows in condition_data:", nrow(condition_data), "\n")
  cat("Column names in count_data_matrix:", colnames(count_data_matrix), "\n")
  
  # Ensure the count_data and condition_data have the same number of columns and rows
  if (ncol(count_data_matrix) != nrow(condition_data)) {
    stop("The number of columns in count_data does not match the number of rows in condition_data")
  }
  
  # Check if count_data_matrix contains only integers
  if (!all(count_data_matrix == round(count_data_matrix))) {
    cat("Warning: Non-integer values detected in count_data. Converting to integers.\n")
    count_data_matrix <- round(count_data_matrix)
  }
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = count_data_matrix,
                                colData = condition_data,
                                design = ~ condition)
  
  dds <- DESeq(dds, fitType = "local")
  res <- results(dds)
  
  res <- as.data.frame(res)
  res$Gene <- count_data$Gene  # Restoring gene labels
  
  return(res)
}

# Running DESeq2 (only for .norm datasets)
count_datasets <- list(
  "counts10_norm" = counts10_norm,
  "counts30_norm" = counts30_norm,
  "counts120_norm" = counts120_norm
)

norm_results_list <- list()

for (name in names(count_datasets)) {
  cat("Running DESeq2 for dataset", name, "\n")
  condition_data <- switch(name,
                           "counts10_norm" = condition_data_10,
                           "counts30_norm" = condition_data_30,
                           "counts120_norm" = condition_data_120)
  
  # Running DESeq2 with dimension checks
  result <- func_DESeq2(count_datasets[[name]], condition_data)
  norm_results_list[[name]] <- result
}




# FOR RAW DATA W/ DESEQ2 INTERNALIZATION

library(DESeq2)


# Condition data
condition_data_10 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts10_raw)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts10_raw)))))))
condition_data_30 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts30_raw)))), 
                                                     rep("Dex", sum(grepl("Dex", colnames(counts30_raw)))))))
condition_data_120 <- data.frame(condition = factor(c(rep("Control", sum(grepl("crtl", colnames(counts120_raw)))), 
                                                      rep("Dex", sum(grepl("Dex", colnames(counts120_raw)))))))

# DESeq2 Function (with integer check)
func_DESeq2 <- function(count_data, condition_data) {
  # Exclude gene labels and convert to matrix
  count_data_matrix <- as.matrix(count_data[-1])  # Exclude first column (gene labels)
  
  # Debug: Print dimensions and column names
  cat("Dimensions of count_data_matrix:", dim(count_data_matrix), "\n")
  cat("Number of rows in condition_data:", nrow(condition_data), "\n")
  cat("Column names in count_data_matrix:", colnames(count_data_matrix), "\n")
  
  # Ensure the count_data and condition_data have the same number of columns and rows
  if (ncol(count_data_matrix) != nrow(condition_data)) {
    stop("The number of columns in count_data does not match the number of rows in condition_data")
  }
  
  # Check if count_data_matrix contains only integers
  if (!all(count_data_matrix == round(count_data_matrix))) {
    cat("Warning: Non-integer values detected in count_data. Converting to integers.\n")
    count_data_matrix <- round(count_data_matrix)
  }
  
  # Create DESeqDataSet with internal normalization
  dds <- DESeqDataSetFromMatrix(countData = count_data_matrix,
                                colData = condition_data,
                                design = ~ condition)
  
  dds <- DESeq(dds)  # This will apply internal normalization
  res <- results(dds)
  
  res <- as.data.frame(res)
  res$Gene <- count_data$Gene  # Restoring gene labels
  
  return(res)
}

# Running DESeq2 for raw counts datasets
count_datasets <- list(
  "counts10_raw" = counts10_raw,
  "counts30_raw" = counts30_raw,
  "counts120_raw" = counts120_raw
)

raw_results_list <- list()

for (name in names(count_datasets)) {
  cat("Running DESeq2 for dataset", name, "\n")
  condition_data <- switch(name,
                           "counts10_raw" = condition_data_10,
                           "counts30_raw" = condition_data_30,
                           "counts120_raw" = condition_data_120)
  
  # Running DESeq2 with dimension checks
  result <- func_DESeq2(count_datasets[[name]], condition_data)
  raw_results_list[[name]] <- result
}
