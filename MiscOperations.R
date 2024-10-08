install.packages("writexl")

# Sort DESeq2 tables by p-value
sorted_DESeq2_10 <- norm_results_list[["counts10_norm"]][order(results_list[["counts10_norm"]]$pvalue), ]
sorted_DESeq2_30 <- norm_results_list[["counts30_norm"]][order(results_list[["counts30_norm"]]$pvalue), ]
sorted_DESeq2_120 <- norm_results_list[["counts120_norm"]][order(results_list[["counts120_norm"]]$pvalue), ]

# Sort limma voom results by p-value (assuming column "P.Value")
sorted_limma_10 <- limma_results_10[order(limma_results_10$P.Value), ]
sorted_limma_30 <- limma_results_30[order(limma_results_30$P.Value), ]
sorted_limma_120 <- limma_results_120[order(limma_results_120$P.Value), ]

# Export to excel

library(writexl)

write_xlsx(list("DESeq2_10min" = sorted_DESeq2_10,
                "limma_10min" = sorted_limma_10,
                "DESeq2_30min" = sorted_DESeq2_30,
                "limma_30min" = sorted_limma_30,
                "DESeq2_120min" = sorted_DESeq2_120,
                "limma_120min" = sorted_limma_120),
           path = "sorted_combined_results.xlsx")


# Percentage Calculations

# Define DESeq2 results list
DESeq2_Percent <- list(
  counts10 = norm_results_list[["counts10_norm"]],
  counts30 = norm_results_list[["counts30_norm"]],
  counts120 = norm_results_list[["counts120_norm"]]
)

# Function to calculate percentage of genes with p-value < 0.05
calculate_percentage <- function(results) {
  total_genes <- nrow(results)  # Total number of genes
  # Count significant genes, excluding NA values
  significant_genes <- sum(results$pvalue < 0.05, na.rm = TRUE)  
  percentage <- (significant_genes / total_genes) * 100  # Calculate percentage
  return(percentage)
}

# Calculate percentages for each time point
percentages_DESeq2 <- sapply(DESeq2_Percent, calculate_percentage)

# Print results
cat("Percentage of significant genes in DESeq2 results:\n")
for (name in names(percentages_DESeq2)) {
  cat(paste0(name, ": ", round(percentages_DESeq2[name], 2), "%\n"))
