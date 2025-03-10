# Load library
library(DESeq2)

# Load count data
count_data = read.csv("counts.filtered.csv", row.names = 1)

# Define metadata (patients and Helicobacter status)
patient <- c("ERR950185", "ERR950169", "ERR950170", "ERR950168")
Helicobactor <- c("Hel", "Hel", "none", "none")

# Create metadata DataFrame, rownames being patients and a column with if they had the bacteria or not
meta <- data.frame(Bacteria = Helicobactor, row.names = patient)


# Ensure the column order in count_data matches row names in meta
count_data = count_data[, colnames(count_data) %in% rownames(meta)]
meta <- meta[colnames(count_data), , drop=FALSE]



# Convert Bacteria column to a factor
meta$Bacteria <- factor(meta$Bacteria, levels = c("none", "Hel"))

# Create DESeq2 dataset object (FIXED: removed `meta$` from design)
deseq <- DESeqDataSetFromMatrix(countData = count_data, 
                                colData = meta,
                                design = ~ Bacteria)

# Run DESeq2 analysis
deseq <- DESeq(deseq)

# Extract results
DESeq_results <- results(deseq)

# Convert results to a DataFrame
log2pvalues_res = as.data.frame(DESeq_results)[, c("log2FoldChange", "pvalue")]

# Load ggplot2 for plotting
library(ggplot2)

# Remove NA values from p-values to avoid errors in log transformation
log2pvalues_res <- na.omit(log2pvalues_res)

# Add a column for the -log10 of the p-value
log2pvalues_res$negLog10pvalue <- -log10(log2pvalues_res$pvalue)

# Thresholds for significance
log2fc_threshold <- 1  
pvalue_threshold <- 0.05  

# Create a new column for significance based on thresholds
log2pvalues_res$significance <- "Not Significant"
log2pvalues_res$significance[log2pvalues_res$pvalue < pvalue_threshold & abs(log2pvalues_res$log2FoldChange) > log2fc_threshold] <- "Significant"


pdf("compare_met_volcano.pdf", width=8, height=6)

# Create the volcano plot
print(
  ggplot(log2pvalues_res, aes(x = log2FoldChange, y = negLog10pvalue, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +  # Scatter plot with some transparency
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "blue")) + 
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  theme(legend.title = element_blank())
)
dev.off()

result_deseq = as.data.frame(DESeq_results)[,c("baseMean", "log2FoldChange", "padj")]
write.csv(result_deseq, "result_deseq.csv")

sorted_results <- result_deseq[order(result_deseq$padj, na.last=TRUE), ]
write.csv(sorted_results, "compare_met_results.csv")
