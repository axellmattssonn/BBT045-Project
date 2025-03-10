#Load library
library(DESeq2)

#load datafiles
count_data = read.csv("counts.filtered.csv", row.names = 1)
metadata = read.csv("patient_list.csv", row.names = 1)

#trim spaces etc in names
rownames(metadata) = trimws(rownames(metadata))

#match the rows and column in each
count_data = count_data[, match(rownames(metadata), colnames(count_data))]
metadata$HealthCondition <- as.factor(metadata$HealthCondition)

#Create DESeq2 dataset object
dds = DESeqDataSetFromMatrix(countData = count_data,
                          colData = metadata,
                         design = ~ HealthCondition)

dds = DESeq(dds)

deseqresults = results(dds)



# Extract p-value and log2FoldChange from DESeq2 results
log2pvalues_res = as.data.frame(deseqresults)[, c("log2FoldChange", "pvalue")]

# Load ggplot2 for plotting
library(ggplot2)

# Add a column for the -log10 of the p-value (for plotting on the y-axis)
log2pvalues_res$negLog10pvalue <- -log10(log2pvalues_res$pvalue)

# Thresholds for significance
log2fc_threshold <- 1  
pvalue_threshold <- 0.05  

# Create a new column for significance based on thresholds
log2pvalues_res$significance <- "Not Significant"
log2pvalues_res$significance[log2pvalues_res$pvalue < pvalue_threshold & abs(log2pvalues_res$log2FoldChange) > log2fc_threshold] <- "Significant"


pdf("/cephyr/NOBACKUP/groups/bbt045_2025/groups/group_gaswar/volcano_plot.pdf", width=8, height=6)
print(
  ggplot(log2pvalues_res, aes(x = log2FoldChange, y = negLog10pvalue, color = significance)) +
    geom_point(alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    labs(title = "Volcano Plot",
         x = "Log2 Fold Change",
         y = "-Log10(p-value)") +
    theme_minimal() +
    theme(legend.title = element_blank())
)
dev.off()

result_deseq = as.data.frame(deseqresults)[,c("baseMean", "log2FoldChange", "padj")]
write.csv(result_deseq, "result_deseq.csv")

sorted_results <- result_deseq[order(result_deseq$padj, na.last=TRUE), ]
write.csv(sorted_results, "sorted_deseq_results.csv")



