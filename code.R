# Load the necessary packages
library("dplyr")
library("ggplot2")

# Read the data
url <- "https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt"
RNAseq_data <- read.table("https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt", sep = " ", header = TRUE)
head(RNAseq_data)
View(RNAseq_data)

RNAseq_data <- RNAseq_data %>%
  mutate(logpvalue = -log10(pvalue))


#	TASK 1: Generate a volcano plot
ggplot(RNAseq_data, aes(x = log2FoldChange, y = logpvalue)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10(p-value)", title = "Volcano Plot") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +  # Example significance threshold
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") 


#	TASK 2: Determine the upregulated genes (Genes with Log2FC > 1 and p-value < 0.01)

# Filter for upregulated genes
upregulated_genes <- RNAseq_data %>%
  filter(log2FoldChange > 1, pvalue < 0.01)

# View the results
print(upregulated_genes)

# Write the results to a CSV file
write.csv(upregulated_genes, "upregulated_genes.csv", row.names = FALSE)


#	TASK 3: Determine the downregulated genes (Genes with Log2FC < -1 and p-value < 0.01)

# Filter for downregulated genes
downregulated_genes <- RNAseq_data %>%
  filter(log2FoldChange < -1, pvalue < 0.01)

# View the results
print(downregulated_genes)

# Write the results to a CSV file
write.csv(downregulated_genes, "downregulated_genes.csv", row.names = FALSE)


# Visualize the significantly expressed genes

# Define thresholds
upregulated_threshold <- 1
downregulated_threshold <- -1
pvalue_threshold <- 0.01

# Filter for upregulated and downregulated genes
RNAseq_data$significance <- "Not Significant"
RNAseq_data$significance[RNAseq_data$log2FoldChange > upregulated_threshold & RNAseq_data$pvalue < pvalue_threshold] <- "Upregulated"
RNAseq_data$significance[RNAseq_data$log2FoldChange < downregulated_threshold & RNAseq_data$pvalue < pvalue_threshold] <- "Downregulated"

# Create the volcano plot
ggplot(RNAseq_data, aes(x = log2FoldChange, y = -log10(pvalue), color = significance)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "Volcano Plot of Significant Genes",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black") + # Significance line
  geom_vline(xintercept = c(upregulated_threshold, downregulated_threshold), linetype = "dashed", color = "black") + # Fold change lines
  theme(legend.position = "top")


write.csv(RNAseq_data, "RNAseq_data.csv", row.names = FALSE)

#-------------------------------------------------------------

# Filter for significant genes
significant_genes <- RNAseq_data %>%
  filter(significance != "Not Significant")

# Count the number of upregulated and downregulated genes
gene_counts <- significant_genes %>%
  group_by(significance) %>%
  summarise(count = n())

# Create the bar plot
ggplot(gene_counts, aes(x = significance, y = count, fill = significance)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Upregulated and Downregulated Genes",
       x = "Gene Regulation",
       y = "Count") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal()
