###############################################################
# Description:
#   - Reads count matrix (gene x sample): 
#   - Performs differential expression analysis
#   - Outputs results
###############################################################

# --- 1. Load libraries ---
library(tidyverse)
library(DESeq2)
library(ggplot2)

# --- 2. Define working directory and input ---
setwd("./results/")  # Relative path (portable)
getwd()
# Input count matrix file (rows = genes, columns = samples) is generated using the prepDE.py script, which converts the output files from StringTie into a gene-level count matrix suitable for differential expression analysis.
# e.g., count_matrix.csv
countData <- as.matrix(read.csv("count_matrix.csv", row.names = "gene_id"))
colData <- read.table("Tran_sample.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
stopifnot(all(rownames(colData) %in% colnames(countData)))
countData <- countData[, rownames(colData)]

# --- 2. Filter lowly expressed genes ---
countData <- countData[rowMeans(countData) > 1, ]

# --- 3. Create DESeq2 object ---
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = unique(colData$condition))

# --- 4. Run DESeq2 ---
dds <- DESeq(dds)
conditions <- levels(dds$condition)
comparisons <- combn(conditions, 2, simplify = FALSE)
dir.create("degresults", showWarnings = FALSE)
for (pair in comparisons) {
  cond1 <- pair[1]  # baseline
  cond2 <- pair[2]  # treatment
  contrast_name <- paste0(cond2, "_vs_", cond1)
  res <- results(dds, contrast = c("condition", cond2, cond1))
  res_df <- as.data.frame(res)
  res_df <- res_df %>%
    dplyr::select(log2FoldChange, pvalue, padj) %>%
    rename(Log2FC = log2FoldChange, Pvalue = pvalue, FDR = padj)

out=cbind(res$log2FoldChange,res$pvalue,res$padj) 
rownames(out)=row.names(res) 
colnames(out)=c('log2FoldChange','Pvalue','FDR') 

# Filter DEGs
FCcut <- 2
FDRcut <- 0.05
diff=out[(!is.na(res$padj) & res$padj<FDRcut) & abs(res$log2FoldChange)>abs(log2(FCcut)),] 

out_deg <- paste0("degresults/DESeq2_DEG_", contrast_name, ".csv")
write.csv(diff, out_deg, quote = FALSE)
out_list <- paste0("degresults/degresults_DEGlist_", contrast_name, ".txt")
write.table(rownames(diff),
              paste0("results/DESeq2_DEGlist_", contrast_name, ".txt"),

              quote = FALSE, row.names = FALSE, col.names = FALSE)

