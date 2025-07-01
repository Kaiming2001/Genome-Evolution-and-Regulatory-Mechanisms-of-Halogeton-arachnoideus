library(tidyverse)
library(DESeq2)
library(ggplot2)
countData <- as.matrix(read.csv("D:/Adata/trans/DEG/2vs4.csv",row.names="gene_id"))
countData <- countData[rowMeans(countData)>1,]
head(countData)
condition <- factor(c(rep("Hara_200mM_30h",3),rep("Hara_400mM_30h",3)))

head(condition)
sampleNames <- colnames(countData)
colData <- data.frame(sampleName = sampleNames, condition = condition)
rownames(colData) <- colData$sampleName
colData$condition = factor(colData$condition)
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ condition)
dds <- DESeq(dds) 
res <- results(dds) 
res
out=cbind(res$log2FoldChange,res$pvalue,res$padj) 
rownames(out)=row.names(res) 
colnames(out)=c('log2FoldChange','Pvalue','FDR') 
#filter condition
FCcut=2 
FDRcut=0.05 
diff=out[(!is.na(res$padj) & res$padj<FDRcut) & abs(res$log2FoldChange)>abs(log2(FCcut)),] 
write.table(out, "D:/Adata/trans/DEG/DESeq2_AllResult.2vs4.csv", sep = ",", quote = F, row.names = T, col.names = T)
write.table(diff, "D:/Adata/trans/DEG/DESeq2_DEG.2vs4.csv", sep = ",", row.names = TRUE, quote = FALSE)
write.table(rownames(diff), "D:/Adata/trans/DEG/DESeq2_DEGlist.2vs4.csv", sep = ",", quote = F, row.names = F, col.names = F)


#volcano plot
res_all <- as.data.frame(res) %>% filter(padj != 'NA')
res_all$type <- ifelse(res_all$padj < 0.05,
                       ifelse(abs(res_all$log2FoldChange) > 1,
                              ifelse(res_all$log2FoldChange < -1,'down','up'),'noSig'),'noSig')
# choose "up" and "down" gene
up_genes <- res_all %>% filter(type == "up")
down_genes <- res_all %>% filter(type == "down")
# save "up" and "down" genes as .CSV file
write.table(up_genes, "DESeq2_up_genes.2vs4.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(down_genes, "DESeq2_down_genes.2vs4.csv", sep = ",", quote = FALSE, row.names = TRUE, col.names = TRUE)
head(res_all,5)
table(res_all$type)
ggplot(res_all,aes(x = log2FoldChange,y = -log10(pvalue))) +
  geom_point(aes(color = type),alpha = 0.5,size = 2) +
  theme_bw(base_size = 16) +
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black')) +
  scale_color_manual(name = '',
                     # color or three types
                     values = c('up'='#DA1212','noSig'='grey','down'='blue'),
                     # legend labels
                     label = c('down'='down (45)','noSig'='noSig (3859)','up'='up (36)')) +
  geom_hline(yintercept = -log10(0.05),lty = 'dashed',linewidth = 0.8) +
  geom_vline(xintercept = c(-1,1),lty = 'dashed',linewidth = 0.8) +
  scale_x_continuous(breaks = c(-6,-4,-2,-1,0,1,2,4,6)) +
  ggtitle(' ')
