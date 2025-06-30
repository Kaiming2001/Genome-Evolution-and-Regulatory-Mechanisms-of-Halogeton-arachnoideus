##Obtain the similarity between GO ids
BiocManager::install("simplifyEnrichment")
BiocManager::install("org.Hs.eg.db")
library(simplifyEnrichment)
library(readxl)
df <- read_excel("D:/hara/result/tree/genefamily/GO/GOID.xlsx")
GOIDs <- as.character(df$GOID)
head(GOIDs)
mat = GO_similarity(GOIDs)
df1=simplifyGO(mat,method = "kmeans")
dist_mat <- as.dist(1 - mat)  
hc <- hclust(dist_mat, method = "average")  
groups <- cutree(hc, k = 5)
#Convert to CytoScape input data
net = reshape2::melt(mat)
net = net[which(net$value > 0.3),] 
net = net[which(as.character(net$Var1)>as.character(net$Var2)),] 
head(net)
head(df1)
setwd("D:/hara/result/tree/genefamily/GO/")
write.table(net,file = "net.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(df1,file = "df.txt",col.names = T,row.names = F,quote = F,sep = "\t")