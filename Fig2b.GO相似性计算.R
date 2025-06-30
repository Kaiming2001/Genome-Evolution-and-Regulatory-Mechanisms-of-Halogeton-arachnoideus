##获取GOid之间的相似性
BiocManager::install("simplifyEnrichment")
BiocManager::install("org.Hs.eg.db")
library(simplifyEnrichment)
library(readxl)
df <- read_excel("D:/盐生草/结果/tree/基因家族扩张和收缩/GO/GOID.xlsx")
GOIDs <- as.character(df$GOID)
head(GOIDs)
mat = GO_similarity(GOIDs)
df1=simplifyGO(mat,method = "kmeans")
dist_mat <- as.dist(1 - mat)  # 相似性转为距离
hc <- hclust(dist_mat, method = "average")  # 层次聚类
groups <- cutree(hc, k = 5)
#转换为cytoscape输入数据
net = reshape2::melt(mat)
net = net[which(net$value > 0.3),]  #删掉一些相似性特别低的点
net = net[which(as.character(net$Var1)>as.character(net$Var2)),]   #删除回路，比如有A到B线和B到A的线，只保留A到B的
head(net)
head(df1)
setwd("D:/盐生草/结果/盐爪爪盐节木/GO/")
write.table(net,file = "net.txt",col.names = T,row.names = F,quote = F,sep = "\t")
write.table(df1,file = "df.txt",col.names = T,row.names = F,quote = F,sep = "\t")