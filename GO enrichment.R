#GO
BiocManager::install("topGO")
BiocManager::install('Rgraphviz')
BiocManager::install('data.table')
library(Rgraphviz)
library(topGO)
library(dplyr)
library(data.table)
library(readxl)
getwd()
setwd("D:/Adata/trans/seeding/root400")
###inport GOid
d <- read.table("D:/hara/result/LTR_overlapgene/restored_GO.pl.out", header=F, sep="\t", fill=TRUE,quote="")
colnames(d)=c("gene_id","go_id")
d$gene_id=as.character(d$gene_id)
d$go_id=as.character(d$go_id)
all_go <- lapply(split(d, sub("//.//d+$", "", d[, 1])), 
                 function(x) unique(x[, 2]))
geneNames=names(all_go)

###inport geneid
gene=read.table("Hara.rapid.ex_all.id",header=F)

head(gene)
colnames(gene)=c("geneid")
outlier_gene = as.vector(gene$geneid)
head(outlier_gene)
geneList=factor(as.integer(geneNames %in% outlier_gene))
names(geneList)=geneNames
head(geneList) 
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,annot=annFUN.gene2GO, gene2GO = all_go)  
restRes=runTest(GOdata,algorithm="classic",statistic="fisher")
restRes
genetable <- GenTable(GOdata, p.value = restRes, orderBy = "p.value")
###choose GO term with pvalue <0.01
gene_table=GenTable(GOdata,Fisher.p=restRes,topNodes=82,numChar=5000000)
###
gene_table_total=GenTable(GOdata,Fisher.p=restRes,topNodes=269)
allGO=usedGO(GOdata)
pvalues=gene_table_total$Fisher.p
gene_table$adjust.p=round(head(p.adjust(pvalues,method="BH"),82),5000)   
write.table(gene_table,file="Hara.rapid.ex.GO.csv",sep=",", quote=F, row.names=F, col.names=T)
