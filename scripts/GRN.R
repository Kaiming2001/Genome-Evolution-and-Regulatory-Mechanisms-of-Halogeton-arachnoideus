# ===============================================================
# Build transcriptional regulatory network using GENIE3
# ===============================================================

library(readxl)
library(dplyr)
library(GENIE3)
library(igraph)
library(readxl)
# ---------- Project paths ----------
setwd("./result/") #The path where the TPM calculation results of the transcripts are stored

# ---------- Input files ----------
expr <- read_excel("200.TPM.xlsx", sheet = 1) # Split transcript TPM expression matrix into two subsets (200 mM and 400 mM NaCl treatments),save them separately as "200.TPM.xlsx" and "400.TPM.xlsx",and construct independent gene regulatory networks (GRNs) for each condition using GENIE3.
expr <- as.data.frame(expr)
rownames(expr) <- expr[[1]]
expr <- expr[, -1]
# ---------- TF list ----------
tf_list <- read_excel("200.TPM.xlsx", sheet = 2) # Sheet2 contains TFs identified based on iTAK annotation.
tf_list <- as.character(tf_list[[1]])
tf_list <- tf_list[!is.na(tf_list)]

weight_matrix <- GENIE3(as.matrix(expr), regulators = tf_list, nTrees = 1000)

link_list <- getLinkList(weight_matrix)

colnames(link_list) <- c("regulator", "target", "weight")

strong_links <- link_list %>% filter(weight > 0.05)

write.csv(strong_links, "200cy_edges_strong.csv")

all_genes_strong <- unique(c(strong_links$regulator, strong_links$target))

nodes_info_strong <- data.frame(
  gene = all_genes_strong,
  is_TF = ifelse(all_genes_strong %in% tf_list, "TF", "Target")
)
write.csv(nodes_info_strong, "200cy_nodes_strong.csv")
# ---------- Differential regulation network ---------- 
# ---------- Input files ----------
expr <- read_excel("200TF.TPM.xlsx")  #Only the values of TF expression levels were selected.
tf_list <- as.character(expr[[1]]) 
expr <- as.data.frame(expr)
rownames(expr) <- expr[[1]]           
expr <- expr[, -1]  
# ----- Extract sample groups -----
cols_6h <- grep("Hara-200mM-6", colnames(expr), value = TRUE)
cols_24h  <- grep("Hara-200mM-24", colnames(expr), value = TRUE)
expr_6h <- expr[, cols_6h]
expr_24h  <- expr[, cols_24h]

weight_6h <- GENIE3(as.matrix(expr_6h), regulators = tf_list, nTrees = 1000)
weight_24h <- GENIE3(as.matrix(expr_24h),  regulators = tf_list, nTrees = 1000)

link_6h <- getLinkList(weight_6h)
link_24h  <- getLinkList(weight_24h)

colnames(link_6h)[3] <- "weight_6h"
colnames(link_24h)[3]  <- "weight_24h"
colnames(link_6h) <- c("regulator", "target", "weight_6h")
colnames(link_24h)  <- c("regulator", "target", "weight_24h")

merged_links <- full_join(link_6h, link_24h, by = c("regulator", "target")) %>%
  mutate(
    weight_6h = ifelse(is.na(weight_6h), 0, weight_6h),
    weight_24h = ifelse(is.na(weight_24h), 0, weight_24h),
    diff = weight_24h - weight_6h
  )

up_regulated <- merged_links %>% filter(diff > 0.03)
down_regulated <- merged_links %>% filter(diff < -0.03)
# Up-regulated edges denote stronger regulatory interactions at 24 h, whereas 'down' indicates stronger interactions at 6 h.
write.csv(up_regulated, "200.24hregulated_TF_edges.csv", row.names = FALSE)
write.csv(down_regulated, "200.6hregulated_TF_edges.csv", row.names = FALSE)

