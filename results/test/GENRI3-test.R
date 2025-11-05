library(readxl)
library(GENIE3)
library(dplyr)
library(tidyr)
library(parallel)

# ---------- 1. Input data ----------
setwd(" ")

expr <- read_excel("200.TPM.xlsx", sheet = 1)
expr <- as.data.frame(expr)
rownames(expr) <- as.character(expr[[1]])
expr <- expr[, -1]

tf_list <- read_excel("200.TPM.xlsx", sheet = 2)
tf_list <- as.character(tf_list[[1]])
tf_list <- tf_list[!is.na(tf_list)]
expr <- as.matrix(expr)


# ---------- 2. Stable selection function ----------
run_genie3_bootstrap <- function(expr_mat, tf_list, n_boot = 100, subsample = 0.8, nCores = 4) {
  genes <- rownames(expr_mat)
  edge_list_all <- list()
  
  cl <- makeCluster(nCores)
  clusterEvalQ(cl, library(GENIE3))
  clusterExport(cl, c("expr_mat", "tf_list", "subsample", "genes"), envir = environment())
  
  edge_list_all <- parLapply(cl, 1:n_boot, function(i) {
    set.seed(i)
    idx <- sample(1:ncol(expr_mat), size = ceiling(ncol(expr_mat) * subsample), replace = FALSE)
    expr_sub <- expr_mat[, idx, drop = FALSE]
    wMat <- GENIE3(expr_sub, regulators = tf_list, nTrees = 1000, verbose = FALSE)
    links <- getLinkList(wMat)
    links
  })
  stopCluster(cl)
  return(edge_list_all)
}

# ---------- 3. Choose stable operation ----------
edge_list_all <- run_genie3_bootstrap(expr, tf_list, n_boot = 100, subsample = 0.8, nCores = 6)

# ---------- 4. Calculate the frequency of the appearance of the edges ----------
all_edges <- bind_rows(edge_list_all) %>%
  group_by(regulatoryGene, targetGene) %>%
  summarise(freq = n()/100, mean_weight = mean(weight, na.rm = TRUE), .groups = "drop")

# ---------- 5. Select high-confidence edges ----------
stable_edges <- all_edges %>% filter(freq >= 0.7) %>% arrange(desc(freq), desc(mean_weight))

# ---------- 6. Export ----------
write.csv(stable_edges, "1GRN_stable_edges_200mM.csv", row.names = FALSE)

# ---------- 7. Summary ----------
summary(stable_edges$freq)
cat("最终稳定边数：", nrow(stable_edges), "\n")

# ---------- Using the entire sample ----------
full_weight <- GENIE3(expr, regulators=tf_list, nTrees=1000)
full_links <- getLinkList(full_weight)
write.csv(full_links, "GRN_edges_200mM.csv", row.names = FALSE)
# ---------- Calculate the correlation ----------
stable_edges <- read.csv("GRN_stable_edges_400mM.csv")
full_links <- read.csv("GRN_edges_400mM.csv")
merged_df <- merge(full_links, stable_edges, by=c("regulatoryGene","targetGene"))
names(merged_df)

merge(full_links, stable_edges, by=c("regulatoryGene","targetGene")) %>%
  summarise(cor = cor(merged_df$weight, merged_df$mean_weight, use = "complete.obs")) #默认method = "pearson"