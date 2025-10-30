# GO Enrichment Analysis with topGO
# --- 1. Load required libraries ---
library(Rgraphviz)
library(topGO)
library(dplyr)
library(data.table)
library(readxl)
# --- 2. Set working directory ---
# Use a temporary path that contains the background GO annotation (result/restored_GO.pl.out) results and the gene list for comparison (e.g., expanded gene sets).
setwd("./result/")
getwd()
# --- 3. Import gene-to-GO annotation table ---
# Expected format: two columns (gene_id, go_id)
# Example file: "result/restored_GO.pl.out"
d <- read.table("restored_GO.pl.out", header = FALSE, sep = "\t", fill = TRUE, quote = "")
colnames(d) <- c("gene_id", "go_id")
# Convert to list format required by topGO
all_go <- lapply(split(d, d$gene_id), function(x) unique(x$go_id))
geneNames <- names(all_go)
# --- 4. Import target gene list (e.g., candidate genes or duplicated genes) ---
# Example file: "Hara.ex_all.id" (one gene ID per line)
gene <- read.table("Hara.ex_all.id", header = FALSE, sep = "\t", fill = TRUE, quote = "")
colnames(gene) <- c("geneid")
outlier_gene <- as.vector(gene$geneid)
# --- 5. Create factor indicating whether each gene is in the target set ---
geneList <- factor(as.integer(geneNames %in% outlier_gene))
names(geneList) <- geneNames
# --- 6. Build topGO data object ---
GOdata <- new("topGOdata",
              ontology = "BP",                # Biological Process ontology
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = all_go)

# --- 7. Run enrichment test (Fisher's exact test) ---
res.classic <- runTest(GOdata,algorithm="classic",statistic="fisher")
res.elim  <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
# --- 8. Extract all GO terms and generate summary table ---
allGO <- usedGO(GOdata)
tab_all <- GenTable(GOdata,
                    classicFisher = res.classic,
                    elimFisher    = res.elim,
                    orderBy       = "classicFisher",
                    ranksOf       = "classicFisher",
                    topNodes      = length(allGO),
                    numChar       = 1000)

# --- 9. Multiple testing correction (Benjamini-Hochberg FDR) ---
tab_all$classicFisher <- as.numeric(tab_all$classicFisher)
tab_all$classicFDR <- p.adjust(tab_all$classicFisher, method = "BH")
tab_all$elimFisher <- as.numeric(tab_all$elimFisher)
tab_all$elimFDR <- p.adjust(tab_all$elimFisher, method = "BH")
# --- 10. Save results ---
write.csv(tab_all, file = "GO_enrichment_results.csv", row.names = FALSE)
