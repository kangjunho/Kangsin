# Step 12. Rank-based GSEA template
if (!requireNamespace("msigdbr", quietly=TRUE)) install.packages("msigdbr")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
if (!requireNamespace("fgsea", quietly=TRUE)) BiocManager::install("fgsea", ask=FALSE, update=FALSE)
library(msigdbr); library(fgsea)
gene_rank <- deg_result$t
names(gene_rank) <- deg_result$gene
gene_rank <- sort(gene_rank[!is.na(gene_rank) & !duplicated(names(gene_rank))], decreasing=TRUE)
msig_h <- msigdbr(species="Homo sapiens", collection="H")
pathways <- split(msig_h$gene_symbol, msig_h$gs_name)
gsea_result <- fgsea(pathways=pathways, stats=gene_rank, minSize=15, maxSize=500)
gsea_result <- gsea_result[order(gsea_result$padj), ]
gsea_result$leadingEdge <- vapply(gsea_result$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
dir.create("results", showWarnings=FALSE)
write.csv(gsea_result, "results/GSEA_Hallmark.csv", row.names=FALSE, fileEncoding="UTF-8-BOM")
