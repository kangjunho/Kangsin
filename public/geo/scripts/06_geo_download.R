# Step 06. GEOquery download
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
library(GEOquery)
gse_id <- "GSE_ACCESSION" # replace before running
gset <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
eset <- gset[[1]]
expr_matrix <- exprs(eset)
sample_info <- pData(eset)
feature_info <- fData(eset)
dir.create("data_processed", showWarnings = FALSE)
saveRDS(eset, "data_processed/expression_set.rds")
write.csv(sample_info, "data_processed/sample_metadata.csv", row.names = FALSE, fileEncoding = "UTF-8")
