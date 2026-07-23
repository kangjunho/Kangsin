# Step 07. Preprocessing template
eset <- readRDS("data_processed/expression_set.rds")
expr_matrix <- Biobase::exprs(eset)
sample_info <- Biobase::pData(eset)
feature_info <- Biobase::fData(eset)
qx <- quantile(expr_matrix, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
log2_needed <- qx[5] > 100 || (qx[6] - qx[1] > 50)
if (log2_needed) { expr_matrix[expr_matrix <= 0] <- NA; expr_matrix <- log2(expr_matrix) }
boxplot(expr_matrix, outline = FALSE, las = 2)
# Set the correct symbol column after inspecting colnames(feature_info).
stop("Inspect feature_info and define probe_map before continuing.")
