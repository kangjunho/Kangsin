# Step 08. limma DEG template
library(limma)
# Required objects: gene_matrix, sample_info with a verified group column
group <- factor(sample_info$group, levels = c("Control", "Disease"))
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
contrast_matrix <- makeContrasts(Disease_vs_Control = Disease - Control, levels = design)
fit <- lmFit(gene_matrix, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)
deg_result <- topTable(fit, coef = "Disease_vs_Control", number = Inf, adjust.method = "BH")
deg_result$status <- "Not significant"
deg_result$status[deg_result$adj.P.Val < .05 & deg_result$logFC >= 1] <- "Up"
deg_result$status[deg_result$adj.P.Val < .05 & deg_result$logFC <= -1] <- "Down"
