# Step 09. Visualization template
library(ggplot2); library(pheatmap)
dir.create("figures", showWarnings = FALSE); dir.create("results", showWarnings = FALSE)
deg_result$gene <- rownames(deg_result)
p_volcano <- ggplot(deg_result, aes(logFC, -log10(adj.P.Val), color = status)) +
  geom_point(alpha = .65) +
  scale_color_manual(values = c(Down="#2166AC", `Not significant`="grey75", Up="#B2182B")) +
  theme_classic()
ggsave("figures/volcano.pdf", p_volcano, width=7, height=6)
write.csv(deg_result, "results/DEG_Disease_vs_Control.csv", row.names=FALSE, fileEncoding="UTF-8-BOM")
