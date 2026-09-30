# VolcanoKang

An R package for generating volcano plots from differential expression results.

## Install from GitHub

```r
install.packages("remotes")
remotes::install_github(
  "kangjunho/Kangsin",
  subdir = "packages/VolcanoKang",
  upgrade = "never"
)
```

## Use

The input data frame must contain `gene`, `log2FoldChange`, and `pvalue`.

```r
library(VolcanoKang)

p <- generate_volcano_plot(
  deg,
  log2fc_cutoff = 1,
  pvalue_cutoff = 0.05,
  label_genes = TRUE,
  label_top_n = 20
)

p
```
