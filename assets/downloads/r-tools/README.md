# Kang's Bioinformatics R tools

This directory contains the downloadable files linked from
`bioinformatics-analysis.html`.

| File | Type | Main function |
|---|---|---|
| `VolcanoKang_0.1.1.tar.gz` | Recommended R source package | `generate_volcano_plot()` |
| `VolcanoKang_0.1.0.tar.gz` | Original R source package | `generate_volcano_plot()` |
| `Volcano_plot_function.R` | Standalone R script | `create_volcano_plot()` |
| `GO_KEGG_eg.R` | Bioconductor analysis script | `run_go_kegg()` |
| `GO_KEGG_Image.R` | Visualization script | `plot_enrichment_dot()` |

See the Bioinformatics Analysis page for installation requirements, input
columns, and example code.

The current VolcanoKang package source is also stored at
`packages/VolcanoKang`, allowing direct installation from GitHub:

```r
remotes::install_github(
  "kangjunho/Kangsin",
  subdir = "packages/VolcanoKang"
)
```
