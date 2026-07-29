# Step 04. GEO2R 코드를 RStudio에서 실행하기
# Dataset: GSE69657

library(GEOquery)
library(limma)
library(umap)

# 다운로드 파일을 임시 폴더가 아닌 프로젝트에 보관합니다.
options(timeout = 600)
options(download.file.method = "libcurl")
dir.create("data_raw", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

gset <- getGEO(
  "GSE69657",
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  destdir = "data_raw"
)

# 여러 platform이 있으면 분석할 GPL을 선택합니다.
if (length(gset) > 1) {
  idx <- grep("GPL570", names(gset))
  if (length(idx) != 1) {
    stop("GPL570을 하나만 선택하지 못했습니다. names(gset)을 확인하세요.")
  }
} else {
  idx <- 1
}
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))

# GEO2R에서 생성된 그룹 정보입니다.
gsms <- "111011100000011110011001011011"
sml <- strsplit(gsms, split = "")[[1]]
if (length(sml) != ncol(gset)) {
  stop("그룹 코드 수와 샘플 수가 다릅니다. GEO2R의 gsms 값을 확인하세요.")
}

# 필요한 경우 log2 변환합니다.
ex <- exprs(gset)
qx <- as.numeric(
  quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
)
LogC <- (qx[5] > 100) ||
  (qx[6] - qx[1] > 50 && qx[2] > 0)

if (LogC) {
  ex[ex <= 0] <- NaN
  exprs(gset) <- log2(ex)
}

# 그룹과 design matrix를 만듭니다.
gs <- factor(sml)
groups <- make.names(c("Responder", "Non_responder"))
levels(gs) <- groups
gset$group <- gs

design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

# 반드시 샘플별 그룹 배정을 확인합니다.
group_check <- data.frame(
  sample = sampleNames(gset),
  code = sml,
  group = as.character(gs)
)
print(group_check)
print(table(group_check$group))

gset <- gset[complete.cases(exprs(gset)), ]

# limma 차등발현 분석
fit <- lmFit(gset, design)
cts <- paste(groups[1], groups[2], sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

# number = Inf로 전체 probe 결과를 반환합니다.
tT <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf
)

# GSE69657 실습에 필요한 열만 선택합니다.
print(colnames(tT))
needed_cols <- c(
  "ID", "adj.P.Val", "P.Value", "t", "B", "logFC", "Gene.symbol"
)
missing_cols <- setdiff(needed_cols, colnames(tT))
if (length(missing_cols) > 0) {
  stop(
    "다음 열이 없습니다: ",
    paste(missing_cols, collapse = ", ")
  )
}
tT <- tT[, needed_cols, drop = FALSE]

write.csv(
  tT,
  "results/GSE69657_limma_all_results.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 실습용 기준: nominal P-value < 0.05, |logFC| > 1
up_degs <- tT[
  !is.na(tT$P.Value) &
    !is.na(tT$logFC) &
    tT$P.Value < 0.05 &
    tT$logFC > 1,
]

down_degs <- tT[
  !is.na(tT$P.Value) &
    !is.na(tT$logFC) &
    tT$P.Value < 0.05 &
    tT$logFC < -1,
]

print(nrow(up_degs))
print(nrow(down_degs))

write.csv(
  up_degs,
  "results/up_degs.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  down_degs,
  "results/down_degs.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# GEO2R 기본 결과 그림
hist(
  tT$adj.P.Val,
  col = "grey",
  border = "white",
  xlab = "P-adj",
  ylab = "Number of probes",
  main = "P-adj value distribution"
)

dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05, lfc = 0)
vennDiagram(dT, circle.col = palette())

t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good], main = "Moderated t statistic")

ct <- 1
volcanoplot(
  fit2,
  coef = ct,
  main = colnames(fit2)[ct],
  pch = 20,
  highlight = length(which(dT[, ct] != 0)),
  names = rep("+", nrow(fit2))
)

plotMD(
  fit2,
  column = ct,
  status = dT[, ct],
  legend = FALSE,
  pch = 20,
  cex = 1
)
abline(h = 0)

# 발현분포와 UMAP
ex <- exprs(gset)
ord <- order(gs)
palette(
  c(
    "#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
    "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"
  )
)

par(mar = c(7, 4, 2, 1))
plot_title <- paste("GSE69657", annotation(gset), sep = " / ")
boxplot(
  ex[, ord],
  boxwex = 0.6,
  notch = TRUE,
  main = plot_title,
  outline = FALSE,
  las = 2,
  col = gs[ord]
)
legend("topleft", groups, fill = palette(), bty = "n")

par(mar = c(4, 4, 2, 1))
density_title <- paste(
  "GSE69657",
  annotation(gset),
  "value distribution",
  sep = " / "
)
plotDensities(ex, group = gs, main = density_title, legend = "topright")

ex_umap <- na.omit(ex)
ex_umap <- ex_umap[!duplicated(ex_umap), ]
ump <- umap(t(ex_umap), n_neighbors = 13, random_state = 123)

par(mar = c(3, 3, 2, 6), xpd = TRUE)
plot(
  ump$layout,
  main = "UMAP plot, nbrs=13",
  xlab = "",
  ylab = "",
  col = gs,
  pch = 20,
  cex = 1.5
)
legend(
  "topright",
  inset = c(-0.15, 0),
  legend = levels(gs),
  pch = 20,
  col = 1:nlevels(gs),
  title = "Group",
  pt.cex = 1.5
)

# maptools와 pointLabel()은 사용하지 않습니다.
plotSA(fit2, main = "Mean variance trend, GSE69657")
