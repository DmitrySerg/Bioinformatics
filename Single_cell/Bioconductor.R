library(scater)
library(SC3)

#source("https://bioconductor.org/biocLite.R")
#biocLite("SC3")

# cell annotation
ann <- data.frame(cell_type1 = colnames(treutlein))
pd <- new("AnnotatedDataFrame", data = ann)
# cell expression

tmp <- pbmc.data



tmp <- treutlein
colnames(tmp) <- rownames(ann)
# SCESEt object
sceset <- newSCESet(fpkmData = tmp, phenoData = pd, logExprsOffset = 1)

sceset <- newSCESet(fpkmData = tmp, logExprsOffset = 1)

sceset <- calculateQCMetrics(sceset)
plotPCA(sceset)

plotPCA(sceset, colour_by = "cell_type1")


# Note that n_cores = 1 is required for compilation of this vignette.
# Please remove this parameter when running on your computer:
# sceset <- sc3(sceset, ks = 2:4, biology = TRUE)
sceset <- sc3(sceset, ks = 2:7, biology = TRUE)

sc3_interactive(sceset)
sc3_export_results_xls(sceset)


p_data <- pData(sceset)
head(p_data[ , grep("sc3_", colnames(p_data))])

plotPCA(
  sceset, 
  colour_by = "sc3_3_clusters", 
  size_by = "sc3_3_log2_outlier_score"
)


f_data <- fData(sceset)
head(f_data[ , grep("sc3_", colnames(f_data))])


plotFeatureData(
  sceset, 
  aes(
    x = sc3_3_markers_clusts, 
    y = sc3_3_markers_auroc, 
    colour = sc3_3_markers_padj
  )
)


sc3_plot_consensus(sceset, k = 3)

sc3_plot_consensus(
  sceset, k = 3, 
  show_pdata = c(
    "cell_type1", 
    "log10_total_features",
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)

sc3_plot_silhouette(sceset, k = 3)

sc3_plot_expression(sceset, k = 3)


sc3_plot_expression(
  sceset, k = 3, 
  show_pdata = c(
    "cell_type1", 
    "log10_total_features",
    "sc3_3_clusters", 
    "sc3_3_log2_outlier_score"
  )
)


sceset <- sc3_prepare(sceset, ks = 2:4, n_cores = 1)
str(sceset@sc3)

sceset <- sc3_estimate_k(sceset)
sceset <- sc3_calc_dists(sceset)
sceset <- sc3_calc_transfs(sceset)
sceset <- sc3_kmeans(sceset)

names(sceset@sc3$kmeans)