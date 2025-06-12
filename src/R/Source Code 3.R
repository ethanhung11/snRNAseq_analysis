##################################################################################################################
# Title: Source Code File 3
# Content: Pseudobulk differential analysis within adipocyte subclusters - continued from the Source Code File 2
##################################################################################################################

# Load all required libraries/packages
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(Matrix)
library(reshape2)
library(S4Vectors)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(limma)
library(magrittr)
library(gridExtra)
library(knitr)
library(limma)



#########################################
##  COUNT DATA EXTRACTION & FILTERING  ##
#########################################

# Extract raw counts and metadata 
counts <- GetAssayData(object = Adipocyte.clustered, slot = "counts", assay = "RNA")
metadata <- Adipocyte.clustered@meta.data  # %>% select(orig.ident, depot, diet)

# Set up metadata as desired for aggregation and DE analysis
Idents(object = Adipocyte.clustered) <- "integrated_snn_res.0.3"
table(Idents(object = Adipocyte.clustered))
metadata$cluster_id <- factor(Adipocyte.clustered@active.ident) 
metadata$orig.ident <- gsub("_", ".", metadata$orig.ident) # for easing the subsetting step 

# Create a single cell experiment object 
sce <- SingleCellExperiment(assay = list(counts = counts), colData = metadata)
dim(colData(sce)) # (number of cells) * (number of meta columns)
# 10750    24

# Additional QC filtering
## 1) Remove undetected/lowly expressed genes
sce <- sce[rowSums(counts(sce) > 1) >= 10, ]
dim(sce)  # 12766 10750
## 2) Remove cells with few or many detected genes 
qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = T)
sce <- sce[, !ol]
dim(sce)  # 12766 9572



##############################################
##  COUNT DATA PREPARATION FOR AGGREGATION  ##
##############################################

# Named vector of cluster names
## 0   1   2   3   4   5 
## 0 =  Ad1, 1 = Ad3, 2 = Ad2, 3 = Ad6, 4 = Ad5, 5 = Ad4
kids <- purrr::set_names(levels(sce$cluster_id))

# Total number of clusters
## 6
nk <- length(kids)

# Named vector of sample names
sids <- purrr::set_names(levels(as.factor(sce$orig.ident)))

# Total number of samples 
## 4
ns <- length(sids)

# Generate sample-level metadata
## Determine the number of cells per sample
table(sce$orig.ident)
## Turn class "table" into a named vector of cells per sample
n_cells <- table(sce$orig.ident) %>%  as.vector()
names(n_cells) <- names(table(sce$orig.ident))

## Match the named vector with metadata to combine it
m <- match(names(n_cells), sce$orig.ident)

## Create the sample level metadata by selecting specific columns
ei <- data.frame(colData(sce)[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("orig.ident", "depot", "diet", "n_cells")
kable(ei)



#########################
##  COUNT AGGREGATION  ##  
#########################

# Aggregate the counts per sample_id and cluster_id
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cluster_id", "orig.ident")]
groups$orig.ident <- factor(groups$orig.ident)

# Aggregate across cluster-sample groups
# Each row corresponds to aggregate counts for a cluster-sample combo
pb <- aggregate.Matrix(t(counts(sce)), 
                       groupings = groups, fun = "sum") 

# class(pb): a sparse Matrix of class "dgCMatrix": sample (row) x gene (col)
# str(pb)
# dim(pb): 6 clusters * 4 samples = 24 rows; 12766 genes as columns
# 24 12766
pb[1:8, 1:8]



########################
##  SPLIT/SUBSETTING  ##  
########################

# Not every cluster is present in all samples; create a vector that represents how to split samples
splitf <- sapply(stringr::str_split(rownames(pb), 
                                    pattern = "_",n = 2), `[`, 1)
splitf <- factor(splitf, levels = c(0,1,2,3,4,5))

# Split data and turn into a list
# Each component corresponds to a cluster; storing associated expression matrix (counts)
# Transform data i.e, so rows are genes and columns are samples 
pb_s <- split.data.frame(pb,splitf) %>%
  lapply(function(u) 
    set_colnames(t(u), gsub(".*_", "", rownames(u))))

# Explore the different components of list
class(pb_s) # list
str(pb_s) 
head(pb_s[[1]]) # each component (per cluster) is a sparse Matrix of class "dgCMatrix": gene (row) x sample (col)

# Print out the table of cells in each cluster-sample group
options(width = 100)
kable(table(sce$cluster_id, sce$orig.ident))
# |   | chow.eWAT| chow.iWAT| HFD.eWAT| HFD.iWAT|
# |:--|---------:|---------:|--------:|--------:|
# |0  |      1378|       835|       56|      378|
# |1  |      1041|       463|      189|      538|
# |2  |       907|       672|      179|      412|
# |3  |        27|        35|      197|      122|
# |4  |       521|       187|       57|      304|
# |5  |       525|       187|      215|      147|



###################################
##  EdgeR DIFFERENTIAL ANALYSIS  ##  
###################################

# Construct design & contrast matrix
design <- model.matrix(~ 0 + ei$orig.ident) %>% 
  set_rownames(ei$orig.ident) %>% 
  set_colnames(levels(factor(ei$orig.ident)))

contrast <- makeContrasts("HFD.eWAT-chow.eWAT", 
                          "HFD.iWAT-chow.iWAT",
                          "chow.eWAT-chow.iWAT",
                          "HFD.eWAT-HFD.iWAT",
                          levels = design)

# In each cluster, run edgeR with default parameters
res_epi <- lapply(names(pb_s), function(k) {
  y <- pb_s[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  #y <- estimateDisp(y, design)
  #fit <- glmQLFit(y, design) 
  #fit <- glmQLFTest(fit, contrast = contrast)
  fit <- glmFit(y, design, dispersion = 0.09) # values from female TRAP & male H3K27ac CnT'
  lrt <- glmLRT(fit, contrast = contrast[,1])
  topTags(lrt, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})

res_ing <- lapply(names(pb_s), function(k) {
  y <- pb_s[[k]]
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  #y <- estimateDisp(y, design)
  #fit <- glmQLFit(y, design) 
  #fit <- glmQLFTest(fit, contrast = contrast)
  fit <- glmFit(y, design, dispersion = 0.09) # values from female TRAP & male H3K27ac CnT'
  lrt <- glmLRT(fit, contrast = contrast[,2])
  topTags(lrt, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)
})

# Filter the results: FDR < 0.05 and |logFC| > 1
res_epi_sig <- lapply(res_epi, 
                      function(u)  u %>% 
                        dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                        dplyr::arrange(logFC))

res_ing_sig <- lapply(res_ing, 
                      function(u)  u %>% 
                        dplyr::filter(p_adj < 0.05, abs(logFC) > 1) %>% 
                        dplyr::arrange(logFC))
str(res_epi_sig)


###########################
##  EXPORT THE RESULTS   ##  
###########################

for(cluster in 1:length(names(pb_s))){
  filePath <- paste0("./Results/pseudobulk/Adipocytes_Epi_re-Cluster", names(pb_s)[cluster])
  out <- res_epi[[cluster]]
  write.table(out, file = paste0(filePath, "_HFD_vs_chow_results.txt"), sep = "\t")
}

for(cluster in 1:length(names(pb_s))){
  filePath <- paste0("./Results/pseudobulk/Adipocytes_Ing_re-Cluster", names(pb_s)[cluster])
  out <- res_ing[[cluster]]
  write.table(out, file = paste0(filePath, "_HFD_vs_chow_results.txt"), sep = "\t")
}



##################################################
##  PREPARATION OF logCPM MATRIX FOR MORPHEUS   ##  
##################################################

# Extract DEG lists 
c0_epi <- rownames(res_epi_sig[[1]]) #1330 
c1_epi <- rownames(res_epi_sig[[2]]) #1827 
c2_epi <- rownames(res_epi_sig[[3]]) #1717
c3_epi <- rownames(res_epi_sig[[4]]) #1197 
c4_epi <- rownames(res_epi_sig[[5]]) #2240  
c5_epi <- rownames(res_epi_sig[[6]]) #1649  

c0_ing <- rownames(res_ing_sig[[1]]) #412
c1_ing <- rownames(res_ing_sig[[2]]) #1180
c2_ing <- rownames(res_ing_sig[[3]]) #720
c3_ing <- rownames(res_ing_sig[[4]]) #310
c4_ing <- rownames(res_ing_sig[[5]]) #118
c5_ing <- rownames(res_ing_sig[[6]]) #435

# Retrieve the union lists of genes
deg_all_union <- unique(c(c0_epi, c1_epi, c2_epi, c3_epi, c4_epi, c5_epi,  
                          c0_ing, c1_ing, c2_ing, c3_ing, c4_ing, c5_ing)) # n=4236

# Make a function for extracting logCPM matrix across all subclusters
cal_logcpm <- function(u) {
  u <- t(u)
  u <- DGEList(u, remove.zeros = TRUE)
  u <- calcNormFactors(u)
  logcpm <- cpm(u, log = TRUE)
}

logcpm <- cal_logcpm(pb)
logcpm <- data.frame(logcpm)
logcpm[1:8, 1:8]
dim(logcpm) # 12766    24


# Filter the logCPM matrix by DEGs
logcpm_sig <- logcpm %>% dplyr::filter(row.names(.) %in% deg_all_union)
dim(logcpm_sig)  # 4236   24
write.table(logcpm_sig, file = "./Results/pseudobulk/Adipocytes_logcpm_deg_ACROSS_re-clusters.txt", sep = "\t", row.names = T)

