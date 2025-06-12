####################################################################################################
# Title: Source Code File 1
# Content: Data integration among our four samples, followed by clustering and cell type annotation
####################################################################################################

# Load all required libraries/packages
library(SingleCellExperiment)
library(Seurat)
library(hdf5r)
library(parallel) 
library(DoubletFinder) 
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scCustomize)


#################
##  LOAD DATA  ##
#################

# Set working directory
setwd(paste0(getwd()))

# Load the CellBender-filtered CellRanger output files and create a Seurat object for each sample
for (file in c("chow_eWAT", "chow_iWAT", "HFD_eWAT", "HFD_iWAT")){
  seurat_data <- Read10X_h5(paste0("./", file, ".output_filtered.h5"), use.names = T, unique.features=T)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

# Merge Seurat objects 
merged_seurat <- merge(x = chow_eWAT, y = c(chow_iWAT, HFD_eWAT, HFD_iWAT), add.cell.id = c("CE", "CI", "HE", "HI"))



#######################
##  QUALITY CONTROL  ##
#######################

# Add quality metrics to the Seurat object metadata
## 1) Number of genes per UMI
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
## 2) Mitochondrial ratio
merged_seurat$percent.mito <- PercentageFeatureSet(object = merged_seurat, pattern = "^mt-")

# Update the metadata with sample grouping conditions
metadata <- merged_seurat@meta.data
## Add cell IDs to metadata
metadata$cells <- rownames(metadata)
## Create columns for grouping conditions
metadata <- metadata %>% 
  mutate(depot = case_when(orig.ident == "chow_eWAT" | orig.ident == "HFD_eWAT" ~ 'eWAT',
                           orig.ident == "chow_iWAT" | orig.ident == "HFD_iWAT" ~ "iWAT"))
metadata <- metadata %>% 
  mutate(diet = case_when(orig.ident == "chow_eWAT" | orig.ident == "chow_iWAT" ~ 'chow',
                          orig.ident == "HFD_eWAT" | orig.ident == "HFD_iWAT" ~ "HFD"))

merged_seurat@meta.data <- metadata

# Filter data at the "cell level"
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nCount_RNA >= 500) & 
                            (nFeature_RNA >= 250) & 
                            (log10GenesPerUMI > 0.80) & 
                            (percent.mito < 15))

# Filter data at the "gene level"
## Extract counts
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
## Filter out mitochondrial genes 
counts_nonmito <- counts[!grepl("mt-", rownames(counts)), ]
## Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts_nonmito > 0
## Sum all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
## Keep the genes expressed in more than 10 cells
filtered_counts <- counts_nonmito[keep_genes, ]

filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)



#######################
##  DOUBLET REMOVAL  ##
#######################

table(filtered_seurat$orig.ident)
split_seurat <- SplitObject(filtered_seurat, split.by = "orig.ident")

# Run the Doublet Finder using a loop
for (i in 1:length(split_seurat)) {
  # Print the sample we are on
  print(paste0("sample ", i))
  
  # Pre-process seurat object with standard seurat workflow
  our.sample <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  our.sample <- FindVariableFeatures(our.sample)
  our.sample <- ScaleData(our.sample)
  our.sample <- RunPCA(our.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- our.sample[["pca"]]@stdev
  sum.stdv <- sum(our.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # Finish pre-processing 
  our.sample <- RunUMAP(our.sample, dims = 1:min.pc)
  our.sample <- FindNeighbors(object = our.sample, dims = 1:min.pc)              
  our.sample <- FindClusters(object = our.sample, resolution = 0.8)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep_v3(our.sample, PCs = 1:min.pc, num.cores = detectCores() - 1)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- our.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(0.075 * nrow(our.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  our.sample <- doubletFinder_v3(seu = our.sample, 
                                 PCs = 1:min.pc, 
                                 pK = optimal.pk,
                                 nExp = nExp.poi.adj)
  metadata <- our.sample@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  our.sample@meta.data <- metadata 
  
  # subset and save
  our.singlets <- subset(our.sample, doublet_finder == "Singlet")
  split_seurat[[i]] <- our.singlets
  remove(our.singlets)
}

# Converge the split objects
singlets_seurat <- merge(x = split_seurat[[1]], 
                         y = c(split_seurat[[2]], split_seurat[[3]], split_seurat[[4]]),
                         project = "Roh adipose snRNAseq")



#####################
##  NORMALIZATION  ##
#####################

# Check if integration is necessary by visualizing condition-specific clustering
singlets_seurat <- SCTransform(singlets_seurat, vars.to.regress = c("percent.mito"))
singlets_seurat <- RunPCA(singlets_seurat)
singlets_seurat <- RunUMAP(singlets_seurat, dims=1:40, reduction="pca")
DimPlot(singlets_seurat, reduction="umap", group.by = "orig.ident")

# Perform SCTransform()
split_seurat <- SplitObject(singlets_seurat, split.by = "orig.ident")

for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("percent.mito"))
}



###################
##  INTEGRATION  ##
###################

# Select the most variable features to use for integration 
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = integ_features)

# Find best buddies
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions
integrated_seurat <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# Visualize UMAP 
integrated_seurat <- RunPCA(object = integrated_seurat)
integrated_seurat <- RunUMAP(integrated_seurat, dims = 1:40, reduction = "pca")
DimPlot(integrated_seurat, group.by = "orig.ident", split.by = "orig.ident", ncol = 2)



##################
##  CLUSTERING  ##
##################

# Determine the K-nearest neighbor graph
integrated_seurat <- FindNeighbors(object = integrated_seurat, dims = 1:40)

# Determine the clusters for various resolutions (res 0.4-1.4 for datasets of 3000-5000 cells)                               
integrated_seurat <- FindClusters(object = integrated_seurat, resolution = c(0.2, 0.4, 0.6))

# Assign cell identity classes
Idents(object = integrated_seurat) <- "integrated_snn_res.0.6"   # choose the resolution that best fits

# Plot UMAP
DimPlot(integrated_seurat, reduction = "umap", label = TRUE, label.size = 4)



################################
##  CELL TYPE IDENTIFICATION  ##
################################

# First, scale data in the RNA assay
DefaultAssay(integrated_seurat) <- "RNA"
integrated_seurat <- NormalizeData(integrated_seurat)
integrated_seurat <- ScaleData(integrated_seurat)

# Method 1) Find top markers for each of the identity classes and visualize them using dot plot
Idents(object = integrated_seurat) <- "integrated_snn_res.0.6"  # choose the resolution that best fits
all.markers <- FindAllMarkers(object = integrated_seurat)
top3.markers <- Extract_Top_Markers(marker_dataframe = all.markers, num_genes = 3, named_vector = FALSE, make_unique = TRUE)
top3.markers.ordered <- c("Tenm4", "Acsl1", "Ghr", "Dcn", "Col3a1", "Fbn1", "Gpm6a", "Il1rapl1", "Nxph1", "Etl4", "Mecom", "Cdh13", "Cacna1c", "Abcc9", "Notch3", 
                          "Lyz2", "F13a1", "Retnla", "Irf8", "Wdfy4", "Cadm1", "Mcpt4", "Cpa3", "Cma1", "Skap1", "Gm2682", "Ccl5", "Ighm", "Bank1", "Igkc")
DotPlot_scCustom(seurat_object = integrated_seurat, features = rev(top3.markers.ordered), x_lab_rotate = F, flip_axes = T)

# Method 2) Visualize the known cell type-specific markers across clusters
# 2-1): Feature Plots (example)
FeaturePlot_scCustom(seurat_object = integrated_seurat, features = "Adipoq", split.by = "orig.ident", num_columns = 2, order=T)  # Adipoq for adipocytes
# 2-2): Stacked violin plots
celltype_features <- c("Adipoq", "Pdgfra", "Upk3b", "Cdh5", "Rgs5", "Adgre1", "Flt3", "Cpa3", "Skap1", "Igkc")
Stacked_VlnPlot(seurat_object = integrated_seurat, features = celltype_features, x_lab_rotate = T , plot_spacing = 0.1)



########################################
##  DATA SUBSETTING BY CELL TYPE  (1) ##
########################################

# Extract UMAP coordinates of individual cells and add them to the metadata for further filtering 
umap_coor <- integrated_seurat[["umap"]]@cell.embeddings
integrated_seurat$umap1 <- umap_coor$UMAP_1
integrated_seurat$umap2 <- umap_coor$UMAP_2


# Rename the cell identity classes as cell types
integrated_seurat <- RenameIdents(object = integrated_seurat,
                                  "0" = "APC",
                                  "1" = "Macrophage",
                                  "2" = "Adipocyte",
                                  "3" = "Adipocyte",
                                  "4" = "Endothelial",
                                  "5" = "Macrophage",
                                  "6" = "Macrophage",
                                  "7" = "APC",
                                  "8" = "APC",
                                  "9" = "Macrophage",
                                  "10" = "Mesothelial",
                                  "11" = "PC/SMC",
                                  "12" = "APC",
                                  "13" = "Adipocyte",
                                  "14" = "T cell",
                                  "15" = "Endothelial",
                                  "16" = "Macrophage",
                                  "17" = "Endothelial",
                                  "18" = "Endothelial",
                                  "19" = "B cell",
                                  "20" = "Endothelial",
                                  "21" = "Dendritic cell",
                                  "22" = "Mast cell",
                                  "23" = "Macrophage",
                                  "24" = "Macrophage") 

levels(integrated_seurat) <- c("Adipocyte", "APC", "Mesothelial", "Endothelial", "PC/SMC", "Macrophage", "Dendritic cell", "Mast cell", "T cell", "B cell")
celltype_colors <- c("#df5612", "#f2a331", "#f8e048", "#70b985", "#c3f1ca", "#b8d0e6", "#6c9dbb", "#256d93", "#ad84d9", "#d54a7e")

Stacked_VlnPlot(seurat_object = integrated_seurat, features = celltype_features, colors_use = celltype_colors, x_lab_rotate = T , plot_spacing = 0.1)
DimPlot_scCustom(seurat_object = integrated_seurat, colors_use = celltype_colors, label = T, figure_plot = T) 
DimPlot_scCustom(seurat_object = integrated_seurat, colors_use = celltype_colors, split.by = "orig.ident", num_columns = 2, label = F) 
table(Idents(integrated_seurat))

# Calculate cell type relative proportions (%)
n_cells_ann <- FetchData(integrated_seurat, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

n_cells_ann <- t(n_cells_ann)
colnames(n_cells_ann) <- n_cells_ann[1, ]
n_cells_ann <- as.data.frame(n_cells_ann[-1, ])
n_cells_ann <- dplyr::mutate_all(n_cells_ann, function(x) as.numeric(x))
n_cells_ann[is.na(n_cells_ann)] <- 0
n_cells_ann <- n_cells_ann %>%
  dplyr::mutate(chow_eWAT_perc = round(chow_eWAT/sum(chow_eWAT)*100, 2)) %>%
  dplyr::mutate(chow_iWAT_perc = round(chow_iWAT/sum(chow_iWAT)*100, 2)) %>%
  dplyr::mutate(HFD_eWAT_perc = round(HFD_eWAT/sum(HFD_eWAT)*100, 2)) %>%
  dplyr::mutate(HFD_iWAT_perc = round(HFD_iWAT/sum(HFD_iWAT)*100, 2))

# Subset data by cell type annotation
Adipocyte <- subset(integrated_seurat, idents = "Adipocyte")
APC <- subset(integrated_seurat, idents = "APC")
Macrophage <- subset(integrated_seurat, idents = "Macrophage")
Endothelial <- subset(integrated_seurat, idents = c("Endothelial", "Pericyte"))
Lymphocyte <- subset(integrated_seurat, idents = c("T cell", "B cell"))



################################################################
##  DATA SUBSETTING BY CELL TYPE  (2): FURTHER FILTERING ##
################################################################

# Filter out presumed doublets from each cell type Seurat object
all <- DimPlot_scCustom(integrated_seurat)+theme(legend.position = "none") + 
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2))
ggplot_build(all)$layout$panel_params[[1]]$x.range #[1] -14.06032  13.00683
ggplot_build(all)$layout$panel_params[[1]]$y.range #[1] -16.05151  12.32887

# 1) Adipocyte
DimPlot_scCustom(Adipocyte) + 
  coord_cartesian(xlim=c(-14.06032, 13.00683),ylim=c(-16.05151, 12.32887)) + 
  scale_x_continuous(breaks=seq(-10,10,5))+
  scale_y_continuous(breaks=seq(-15,10,5))+
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.y = element_line(color = "gray50", size = 0.25, linetype = 2))

Adipocyte.filtered <- subset(x = Adipocyte, subset = (umap2 > 1))

# 2) APC - no filtering needed
DimPlot_scCustom(APC) + 
  coord_cartesian(xlim=c(-14.06032, 13.00683),ylim=c(-16.05151, 12.32887)) + 
  scale_x_continuous(breaks=seq(-10,10,5))+
  scale_y_continuous(breaks=seq(-15,10,5))+
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2)) +  
  theme(panel.grid.minor.y = element_line(color = "gray50", size = 0.25, linetype = 2)) 


# 3) Macrophage
DimPlot_scCustom(Macrophage) + 
  coord_cartesian(xlim=c(-14.06032, 13.00683),ylim=c(-16.05151, 12.32887)) + 
  scale_x_continuous(breaks=seq(-10,10,5))+
  scale_y_continuous(breaks=seq(-15,10,5))+
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.minor.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.y = element_line(color = "gray50", size = 0.25, linetype = 2)) 

Mac.filtered <- subset(x = Macrophage, subset = (umap1 > -2.5) & (umap2 < -2.5))

# 4) Lymphocyte
DimPlot_scCustom(Lymphocyte) + 
  coord_cartesian(xlim=c(-14.06032, 13.00683),ylim=c(-16.05151, 12.32887)) + 
  scale_x_continuous(breaks=seq(-10,10,5))+
  scale_y_continuous(breaks=seq(-15,10,5))+
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.minor.y = element_line(color = "gray50", size = 0.25, linetype = 2)) 

Lym.filtered <- subset(x = Lymphocyte, subset = (umap1 < -5))

# 4) Endothelial cells
DimPlot_scCustom(Endothelial) + 
  coord_cartesian(xlim=c(-14.06032, 13.00683),ylim=c(-16.05151, 12.32887)) + 
  scale_x_continuous(breaks=seq(-10,10,5))+
  scale_y_continuous(breaks=seq(-15,10,5))+
  theme(legend.position = "none") +
  theme(panel.grid.major.x = element_line(color = "gray50", size = 0.25, linetype = 2)) + 
  theme(panel.grid.minor.x = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.major.y = element_line(color = "gray50", size = 0.25, linetype = 2)) +
  theme(panel.grid.minor.y = element_line(color = "gray50", size = 0.25, linetype = 2)) 

Endo.filtered <- subset(x = Endothelial2, subset = (umap1 > -7.5 & umap1 < 5) & (umap2 > -6.5 & umap2 <5))
Endo.filtered <- subset(x = Endo2.filtered, subset = (umap1 < 2.5) | (umap1 >= 2.5 & umap2 > 2.5))
