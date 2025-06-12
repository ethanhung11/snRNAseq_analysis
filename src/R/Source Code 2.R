####################################################################################################
# Title: Source Code File 2
# Content: Adipocyte-specific analysis - continued from the Source Code File 1
####################################################################################################

# Load all required libraries/packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scCustomize)


# Visualize adipocyte clusters 
Idents(object = Adipocyte.filtered) <- "integrated_snn_res.0.6"
## Plot 1) Dimension/UMAP Plot
DimPlot_scCustom(Adipocyte.filtered, colors_use = c("#f6252f", "#fe3ffa", "#f8a19f"), split.by = "orig.ident", num_columns = 4, label = F, split_seurat = T)
## Plot 2) Feature Plot
DefaultAssay(Adipocyte.filtered) <- "RNA"
Adipocyte.filtered <- NormalizeData(Adipocyte.filtered)
Adipocyte.filtered <- ScaleData(Adipocyte.filtered)
FeaturePlot_scCustom(seurat_object = Adipocyte.filtered, features = "Pparg", order = T) # markers of interest can be visualized



#####################
##  RE-CLUSTERING  ##
#####################
# Re-cluster adipocytes into additional clusters
DefaultAssay(Adipocyte.filtered) <- "integrated"
Adipocyte.clustered <- FindNeighbors(object = Adipocyte.filtered, dims = 1:40)
Adipocyte.clustered <- FindClusters(object = Adipocyte.clustered, resolution = c(0.1, 0.2, 0.3, 0.4))

# Visualize using a UMAP plot
Idents(object = Adipocyte.clustered) <- "integrated_snn_res.0.3"
Adipocyte.clustered <- RenameIdents(object = Adipocyte.clustered,
                                    "0" = "Ad1",
                                    "1" = "Ad3",
                                    "2" = "Ad2",
                                    "3" = "Ad6", 
                                    "4" = "Ad5", 
                                    "5" = "Ad4") 
levels(Adipocyte.clustered) <- c("Ad1", "Ad2", "Ad3", "Ad4", "Ad5", "Ad6")
adtype_colors <- c("#ffcc99", "#f98b88", "#f6252f", "#fe3ffa" ,"#cdb0ee","#683b7a")

DimPlot_scCustom(Adipocyte.clustered, figure_plot = T, colors_use = adtype_colors, label = F)
DimPlot_scCustom(Adipocyte.clustered, colors_use = adtype_colors, split.by = "orig.ident", num_columns = 2, label = F, split_seurat = T) + theme(legend.position = "none")



################################################
##  RELATIVE PROPORTIONS BY CELL SUBCLUSTERS  ##
################################################
n_cells <- FetchData(Adipocyte.clustered, vars = c("ident", "orig.ident")) %>%
  dplyr::count(ident, orig.ident) %>%
  tidyr::spread(ident, n)

n_cells <- t(n_cells)
colnames(n_cells) <- n_cells[1, ]
n_cells <- as.data.frame(n_cells[-1, ])
n_cells <- dplyr::mutate_all(n_cells, function(x) as.numeric(x))

n_cells[is.na(n_cells)] <- 0

n_cells <- n_cells %>%
  dplyr::mutate(chow_eWAT_perc = round(chow_eWAT/sum(chow_eWAT)*100, 2)) %>%
  dplyr::mutate(chow_iWAT_perc = round(chow_iWAT/sum(chow_iWAT)*100, 2)) %>%
  dplyr::mutate(HFD_eWAT_perc = round(HFD_eWAT/sum(HFD_eWAT)*100, 2)) %>%
  dplyr::mutate(HFD_iWAT_perc = round(HFD_iWAT/sum(HFD_iWAT)*100, 2))



#############################################
##  MARKER IDENTIFICATION FOR SUBCLUSTERS  ##
#############################################
# Identify top markers of subclusters
Ad.markers <- FindAllMarkers(Adipocyte.clustered)

# Visualize top/signature markers 
DefaultAssay(Adipocyte.clustered) <- "RNA"
Adipocyte.clustered <- NormalizeData(Adipocyte.clustered)
Adipocyte.clustered <- ScaleData(Adipocyte.clustered)

## 1) Feature plot
FeaturePlot_scCustom(seurat_object = Adipocyte.clustered, features = "Adipoq", order = T)  # Adipoq, Fam13a, Fgf14, Mageb18, Cacna1e, Cfd, B2m
## 2) Stacked violin plot
Ad_features <- c("Adipoq", "Fam13a", "Fgf14", "Mageb18", "Cacna1e", "Cfd", "B2m")
Stacked_VlnPlot(seurat_object = Adipocyte.clustered, features = Ad_features, colors_use = adtype_colors, x_lab_rotate = T , plot_spacing = 0.05)



####################################################
##  VISUALIZATION OF PSEUDOBULK ANALYSIS RESULTS  ##  # Pseudobulk analysis was performed in a separate R script
####################################################
# Generate dot plots for pathways enriched in the top markers (DEGs) of each subcluster
## 1) Ad1 & Ad2
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Sh2b2", "Pik3r1", "Insr", "Irs1")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Pparg", "Meis2", "Zeb1", "Heg1", "Sox6", "Sox5")), x_lab_rotate = F, flip_axes = T)
## 2) Ad3 & Ad4
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Vcl", "Sntb1", "Magi1", "Nedd9", "Bcar3")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Cacna1e", "Sbno2", "Ddr2", "Chrdl1", "Col5a2", "Runx2", "Runx1")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Trhde", "Adam12", "Col15a1", "Col5a3", "Adamts9", "Adamts2")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Thbs1", "Svep1", "Col5a2", "Ccn1", "Adamts5")), x_lab_rotate = F, flip_axes = T)
## 2) Ad5 & Ad6
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Hspe1", "Hspa8", "Hsp90ab1")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Fth1", "Sod2", "Prdx5", "Gpx1", "Gpx4")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Acta2", "Arpc1b", "Actb", "Actg1", "Vim", "Tuba1b", "Tuba1a")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Cd74", "H2-D1", "H2-K1", "B2m")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = Adipocyte.clustered, features = rev(c("Hp", "Rbp4", "Retn", "Lep")), x_lab_rotate = F, flip_axes = T)

