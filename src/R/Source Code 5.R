####################################################################################################
# Title: Source Code File 5
# Content: Other cell type-specific analyses - continued from the Source Code File 1
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


##############
##  1) APC  ##
##############
Idents(object = APC) <- "integrated_snn_res.0.6"
APC <- RenameIdents(object = APC,
                    "0" = "APC4",
                    "7" = "APC1",
                    "8" = "APC2",
                    "12" = "APC3") 

levels(APC) <- c("APC1", "APC2", "APC3", "APC4")
APC_colors <- c("#b01e68", "#1dfbce", "#aa38fe", "#5a5156")
DimPlot_scCustom(APC, colors_use = APC_colors, figure_plot = T)


##  RELATIVE PROPORTIONS BY CELL SUBCLUSTERS  ##
n_cells <- FetchData(APC, vars = c("ident", "orig.ident")) %>%
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


##  UNDERSTANDING SUBCLUSTERS  ##
DefaultAssay(APC) <- "RNA"
APC <- NormalizeData(APC)
APC <- ScaleData(APC)

# Identify top markers of subclusters
APC.markers <- FindAllMarkers(APC)

## 1) Feature plot
FeaturePlot_scCustom(seurat_object = APC, features = "Pdgfra", order = T)  

## 2) Stacked violin plot
APC_features <- c("Pdgfra", "Dpp4", "Cd55", "F3", "Lpl")
Stacked_VlnPlot(seurat_object = APC, features = APC_features, colors_use = APC_colors, x_lab_rotate = T , plot_spacing = 0.05)

## 3) Dot plots for genes of interest 
DotPlot_scCustom(seurat_object = APC, features = rev(c("Ccl2", "Tgfb2", "Loxl1", "Fn1")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = APC, features = rev(c("Ltbp4", "Chrdl1", "Vegfc", "Sfrp1", "Gdf10", "Gas6", "Ccn2")), x_lab_rotate = F, flip_axes = T)
DotPlot_scCustom(seurat_object = APC, features = rev(c("Car3", "Pparg", "Cd36", "Fabp5", "Apoe", "Fabp4", "Hsd11b1", "Lpl")), x_lab_rotate = F, flip_axes = T)




######################
##  2) MACROPHAGES  ##
######################
Mac.filtered <- RenameIdents(object = Mac.filtered,
                             "1" = "MC1",
                             "5" = "MC2",
                             "6" = "MC3",
                             "9" = "MC4", 
                             "16" = "MC5", 
                             "23" = "MC6", 
                             "24" = "MC7") 
levels(Mac.filtered) <- c("MC1", "MC2", "MC3", "MC4", "MC5", "MC6", "MC7")
Mac_colors <- c("#e4e1e3", "#3283fe", "#feaf16", "#90ad1d", "#1d8356", "#f7e1a0", "#c075a6")
DimPlot_scCustom(Mac.filtered, colors_use = Mac_colors, figure_plot = T)


##  RELATIVE PROPORTIONS BY CELL SUBCLUSTERS  ##
n_cells <- FetchData(Mac.filtered, vars = c("ident", "orig.ident")) %>%
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


##  UNDERSTANDING SUBCLUSTERS  ##
DefaultAssay(Mac.filtered) <- "RNA"
Mac.filtered <- NormalizeData(Mac.filtered)
Mac.filtered <- ScaleData(Mac.filtered)

# Identify top markers of subclusters
Mac.markers <- FindAllMarkers(Mac.filtered)

## 1) Feature plot
FeaturePlot_scCustom(seurat_object = Mac.filtered, features = "Clec10a", order = T)  # check the known markers of macrophage subtypes

## 2) Stacked violin plot
Mac_features <- c("Adgre1", "Cd209f", "Creb5", "Gpnmb", "Cd226", "Plac8", "Slpi", "Col1a2")
Stacked_VlnPlot(seurat_object = Mac.filtered, features =Mac_features, colors_use = Mac_colors, x_lab_rotate = T , plot_spacing = 0.05)

## 3) Dot plots for genes of interest 
DotPlot_scCustom(seurat_object = Mac.filtered, features = rev(cc("Pparg", "Slc37a2", "Cd36", "Plin2", "Fabp5", "Lpl")), x_lab_rotate = F, flip_axes = T)


##  IDENTIFICATION OF ENRICHED PATHWAYS IN EACH SUBCLUSTER  ##
Mac.top100 <- Mac.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
Mac.top100pval <- subset(Mac.top100, rowSums( Mac.top100[5] < 0.05) > 0)
Mac.df <- Mac.top100pval[ , 7:6]
Mac.df <- split(Mac.df$gene, Mac.df$cluster)

# Convert gene symbol to entrez ID
Mac.df$`1` <- bitr(Mac.df$`1`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`5` <- bitr(Mac.df$`5`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`6` <- bitr(Mac.df$`6`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`9` <- bitr(Mac.df$`9`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`16` <- bitr(Mac.df$`16`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`23` <- bitr(Mac.df$`23`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")
Mac.df$`24` <- bitr(Mac.df$`24`, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db")

# Make a list of genes: each element should be named!
Mac.genelist <- list(PVM1 = Mac.df$`1`$ENTREZID,
                     PVM2 = Mac.df$`5`$ENTREZID,
                     LAM = Mac.df$`6`$ENTREZID,
                     NPVM = Mac.df$`9`$ENTREZID,
                     Mono = Mac.df$`16`$ENTREZID,
                     Eff = Mac.df$`23`$ENTREZID,
                     CEM = Mac.df$`24`$ENTREZID)

Mac.GOplot <- compareCluster(geneCluster = Mac.genelist, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
dotplot(Mac.GOplot)



######################
##  3) LYMPHOCYTES  ##
######################
# Visualize lymphocyte clusters 
Idents(object = Lym.filtered) <- "integrated_snn_res.0.6"
## Plot 1) Dimension/UMAP Plot
DimPlot_scCustom(Lym.filtered, colors_use = c("#315a9b", "#fbe426"), split.by = "orig.ident", num_columns = 4, label = F, split_seurat = T)
## Plot 2) Feature Plot
DefaultAssay(Lym.filtered) <- "RNA"
Lym.filtered <- NormalizeData(Lym.filtered)
Lym.filtered <- ScaleData(Lym.filtered)
FeaturePlot_scCustom(seurat_object = Lym.filtered, features = "Cd4", order = T) # markers of interest can be visualized

# Re-cluster lymphocytes into additional clusters
DefaultAssay(Lym.filtered) <- "integrated"
Lym.clustered <- FindNeighbors(object = Lym.filtered, dims = 1:40)
Lym.clustered <- FindClusters(object = Lym.clustered, resolution = c(0.1, 0.2, 0.4, 0.6, 0.8))

# Visualize using a UMAP plot
Idents(object = Lym.clustered) <- "integrated_snn_res.0.4"
Lym.clustered <- RenameIdents(object = Lym.clustered,
                              "0" = "BC1",
                              "1" = "Treg",
                              "2" = "CD4+ T",
                              "3" = "NKT", 
                              "4" = "CD8+ T", 
                              "5" = "BC2") 
levels(Lym.clustered) <- c("CD4+ T", "CD8+ T", "Treg", "NKT", "BC1", "BC2")
Lym_colors <- c("#315a9b","#6b8cd2", "#ddf9ff","#a3c1ff", "#fbe426", "#cb7600")
DimPlot_scCustom(Lym.clustered, figure_plot = T, colors_use = Lym_colors, label = F)
DimPlot_scCustom(Lym.clustered, colors_use = Lym_colors, split.by = "orig.ident", num_columns = 2, label = F, split_seurat = T) + theme(legend.position = "none")


##  RELATIVE PROPORTIONS BY CELL SUBCLUSTERS  ##
n_cells <- FetchData(Lym.clustered, vars = c("ident", "orig.ident")) %>%
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


##  UNDERSTANDING SUBCLUSTERS  ##
# Visualize top/signature markers 
DefaultAssay(Lym.clustered) <- "RNA"
Lym.clustered <- NormalizeData(Lym.clustered)
Lym.clustered <- ScaleData(Lym.clustered)

## 1) Feature plot
FeaturePlot_scCustom(seurat_object = Lym.clustered, features = "Prkcq", order = T)
## 2) Stacked violin plot
Lym_features <- c("Cd4", "Cd8a", "Il2ra", "Klra8", "Bcl11a", "Prkar2b") 
Stacked_VlnPlot(seurat_object = Lym.clustered, features = Lym_features, colors_use = Lym_colors, x_lab_rotate = T , plot_spacing = 0.05)



#########################
##  4) VASCULAR CELLS  ##
#########################
# Visualize lymphocyte clusters 
Idents(object = Lym.filtered) <- "integrated_snn_res.0.6"


## Plot 1) Dimension/UMAP Plot
DimPlot_scCustom(Endo.filtered, colors_use = c("#19f932", "#dea0fd", "#c03d35", "#85660d","#b127a1", "#1cbe4f"), 
                 split.by = "orig.ident", num_columns = 4, label = F, split_seurat = T)
## Plot 2) Feature Plot
DefaultAssay(Endo.filtered) <- "RNA"
Endo.filtered <- NormalizeData(Endo.filtered)
Endo.filtered <- ScaleData(Endo.filtered)
FeaturePlot_scCustom(seurat_object = Endo.filtered, features = "Pecam1", order = T) # markers of interest can be visualized

# Re-cluster lymphocytes into additional clusters
DefaultAssay(Endo.filtered) <- "integrated"
Endo.clustered <- FindNeighbors(object = Endo.filtered, dims = 1:40)
Endo.clustered <- FindClusters(object = Endo.clustered, resolution = c(0.1, 0.15, 0.2, 0.3, 0.4))

# Visualize using a UMAP plot
Idents(object = Endo.clustered) <- "integrated_snn_res.0.3"
Endo.clustered <- RenameIdents(object = Endo.clustered,
                               "0" = "VcapEC",
                               "1" = "AcapEC",
                               "2" = "PC",
                               "3" = "AdEC", 
                               "4" = "AngEC", 
                               "5" = "ArtEC",
                               "6" = "VenEC",
                               "7" = "SMC",
                               "8" = "LEC") 
levels(Endo.clustered) <- c("VenEC", "ArtEC", "VcapEC", "AcapEC", "AngEC", "AdEC", "PC", "SMC", "LEC")
Endo_colors <- c("#b127a1","#1cbe4f", "#19f932", "#ffb676", "#85660d", "#c03d35", "#dea0fd", "#f884c6", "#8db4dd")
DimPlot_scCustom(Endo.clustered, figure_plot = T, colors_use = Endo_colors, label = F)


##  RELATIVE PROPORTIONS BY CELL SUBCLUSTERS  ##
n_cells <- FetchData(Endo.clustered, vars = c("ident", "orig.ident")) %>%
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


##  UNDERSTANDING SUBCLUSTERS  ##
# Visualize top/signature markers 
DefaultAssay(Endo.clustered) <- "RNA"
Endo.clustered <- NormalizeData(Endo.clustered)
Endo.clustered <- ScaleData(Endo.clustered)

## 1) Feature plot
FeaturePlot_scCustom(seurat_object = Endo.clustered, features = "Prkcq", order = T)
## 2) Stacked violin plot
Endo_features <- c("Pecam1", "Selp", "Gja5", "Car4", "Sema3g", "Pgf", "Scd1", "Kcnq5", "Dgkb", "Mmrn1")
Stacked_VlnPlot(seurat_object = Endo.clustered, features = Endo_features, colors_use = Endo_colors, x_lab_rotate = T , plot_spacing = 0.05)

