########################################################################################################################
# Title: Source Code File 6
# Content: Comparison to previously published sc/snRNA-seq dataset (Emont et al.) - continued from Source Code File 1
########################################################################################################################

# Load all required libraries/packages
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(scCustomize)
library(harmony)
library(glmGamPoi)
library(plyr)


#################
##  LOAD DATA  ##
#################

# Set working directory
setwd(paste0(getwd()))

# Load the subsetted Seurat objects for individual cell types
emont_global <- readRDS("mouse_all.rds") # downloaded from here: https://singlecell.broadinstitute.org/single_cell/study/SCP1376/a-single-cell-atlas-of-human-and-mouse-white-adipose-tissue#study-summary
so_global <- integrated_seurat



################################
##  SCTransform NORMALIZATION ##
################################

list <- list(so_global, emont_global)
list <- lapply(X = list, FUN = SCTransform)



#########################
## HARMONY INTEGRATION ##
#########################

# Select the most variable features to use for integration (default = the top 2000 genes)
integ_features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000) 

# Merge split (normalized) objects
integrated_global <- merge(x = list[[1]], 
                           y = list[[2]], 
                           merge.data = T)
DefaultAssay(integrated_global) <- "SCT"

# Manually set variable features of merged Seurat object
VariableFeatures(integrated_global) <- integ_features

# Calculate PCs using manually set variable features 
integrated_global <- RunPCA(integrated_global, assay = "SCT", npcs = 50)

#... Add dataset info as a metadata column
metadata <- integrated_global@meta.data
metadata <- metadata %>% mutate(dataset = if_else(orig.ident == "Mm", "Emont", "So"))
integrated_global@meta.data <- metadata

# Make sure if seurat obj contains one (or more) variable(s) describing the "factor(s)" we want to integrate on
integrated_global$dataset <- factor(integrated_global$dataset, levels = c("So", "Emont"))
integrated_global$orig.ident <- factor(integrated_global$orig.ident, levels = c("epi_WT", "epi_KO", "ing_WT", "ing_KO", "Mm"))

# Run Harmony!
harmonized_global <- RunHarmony(integrated_global,
                                group.by.vars = "dataset",
                                reduction = "pca",
                                assay.use = "SCT",
                                reduction.save = "harmony")


# UMAP visualization: PCA or UMAP
# Run PCA and UMAP
harmonized_global <- RunUMAP(harmonized_global, dims = 1:40, reduction = "harmony", assay = "SCT")
harmonized_global <- RunTSNE(object = harmonized_global)



#######################
## METADATA REFINING ##
#######################

metadata <- harmonized_global@meta.data   # 251,116 x 52
metadata_so <- so_global@meta.data        #  53,395 x 22
metadata_emont <- emont_global@meta.data  # 197,721 x 37


# Refine metadata columns of each dataset (to make them consistent between the datasets)
metadata_so_short <- metadata_so[c(1:3, 5:8, 19, 22)]
metadata_so_short <- metadata_so_short %>% mutate(depot = recode(depot, 'eWAT' = 'PG', 'iWAT' = 'ING'))
metadata_so_short <- metadata_so_short %>% mutate(diet = recode(diet, 'chow' = 'Chow'))
metadata_so_short <- metadata_so_short %>% rename(sample = orig.ident)
metadata_so_short <- metadata_so_short %>% rename(mt.percent = percent.mito)
metadata_so_short <- metadata_so_short %>% rename(celltype_so = celltype_res0.6)

metadata_emont_short <- metadata_emont[c(2:3, 8, 10, 12, 24, 34:35)]
metadata_emont_short$cells <- rownames(metadata_emont_short)
metadata_emont_short <- metadata_emont_short %>% rename(celltype_emont1 = cell_type,
                                                        celltype_emont2 = cell_type2)

# Combine two metadata dataframes
metadata_global_short <- rbind.fill(metadata_emont_short, metadata_so_short)

# Add a dataset column
metadata_global_short <- metadata_global_short %>% mutate(dataset = if_else(sample == "chow_eWAT"| sample == "chow_iWAT" | sample == "HFD_eWAT" | sample == "HFD_iWAT", "So", "Emont"))

# Convert grouping variables from character to factor 
cols <- c("sample", "depot", "diet", "celltype_emont1", "celltype_emont2", "celltype_so", "ad_type", "dataset")
metadata_global_short[ ,cols] <- lapply(metadata_global_short[ ,cols], factor)
metadata_global_short$dataset <- factor(metadata_global_short$dataset, levels = c("So", "Emont"))

# Set the combined DF as metadata of the integrated Seurat object
harmonized_global@meta.data <- metadata_global_short
rownames(harmonized_global@meta.data) <- harmonized_global@meta.data$cells



################
## CLUSTERING ##
################

# Cluster cells
harmonized_global <- FindNeighbors(object = harmonized_global, dims = 1:40)
harmonized_global <- FindClusters(object = harmonized_global, resolution = c(0.1, 0.2, 0.3, 0.4))
save(harmonized_global, file= "./data/RData_chow_HFD_new_CellBender/Global_so_emont_harmony-integrated_clustered.RData")

# Assign cell identity classes
Idents(object = harmonized_global) <- "SCT_snn_res.0.3"

# Plot UMAP
DimPlot_scCustom(seurat_object = harmonized_global, pt.size = 1) + theme(legend.position = "none")
DimPlot_scCustom(seurat_object = harmonized_global, split.by = "dataset") 



##########################
## CELL TYPE ANNOTATION ##
##########################

# Visualize cell type-specific markers using feature plots
DefaultAssay(harmonized_global) <- "RNA"
harmonized_global <- NormalizeData(harmonized_global)
harmonized_global <- ScaleData(harmonized_global)

FeaturePlot_scCustom(seurat_object = harmonized_global, features = "Adipoq", order=T)
# Cell type markers: Adipoq, Pdgfra, Upk3b, Cdh5, Rgs5, Adgre1, Flt3, Cpa3, Skap1, Igkc

# Annotate clusters 
Idents(object = harmonized_global) <- "SCT_snn_res.0.3"
harmonized_global <- RenameIdents(object = harmonized_global,
                                      "0" = "Macrophage", 
                                      "1" = "Adipocyte",
                                      "2" = "APC",
                                      "3" = "APC",
                                      "4" = "Macrophage",
                                      "5" = "Adipocyte",
                                      "6" = "Endothelial",
                                      "7" = "Macrophage",
                                      "8" = "Monocyte_DC",
                                      "9" = "Macrophage",
                                      "10" = "Mesothelial",
                                      "11" = "Mesothelial",
                                      "12" = "Adipocyte",
                                      "13" = "APC",
                                      "14" = "Epithelial-female",
                                      "15" = "Macrophage",
                                      "16" = "APC", 
                                      "17" = "APC",
                                      "18" = "T cell",
                                      "19" = "Pericyte_SMC",
                                      "20" = "Epithelial-male",
                                      "21" = "Epithelial-female",
                                      "22" = "B cell",
                                      "23" = "Epithelial-female",
                                      "24" = "Adipocyte",
                                      "25" = "LEC", 
                                      "26" = "Epithelial-male",
                                      "27" = "Macrophage",
                                      "28" = "Mast cell",
                                      "29" = "Epithelial-male",
                                      "30" = "APC",
                                      "31" = "APC") 

harmonized_global@active.ident <- factor(harmonized_global_ref@active.ident, 
                                             levels = c("Adipocyte", "APC", "Mesothelial", "Endothelial", "Pericyte_SMC", "LEC", 
                                                        "Macrophage", "Monocyte_DC", "Mast cell", "T cell", "B cell", "Epithelial-female", "Epithelial-male"))

DimPlot_scCustom(seurat_object = harmonized_global, pt.size = 1, split.by = "dataset", 
                 colors_use = c("#df5610", "#f2a431", "#f8e049", "#70b985", "#bbe8c2", "#d1eb30" ,
                                "#88c9f6","#256d92", "#011993", "#ae83d9", "#d64a7e", "#929292", "#424242"),
                 label = F, figure_plot = T) 


# Calculate relative proportions
n_cells <- FetchData(harmonized_global, 
                     vars = c("ident", "dataset")) %>%
  dplyr::count(ident, dataset) %>%
  tidyr::spread(ident, n)

n_cells <- t(n_cells)
colnames(n_cells) <- n_cells[1, ]
n_cells <- as.data.frame(n_cells[-1, ])
n_cells <- dplyr::mutate_all(n_cells, function(x) as.numeric(x))

n_cells[is.na(n_cells)] <- 0

n_cells <- n_cells %>%
  dplyr::mutate(So_perc = round(So/sum(So)*100, 2)) %>%
  dplyr::mutate(Emont_perc = round(Emont/sum(Emont)*100, 2)) 



#######################################
## QC METRICS BY DATASET & CELL TYPE ##
#######################################
harmonized_global@meta.data$int_celltype <- Idents(object = harmonized_global) 
metadata <- as.data.frame(harmonized_global@meta.data)
metadata <- metadata[, c(1:2, 6, 12, 18)]
metadata$int_celltype <- factor(metadata$int_celltype, levels = c("Adipocyte", "APC", "Mesothelial", 
                                                                  "Endothelial", "Pericyte_SMC", "LEC", 
                                                                  "Macrophage", "Monocyte_DC", "Mast cell",
                                                                  "T cell", "B cell", "Epithelial-female", "Epithelial-male"))

## 1) nCount_RNA: nUMI 
p1 <- ggplot(metadata, aes(x = int_celltype, y = nCount_RNA, fill = dataset)) + 
  geom_boxplot(outlier.size = 0.3) + 
  facet_wrap(~int_celltype, scale = "free", nrow = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  labs(y = "nUMI")

## 2) nFeature_RNA: nGene
p2 <- ggplot(metadata, aes(x = int_celltype, y = nFeature_RNA, fill = dataset)) + 
  geom_boxplot(outlier.size = 0.3) + 
  facet_wrap(~int_celltype, scale = "free", nrow = 2) + 
  theme_bw() + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) + 
  labs(y = "nGene")



###################################
## ADIPOCYTE-SPECIFIC COMPARISON ##
###################################

# Subset adipocytes
Adipocyte <- subset(harmonized_global, idents = "Adipocyte")

# Visualize adipocyte subclusters (So's) only
Idents(object = Adipocyte) <- "ad_type"
Ad1 <- WhichCells(object = Adipocyte, idents = "Ad1")
Ad2 <- WhichCells(object = Adipocyte, idents = "Ad2")
Ad3 <- WhichCells(object = Adipocyte, idents = "Ad3")
Ad4 <- WhichCells(object = Adipocyte, idents = "Ad4")
Ad5 <- WhichCells(object = Adipocyte, idents = "Ad5")
Ad6 <- WhichCells(object = Adipocyte, idents = "Ad6")

cells_so <- list("Ad1" = Ad1,
                 "Ad2" = Ad2,
                 "Ad3" = Ad3,
                 "Ad4" = Ad4,
                 "Ad5" = Ad5,
                 "Ad6" = Ad6)
Cell_Highlight_Plot(seurat_object = Adipocyte, 
                    #split.by = "dataset",
                    cells_highlight = cells_so, pt.size = 0.1, 
                    highlight_color = c("#ffcc99", "#f98b88", "#f6252f", "#fe3ffa" ,"#cdb0ee","#683b7a"), 
                    background_color = c(rep("lightgray", 6)))


# Visualize adipocyte subclusters (Emont's) only
Idents(object = Adipocyte) <- "celltype_emont1"
levels(droplevels(Adipocyte$celltype_emont1))
mAd1 <- WhichCells(object = Adipocyte, idents = "mAd1")
mAd2 <- WhichCells(object = Adipocyte, idents = "mAd2")
mAd3 <- WhichCells(object = Adipocyte, idents = "mAd3")
mAd4 <- WhichCells(object = Adipocyte, idents = "mAd4")
mAd5 <- WhichCells(object = Adipocyte, idents = "mAd5")
mAd6 <- WhichCells(object = Adipocyte, idents = "mAd6")

cells <- list("mAd1" = mAd1,
              "mAd2" = mAd2,
              "mAd3" = mAd3,
              "mAd4" = mAd4,
              "mAd5" = mAd5,
              "mAd6" = mAd6)
Cell_Highlight_Plot(seurat_object = Adipocyte, 
                    #split.by = "dataset",
                    cells_highlight = cells, pt.size = 0.1,
                    highlight_color = c("#fee005", "#feaf07", "#fe7f0c", "#d73832" ,"#ba1d4d","#97076b"), 
                    background_color = c(rep("lightgray", 23)))


## QC METRICS BY DATASET & ADIPOCYTE SUBCLUSTER ##
metadata_ad <- as.data.frame(Adipocyte@meta.data)
metadata_ad <- metadata_ad[, c(1:2, 6:7, 11:12)]
metadata_ad_so <- metadata_ad %>% filter(dataset == "So")
metadata_ad_emont <- metadata_ad %>% filter(dataset == "Emont")

# 1) nCount_RNA: nUMI 
table(metadata_ad_so$ad_type, useNA = "ifany")
ggplot(metadata_ad_so, aes(x = ad_type, y = nCount_RNA, fill = ad_type)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#ffcc99", "#f98b88", "#f6252f", "#fe3ffa" ,"#cdb0ee","#683b7a"))+
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(y = "nUMI") + 
  scale_y_continuous(limits = c(0, 10000), breaks=seq(0,9000,3000))

table(metadata_ad_emont$celltype_emont1, useNA = "ifany")
metadata_ad_emont$celltype_emont1 <- droplevels(metadata_ad_emont$celltype_emont1)
metadata_ad_emont <- metadata_ad_emont %>% mutate(celltype_emont1 = if_else(!celltype_emont1 %in% c("mAd1", "mAd2", "mAd3", "mAd4","mAd5", "mAd6"), NA_real_, celltype_emont1))

levels(metadata_ad_emont$celltype_emont1)[levels(metadata_ad_emont$celltype_emont1)%in% c("mASPC1", "mASPC2", "mASPC4", "mASPC5", 
                                                                                          "mBcell","mcDC1", "mcDC2", "mDC3", 
                                                                                          "mMac1", "mMac3", "mMac4", 
                                                                                          "mMast", "mMes3", "mMono1", "mMono2", 
                                                                                          "mSMC", "mTcell2")] <- NA
ggplot(metadata_ad_emont, aes(x = celltype_emont1, y = nCount_RNA, fill = celltype_emont1)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#fee005", "#feaf07", "#fe7f0c", "#d73832" ,"#ba1d4d","#97076b"))+
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(y = "nUMI") + 
  scale_y_continuous(limits = c(0, 10000), breaks=seq(0,9000,3000))

# 2) nFeature_RNA: nGene
ggplot(metadata_ad_so, aes(x = ad_type, y = nFeature_RNA, fill = ad_type)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#ffcc99", "#f98b88", "#f6252f", "#fe3ffa" ,"#cdb0ee","#683b7a"))+
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(y = "nGene") + 
  scale_y_continuous(limits = c(0, 8000), breaks=seq(0,8000,2000))

ggplot(metadata_ad_emont, aes(x = celltype_emont1, y = nFeature_RNA, fill = celltype_emont1)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values = c("#fee005", "#feaf07", "#fe7f0c", "#d73832" ,"#ba1d4d","#97076b"))+
  theme_bw() + 
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") + 
  labs(y = "nGene") + 
  scale_y_continuous(limits = c(0, 8000), breaks=seq(0,8000,2000))


## RELATIVE PROPORTION CALCULATION BY DATASET ##
# 1) So et al.
Idents(object = Adipocyte.filtered) <- "dataset"
Adipocyte.filtered.so <- subset(Adipocyte.filtered, idents = "So")
Idents(object = Adipocyte.filtered.so) <- "ad_type"

n_cells <- FetchData(Adipocyte.filtered.so, 
                     vars = c("ident", "sample")) %>%
  dplyr::count(ident, sample) %>%
  tidyr::spread(ident, n)
n_cells <- t(n_cells)
colnames(n_cells) <- n_cells[1, ]
n_cells <- as.data.frame(n_cells[-1, ])
n_cells <- dplyr::mutate_all(n_cells, function(x) as.numeric(x))

n_cells[is.na(n_cells)] <- 0

n_cells <- n_cells %>%
  dplyr::mutate(chow_eWAT_perc = round(chow_eWAT/sum(chow_eWAT)*100, 2)) %>%
  dplyr::mutate(chow_iWAT_perc = round(chow_iWAT/sum(chow_iWAT)*100, 2)) %>%
  dplyr::mutate(HFD_eWAT_perc = round(HFD_eWAT/sum(chow_eWAT)*100, 2)) %>%
  dplyr::mutate(HFD_iWAT_perc = round(HFD_iWAT/sum(chow_iWAT)*100, 2))

# 2) Emont et al.
Adipocyte.filtered.emont <- subset(Adipocyte.filtered, idents = "Emont")
Idents(object = Adipocyte.filtered.emont) <- "celltype_emont1"
table(Adipocyte.filtered.emont$celltype_emont1, useNA = "ifany")

library(stringr)
metadata_ad_emont_flt <- Adipocyte.filtered.emont@meta.data
metadata_ad_emont_flt <- metadata_ad_emont_flt %>% mutate(group = str_c(diet, depot))
metadata_ad_emont_flt <- metadata_ad_emont_flt %>% mutate(celltype_emont1 = if_else(!celltype_emont1 %in% c("mAd1", "mAd2", "mAd3", "mAd4","mAd5", "mAd6"), NA_real_, celltype_emont1))
table(metadata_ad_emont_flt$celltype_emont1, useNA = "ifany")
levels(metadata_ad_emont_flt$celltype_emont1)[levels(metadata_ad_emont_flt$celltype_emont1)%in% c("luminal_epithelial_AV", "luminal_epithelial_HS", "male_epithelial_1", "male_epithelial_2", "male_epithelial_3", "mammary_fibroblast",
                                                                                                  "mASPC1", "mASPC2", "mASPC3", "mASPC4", "mASPC5", "mASPC6", 
                                                                                                  "mBcell", "mcDC1", "mcDC2", "mDC3", "mEndoA1", "mEndoA2", "mEndoS1", "mEndoS2", "mEndoV", 
                                                                                                  "mLEC1", "mLEC2", "mMac1", "mMac2", "mMac3", "mMac4", "mMac5", "mMast",
                                                                                                  "mMes1", "mMes2", "mMes3", "mMono1", "mMono2", "mNeu", "mNK", "mPeri", "mSMC", "mTcell1", "mTcell2", "mTcell3", "myoepithelial")] <- NA
Adipocyte.filtered.emont@meta.data <- metadata_ad_emont_flt
Idents(object = Adipocyte.filtered.emont) <- "celltype_emont1"
n_cells <- FetchData(Adipocyte.filtered.emont, 
                     vars = c("ident", "group")) %>%
  dplyr::count(ident, group) %>%
  tidyr::spread(ident, n)
n_cells <- t(n_cells)
colnames(n_cells) <- n_cells[1, ]
n_cells <- as.data.frame(n_cells[-1, ])
n_cells <- dplyr::mutate_all(n_cells, function(x) as.numeric(x))

n_cells[is.na(n_cells)] <- 0

n_cells <- n_cells %>%
  dplyr::mutate(ChowING_perc = round(ChowING/sum(ChowING)*100, 2)) %>%
  dplyr::mutate(ChowPG_perc = round(ChowPG/sum(ChowPG)*100, 2)) %>%
  dplyr::mutate(HFDING_perc = round(HFDING/sum(HFDING)*100, 2)) %>%
  dplyr::mutate(HFDPG_perc = round(HFDPG/sum(HFDPG)*100, 2))

