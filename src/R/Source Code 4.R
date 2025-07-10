####################################################################################################
# Title: Source Code File 4
# Content: Pseudotime analysis using Monocle3 - continued from the Source Code File 2
####################################################################################################

# Load all required libraries/packages
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(scCustomize)
library(viridisLite)



# MONOCLE3 WORKFLOW:
##########################################
##  CONVERSION TO cell_data_set OBJECT  ##
##########################################
cds <- as.cell_data_set(Adipocyte.clustered)
cds

# Get cell metadata
colData(cds)

# Get gene data
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))

# Get count data
counts(cds)
dim(counts(cds)) # 22740 10750



#######################
##  CELL CLUSTERING  ##
#######################
# Assign all cells as a single partition
reacreate.partition <- c(rep(1,length(cds@colData@rownames))) # make a vector
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- Adipocyte.clustered@active.ident
cds@clusters$UMAP$clusters <- list_cluster

# Assign UMAP coordinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- Adipocyte.clustered@reductions$umap@cell.embeddings

# Plot
cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        show_trajectory_graph = F,
                                        group_label_size = 5) + theme(legend.position = "right")
cluster.before.trajectory



###########################
##  LEARNING TRAJECTORY  ##
###########################
cds <- learn_graph(cds, use_partition = FALSE)

trajectory <- plot_cells(cds,
                         color_cells_by = 'cluster',
                         label_groups_by_cluster = FALSE,
                         label_branch_points = FALSE,
                         label_roots = FALSE,
                         label_leaves = FALSE,
                         group_label_size = 5,
                         cell_size = 0,
                         cell_stroke = 1) +
  coord_cartesian(xlim=c(-13.5, -2.5),ylim=c(1.5, 12.5)) + 
  scale_x_continuous(breaks=seq(-12.5, 2.5, 2.5))+
  scale_y_continuous(breaks=seq(2.5, 12.5, 2.5)) 



###################################
##  CELL ORDERING BY PSEUDOTIME  ##
###################################
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == c(0,2)])) # assign root cells with cluster 0 (i.e., Ad1)

pseudoplot <- plot_cells(cds,
                         color_cells_by = 'pseudotime',
                         show_trajectory_graph = TRUE,
                         label_groups_by_cluster = TRUE,
                         label_branch_points = FALSE,
                         label_roots = FALSE,
                         label_leaves = FALSE,
                         cell_size = 0,
                         cell_stroke = 1) + 
  viridis::scale_color_viridis(option = "plasma", direction = -1) + 
  coord_cartesian(xlim=c(-13.5, -2.5),ylim=c(1.5, 12.5)) + 
  scale_x_continuous(breaks=seq(-12.5, 2.5, 2.5))+
  scale_y_continuous(breaks=seq(2.5, 12.5, 2.5)) +
  theme(legend.position = "none")

# cells ordered by monocle3 pseudotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(integrated_snn_res.0.3, monocle3_pseudotime, median), fill = integrated_snn_res.0.3)) +
  geom_boxplot()



