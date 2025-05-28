library(here)
library(dplyr)
library(ggplot2)

library(hdf5r)
library(patchwork)
library(Seurat)

#' Integrate Data, across conditions
#' v5 integration: <https://satijalab.org/seurat/articles/seurat5_integration>, compiled Oct 2023
#' Archived SCTransform integration: <https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette>, compiled Jan 2024

# new version