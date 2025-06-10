#!/usr/bin/env Rscript

library(glue)
library(here)
library(dplyr)
library(parallel)
library(ggplot2)

library(hdf5r)
library(patchwork)
library(Seurat)
library(DoubletFinder)
options(Seurat.object.assay.version = 'v3')

retrieveData <- function(directory,inputtype) {
  sample <- strsplit(directory,'/')[[1]] %>% nth(-2)
  if (inputtype=="folder"){
    dataset <- here(directory) %>%
      Read10X() %>%
      CreateSeuratObject(min.cells=10, project = sample)
  } else if (inputtype=="h5"){
    dataset <- here(directory) %>%
      Read10X_h5() %>%
      CreateSeuratObject(min.cells=10, project = sample)
  } else {
    stop("Not a valid data type!")
  }
  
  return(dataset)
}

QCFilter <- function(SCdata,visualize=F,verbose=T) {
  #' <https://satijalab.org/seurat/articles/sctransform_vignette>

  # QC scores
  # ----------------------------------------------------------------------------
  if (verbose) {printData(SCdata)}
  SCdata[["percent.mt"]] <- PercentageFeatureSet(SCdata, pattern = "^Mt")/100
  SCdata[["percent.rb"]] <- PercentageFeatureSet(SCdata, pattern = "^Rb")/100
  SCdata[["CountPerFeature"]] <- log10(SCdata$nFeature_RNA) / log10(SCdata$nCount_RNA)
  # Plotting QC
  if (visualize) {
    plots <- VlnPlot(SCdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2) +
      FeatureScatter(SCdata, feature1 = "nCount_RNA", feature2 = "percent.mt") +
      FeatureScatter(SCdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  } else {
    plots <- c()
  }
  # each sample was further filtered...based on the number of UMIs (≥500) and genes (≥250), mitochondrial ratio (<15%), and the number of genes detected per UMI (>0.8)
  if (verbose) {
    cat("Reads Filter: ", sum(SCdata$nCount_RNA<=500)%>%as.character(),
        "\nGene Filter: ", sum(SCdata$nFeature_RNA<=250)%>%as.character(),
        "\nMT Filter: ", sum(SCdata$percent.mt>=0.15)%>%as.character(),
        "\nRB Filter: ", sum(SCdata$percent.rb>=0.15)%>%as.character(),
        "\nReads/Gene Filter: ", sum(SCdata$CountPerFeature<=0.8)%>%as.character(),
        "\n")
  }
  SCdata <- subset(SCdata, subset = nCount_RNA>500 & nFeature_RNA>250 & percent.mt<0.15 & percent.rb<0.15 & CountPerFeature>0.8)
  if (verbose) printData(SCdata)
  
  # filter out mitochondrial & ribosomal genes
  # ----------------------------------------------------------------------------
  SCdata <- SCdata[!grepl("^Mt", rownames(SCdata)), ]
  SCdata <- SCdata[!grepl("^Rp", rownames(SCdata)), ]
  if (verbose) printData(SCdata)
  
  # SCTransform for Normalization
  # ----------------------------------------------------------------------------
  nPCs <- 50
  SCdata <- SCTransform(SCdata, verbose=verbose) %>%
    RunPCA(npcs = nPCs, verbose=verbose) %>%
    RunUMAP(reduction = "pca", dims = 1:nPCs, verbose=verbose) %>%
    FindNeighbors(reduction = "pca", dims = 1:nPCs, verbose=verbose) %>%
    FindClusters(resolution = 0.7, verbose=verbose)
  
  if (verbose) printData(SCdata)
  
  ## Find relevant PCS
  stdv <- SCdata[["pca"]]@stdev
  percent_stdv <- (stdv/sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  if (visualize) {
    df <- data.frame(r2=percent_stdv,comp=c(1:length(percent_stdv)))
    plots <- plots + ggplot(df, aes(x=comp, y=r2)) +
      geom_bar(stat = "identity")
  }
  nPCs <- which(cumulative > 90 & percent_stdv < 1)[1]
  if (verbose) cat("PCS:",nPCs,'\n')
  
  if (verbose) {print("QCFilter Complete!")}
  return(list(data=SCdata,PCs=nPCs,plots=plots))
}

DoubletFilter <- function(input,visualize=F,verbose=T,multipletRate=0.075) {
  #'https://biostatsquid.com/doubletfinder-tutorial/
  #'https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #'https://rpubs.com/kenneditodd/doublet_finder_example
  
  SCdata <- input$data
  nPCs <- input$PCs
  plots <- input$plots
  ## 1. pK Identification (no ground-truth)
  #----------------------------------------------------------------------------
  sweep.res.list <- paramSweep(SCdata, PCs = 1:nPCs, sct=TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  optimal.pk <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>%
    select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))
  
  # 2. Homotypic Doublet Proportion Estimate
  # see https://kb.10xgenomics.com/hc/en-us/articles/360054599512-What-is-the-cell-multiplet-rate-when-using-the-3-CellPlex-Kit-for-Cell-Multiplexing
  #----------------------------------------------------------------------------
  annotations <- SCdata@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(multipletRate*nrow(SCdata@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## 3. Run DoubletFinder
  #----------------------------------------------------------------------------
  SCdata <- doubletFinder(SCdata, PCs = 1:nPCs, pK = optimal.pk, nExp = nExp_poi, sct=TRUE)
  colnames(SCdata@meta.data)[grepl('DF.classifications.*', colnames(SCdata@meta.data))] <- "doublet_finder"
  if (visualize){
    plots <- plots + VlnPlot(seu, group.by = 'SampleID', split.by = "doublet_finder",
            features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_ribo", "percent_hb"), 
            ncol = 3, pt.size = 0) + theme(legend.position = 'right')
  }
  SCdata <- subset(SCdata, subset = doublet_finder=='Singlet')
  if (verbose) {printData(SCdata)}
  if (verbose) {print("DoubletFinder Complete!")}
  return(list(data=SCdata,plots=plots))
}

preprocessData <- function(dataset,multipletRate,visualize,verbose){
  processed <- dataset %>%
    QCFilter(visualize=visualize,verbose=verbose) %>%
    DoubletFilter(visualize=visualize,verbose=verbose,multipletRate=multipletRate)
    
  plots <- processed$plots
  data <- processed$data %>% SCTransform(min_cells=10,vars.to.regress = c('percent.mt'))
  return(data)
}

#' CellCycleFilter <- function(SCdata) {
#'   #'https://satijalab.org/seurat/articles/cell_cycle_vignette
#'   SCdata <- CellCycleScoring(object = SCdata, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
#'   datasets <- lapply(datasets,SCTransform,min_cells=10,vars.to.regress = c('percent.mt',g2m.features,cc.genes$s.genes))
#' }

printData <- function(seurat.obj) {
  print(glue("{c('cells','genes')}: {dim(seurat.obj)}"))
}


#'  V5 Version -- Not sure how to get it to work, but not necessary for replicating old analysis
# integrateData <- function(seurat.list, name="integrated", visualize=F, verbose=T) {
#   #' https://satijalab.org/seurat/articles/integration_introduction#perform-integration-with-sctransform-normalized-datasets
#   #' https://github.com/satijalab/seurat/issues/9693
#   data_labels <- seurat.list %>% lapply(attr,which="project.name") %>% unlist()
  
#   merged_data <- merge(seurat.list[[1]],
#                        y = seurat.list[-1],
#                        add.cell.ids = data_labels, 
#                        project = name)
  
#   integrated_data <- IntegrateLayers(object = merged_data,
#                                      normalization.method = "SCT",
#                                      method = CCAIntegration,
#                                      features=3000,
#                                      verbose=verbose)
  
#   return(integrated_data)
# }

integrateData <- function(data.list) {
  features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = features)
  anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",
      anchor.features = features)
  integrated_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  return(integrated_data)
}


if (!interactive()) {
  library(optparse)
  print("imported all libraries")
  
  # Define command line options
  opt_list <- c(make_option(c("-i", "--input"), type="character", action="store", dest="directories", help="Directory to process [required]"),
                make_option(c("-n", "--name"), type="character", action="store", help="Processing Run Name [required]"),
                make_option(c("-t", "--inputtype"), type="character", action="store", default="h5", help="Processing Run Name [required]"),
                make_option(c("--multipletRate"), type="double", action="store", default=0.075, help="Estimated doublet rate, defaults to 0.075"),
                make_option(c("--cores"), type="integer", action="store", default=10, help="Number of cores to use, defaults to 10"),
                make_option(c("--verbose"), type="logical", action="store", default=TRUE, help="verbosity setting, defaults to TRUE"),
                make_option(c("--visualize"), type="logical", action="store", default=FALSE, help="show visualizations in interactive mode, defaults to FALSE"))
  
  opt <- parse_args(OptionParser(option_list=opt_list))
  
  if (opt$directories == "default") {
    opt$directories <- paste("data","preprocessed",here("data","preprocessed") %>% list.files(),sep='/')
    opt$name <- "testing"
    opt$inputtype <- "folder"
    opt$multipletRate <- 0.075
    opt$cores <- 10
    opt$verbose <- TRUE
    opt$visualize <- FALSE
  }
  
  opt$directories <- unlist(strsplit(opt$directories, ","))
  
  print(opt)
  
  
  if (opt$verbose)
    print("begin accessing datasets...")
  system.time(
      datasets <- mclapply(opt$directories,retrieveData,inputtype=opt$inputtype,mc.cores=opt$cores)
  )
  
  if (opt$verbose)
    print(datasets)
    print("begin preprocessing datasets...")
  system.time(
    processed_datasets <- mclapply(datasets,preprocessData,mc.cores=opt$cores,
                                   multipletRate=opt$multipletRate,
                                   visualize=opt$visualize,
                                   verbose=opt$verbose)
  )
  
  savename <- here("data","processed",glue("{opt$name}-processed.Rds"))
  if (opt$verbose)
    print("save preprocessed datasets...")
    print(savename)
  saveRDS(processed_datasets,file=savename)
}
