from typing import Literal, List
from anndata import AnnData

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scvi
from doubletdetection import BoostClassifier
from pacmap import LocalMAP
import rpy2.robjects as ro
from .R import get_converter


def Filter_QC(
    adata: AnnData,
    GenePerCell: int = 250,
    CountPerCell: int = 500,
    CellPerGene: int = 10,
    mt_perc: float = 5,
    ribo_perc: float = 5,
    verbose: bool = False,
):

    start = adata.shape[0]
    sc.pp.filter_cells(adata, min_genes=GenePerCell)
    sc.pp.filter_cells(adata, min_counts=CountPerCell)
    sc.pp.filter_genes(adata, min_cells=CellPerGene)
    adata.var_names_make_unique()
    if verbose is True:
        print(f"Cells removed by gene/count filters: {start - adata.shape[0]}")

    adata.var["mt"] = adata.var_names.str.startswith("Mt")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=[20], log1p=True
    )

    remove = [
        "total_counts_mt",
        "log1p_total_counts_mt",
        "total_counts_ribo",
        "log1p_total_counts_ribo",
    ]
    adata.obs = adata.obs[[x for x in adata.obs.columns if x not in remove]]

    if verbose is True:
        print(
            f"Cells with >{mt_perc}% mt genes: {np.sum(adata.obs["pct_counts_mt"] >= mt_perc)}"
        )
        print(
            f"Cells with >{ribo_perc}% ribo genes: {np.sum(adata.obs["pct_counts_ribo"] >= ribo_perc)}"
        )

    adata = adata[adata.obs["pct_counts_mt"] <= mt_perc, :]
    adata = adata[adata.obs["pct_counts_ribo"] <= ribo_perc, :]

    return adata


def Filter_Doublet(
    adata: AnnData,
    method: Literal[
        "doubletdetection", "scrublet", "scDblFinder", "DoubletFinder", "SOLO"
    ] = "doubletdetection",
    remove: bool = True,
    multipletRate: float = 0.075,
):

    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    # https://github.com/mousepixels/sanbomics_scripts/blob/main/sc2024/preprocessing.ipynb
    if method == "doubletdetection":
        clf = BoostClassifier(
            n_iters=10,
            clustering_algorithm="louvain",
            standard_scaling=True,
            pseudocount=0.1,
            n_jobs=10,
        )

        clf.fit(adata.X)
        doublets = clf.predict(p_thresh=1e-3, voter_thresh=0.5)
        doublet_score = clf.doublet_score()

        adata.obs["predicted_doublet"] = doublets.astype(bool)
        adata.obs["doublet_score"] = doublet_score
    # https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html#doublet-detection
    elif method == "scrublet":
        sc.pp.scrublet(adata)
    elif method == "scDblFinder":
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.r(
                """
            library(scater)
            library(scDblFinder)
            library(BiocParallel)
                 
            sce <- scDblFinder(sce)
            doublet_score <- sce$scDblFinder.score
            doublet_class <- sce$scDblFinder.class
            """
            )
            print(pd.Series(ro.globalenv["doublet_class"]))
            adata.obs["predicted_doublet"] = (
                pd.Series(ro.globalenv["doublet_class"])
                .replace({"singlet": 0, "doublet": 1})
                .to_list()
            )
            adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
            adata.obs["doublet_score"] = ro.globalenv["doublet_score"]
    elif method == "DoubletFinder":
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.globalenv["multipletRate"] = multipletRate
            ro.r(
                """
            library(scater)
            library(Seurat)
            library(dplyr)
            library(DoubletFinder)
            options(Seurat.object.assay.version = "v3")
                 
            data <- as.Seurat(sce,data = NULL)
            data <- RenameAssays(data,"originalexp","RNA")
            data <- NormalizeData(data) %>%
                FindVariableFeatures() %>%
                ScaleData()
                 
            data <- RunPCA(data, nfeatures.print = 10)
            stdv <- data[["pca"]]@stdev
            sum.stdv <- sum(data[["pca"]]@stdev)
            percent.stdv <- (stdv / sum.stdv) * 100
            cumulative <- cumsum(percent.stdv)
            co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
            co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                                percent.stdv[2:length(percent.stdv)]) > 0.1), 
                        decreasing = T)[1] + 1
            min.pc <- min(co1, co2)
            min.pc
            
            # Finish pre-processing 
            data <- RunUMAP(data, dims = 1:min.pc) %>%
                FindNeighbors(dims = 1:min.pc) %>%
                FindClusters(resolution = 0.8)

            ## 1. pK Identification (no ground-truth)
            #----------------------------------------------------------------------------
            sweep.res.list <- paramSweep(data, PCs = 1:min.pc)
            sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
            bcmvn <- find.pK(sweep.stats)
            bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
            optimal.pk <- bcmvn.max$pK
            optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
            
            # 2. Homotypic Doublet Proportion Estimate
            #----------------------------------------------------------------------------
            annotations <- data@meta.data$seurat_clusters
            homotypic.prop <- modelHomotypic(annotations)
            nExp_poi <- round(multipletRate * nrow(data@meta.data))
            nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
            
            ## 3. Run DoubletFinder
            #----------------------------------------------------------------------------
            data <- doubletFinder(data, PCs = 1:min.pc, pK = optimal.pk, nExp = nExp_poi.adj)
            colnames(data@meta.data)[grepl('DF.classifications.*', colnames(data@meta.data))] <- "doublet_finder"
            
            doublets <- data@meta.data$doublet_finder
            """
            )

            adata.obs["predicted_doublet"] = (
                pd.Series(ro.globalenv["doublets"])
                .replace({"Singlet": 0, "Doublet": 1})
                .to_list()
            )
            adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
            adata.obs["doublet_score"] = np.nan
    elif method == "SOLO":
        scvi.model.SCVI.setup_anndata(adata)
        vae = scvi.model.SCVI(adata)
        vae.train()
        solo = scvi.external.SOLO.from_scvi_model(vae)
        solo.train(early_stopping=True)
        adata.obs["prediction"] = solo.predict(soft=True)
        adata.obs["prediction"] = (
            solo.predict(soft=False)
            .prediction.map({"singlet": 0, "doublet": 1})
            .astype(bool)
        )
    else:
        raise ValueError("Not a valid doublet detection methodtype!")

    if remove is True:
        adata.uns["methods"]["doublet_remover"] = method
        adata.uns["doublets_removed"] = adata.obs["predicted_doublet"].sum()
        adata = adata[adata.obs["predicted_doublet"] == 0]

    else:
        if "doublet_remover" in adata.uns:
            adata.uns["doublet_remover"].append(method)
        else:
            adata.uns["doublet_remover"] = [method]
        adata.obs[f"predicted_doublet-{method}"] = adata.obs["predicted_doublet"]
        adata.obs[f"doublet_score-{method}"] = adata.obs["doublet_score"]
        del adata.obs["predicted_doublet"], adata.obs["doublet_score"]

    return adata


def Normalize(
    adata: AnnData | List[AnnData], kind: Literal["log1p", "mnn"] = "log1p", **kwargs
) -> AnnData:
    if type(adata) is list:
        adata = sc.concat(adata, join="outer")
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    print("Saving normalized counts in layer 'normalized'.")
    adata.layers["normalized"] = adata.layers["counts"].copy()

    if kind == "log1p":
        sc.pp.normalize_total(adata, layer="normalized")
        sc.pp.log1p(adata, layer="normalized")
        adata.uns["methods"]["normalization"] = "log1p"

    elif kind == "mnn":
        tmp = adata.copy()
        sce.pp.mnn(tmp)
        adata.layers["normalized"] = tmp.X
        adata.uns["methods"]["normalization"] = "mnn"

    else:
        raise ValueError("Not a valid normalization type!")

    return adata


def FindVariableGenes(
    adata: AnnData | List[AnnData],
    kind: Literal[
        "seurat", "seurat_v3", "seurat_v3_paper", "cell_ranger", "pearson", "deviant"
    ] = "seurat",
    batch_column: str = None,
    n_features=2000,
) -> AnnData:

    print("Finding HVGs...")

    if type(adata) is list:
        adata = sc.concat(adata, join="outer")
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    if kind == "seurat" or kind == "cell_ranger":
        sc.pp.highly_variable_genes(
            adata,
            batch_key=batch_column,
            flavor=kind,
            layer="normalized",
            n_top_genes=n_features,
        )

    elif kind == "seurat_v3" or kind == "seurat_v3_paper":
        sc.pp.highly_variable_genes(
            adata,
            batch_key=batch_column,
            flavor=kind,
            layer="counts",
            n_top_genes=n_features,
        )

    elif kind == "pearson":
        sc.experimental.pp.highly_variable_genes(
            adata,
            batch_key=batch_column,
            flavor=kind,
            layer="counts",
            n_top_genes=n_features,
        )

    elif kind == "deviant":
        with ro.conversion.localconverter(get_converter()):
            tmp = adata.copy()
            del tmp.uns
            ro.globalenv["sce"] = tmp
            ro.r(
                """
            library(scry)
            devianceFeatureSelection(sce, assay='X')
            """
            )
            result = ro.r("rowData(sce)$binomial_deviance").T

        idx = result.var["binomial_deviance"].argsort()[-n_features, :]
        mask = np.zeros_like(result, dtype=bool)
        mask[idx] = True
        adata.var["binomial_deviance"] = result.var["binomial_deviance"]
        adata.var["highly_variable"] = mask

    else:
        raise ValueError("Not a valid gene selection type!")

    adata.uns["methods"]["hvg"] = kind

    return adata


def Integrate(
    adata: AnnData,
    batch_column: str,
    use_var_genes: bool = True,
    kind: Literal["harmony", "scvi", "seurat"] = "harmony",
    **kwargs,
) -> AnnData:

    print(f"Integrating by Column {batch_column}: {adata.obs[batch_column].unique()}")
    if use_var_genes is True:
        assert "hvg" in adata.uns
    if "methods" not in adata.uns:
        adata.uns["methods"] = {}

    if kind == "harmony":
        sc.pp.pca(adata, use_highly_variable=use_var_genes, layer="normalized")
        sce.pp.harmony_integrate(
            adata,
            key=batch_column,
            max_iter_harmony=50,
            adjusted_basis="integrated",
            **kwargs,
        )

    elif kind == "scvi":
        if use_var_genes is True:
            input_adata = adata[:, adata.var["highly_variable"]].copy()
        else:
            input_adata = adata.copy()

        scvi.model.SCVI.setup_anndata(
            input_adata, layer="counts", batch_key=batch_column
        )
        model_scvi = scvi.model.SCVI(input_adata)
        model_scvi.train(early_stopping=True, **kwargs)
        adata.obsm["integrated"] = model_scvi.get_latent_representation()

    elif kind in "seurat":
        with ro.conversion.localconverter(get_converter()):
            assert use_var_genes is True

            subset = adata.copy()
            del subset.uns

            if "hvg" in adata.uns:
                ro.globalenv["sce"] = subset[:, subset.var["highly_variable"]]
                ro.r(
                    """
                    seurat <- as.Seurat(sce,data = NULL)
                    seurat <- RenameAssays(seurat,"originalexp","RNA")
                    """
                )
            else:
                ro.globalenv["sce"] = subset
                ro.r(
                    """
                    seurat <- as.Seurat(sce,data = NULL)
                    seurat <- RenameAssays(seurat,"originalexp","RNA")
                    seurat <- SCTransform(seurat, vars.to.regress = "pct_counts_mt")
                        """
                )

            ro.globalenv["batch_key"] = batch_column
            ro.r(
                f"""
            batch_list <- SplitObject(seurat, split.by = batch_key)
            anchors <- FindIntegrationAnchors(batch_list, anchor.features = rownames(seurat))
            integrated <- IntegrateData(anchors)
            integrated_expr <- GetAssayData(integrated)
            integrated_expr <- integrated_expr[rownames(seurat), colnames(seurat)]
            integrated_expr <- t(integrated_expr)
            """
            )
            subset = ro.globalenv["integrated_expr"]
            sc.pp.pca(subset, use_highly_variable=use_var_genes)

            adata.obsm["integrated"] = subset.obsm["X_pca"]
            adata.obsm["X_pca"] = subset.obsm["X_pca"]
            adata.uns["pca"] = subset.uns["pca"]

    elif kind == "scanvi":
        raise NotImplementedError

    adata.uns["methods"]["integration"] = kind
    return adata


def Visualize(adata, key=""):
    # UMAP
    print("Starting UMAP...")
    sc.pp.neighbors(adata, use_rep="integrated")
    sc.tl.umap(adata, key_added=f"UMAP{key}")

    # LocalMAP
    print("Starting LocalMAP...")
    lm = LocalMAP()
    adata.obsm[f"LocalMAP{key}"] = lm.fit_transform(adata.obsm["integrated"])

    return adata
