import os
import scanpy as sc
import rpy2.robjects as ro
from .R import R_preload, get_converter
from scanpy import AnnData


def sc_get_data(directory: str, name: str, filetype: str):
    if filetype == ".h5":
        adata = sc.read_10x_h5(directory + filetype)
    elif filetype == ".h5ad":
        adata = sc.read_h5ad(directory + filetype)
    elif filetype == ".mtx":
        adata = sc.read_10x_mtx(directory)
    else:
        raise ValueError(
            f"{filetype} is not a valid filetype. Must be h5, h5ad, or mtx."
        )

    adata.obs["Identifier"] = name
    adata.obs.index = adata.obs.index + "_" + name
    adata.layers["counts"] = adata.X.copy()

    return adata


def get_data(
    data_dir: str,
    experiment_dir: str,
    filename: str,
    filetype: str,
    convertR: bool = False,
):
    if convertR is True:
        basedir = os.path.join(data_dir, experiment_dir)
        with ro.conversion.localconverter(get_converter()):
            ro.globalenv["datadir"] = basedir
            ro.globalenv["filename"] = filename
            ro.globalenv["filetype"] = filetype

            print(basedir, filename, filetype)

            R_preload()
            ro.r(
                """
            datadir <- here(datadir)
            experiments <- list.files(datadir)
            adatas <- c()
            sces <- c()
            for (exp in experiments) {
                if (filetype==".h5") {
                    filedir <- paste(datadir, exp, paste0(filename,filetype), sep="/")
                    print(filedir)
                    seurat_data <- Read10X_h5(filedir, use.names = T, unique.features=T)
                } else if (filetype==".mtx") {
                    filedir <- paste(datadir, exp, filename, sep="/")
                    print(filedir)
                    seurat_data <- Read10X(filedir, unique.features=T)
                }
                
                seurat_obj <- CreateSeuratObject(counts = seurat_data, project = exp)
                sce <- convert_seurat_to_sce(seurat_obj)
                sces <- c(sces,sce)
            }
            """
            )
            adatas = list(ro.globalenv["sces"])

        for ad in adatas:
            ad.obs["Identifier"] = ad.obs["orig.ident"]
            ad.obs.index = ad.obs.index + "_" + ad.obs["orig.ident"].astype(str)
            ad.layers["counts"] = ad.X.copy()
            del ad.uns

    else:
        basedir = os.path.join(data_dir, experiment_dir)
        samples = os.listdir(basedir)
        filt_files = [os.path.join(basedir, s, filename) for s in samples]
        adatas = [
            sc_get_data(file, samples[n], filetype) for n, file in enumerate(filt_files)
        ]

    return adatas


def create_cloupe(adata: AnnData):
    with ro.conversion.localconverter(get_converter()):
        R_preload()
        ro.globalenv["sce"] = adata
        ro.r("""
        library("loupeR")
        clust <- as.list(colData(sce)) %>% lapply(as.factor)
        proj <- lapply(as.list(reducedDims(sce)), function(df) df %>% select(1:2))

        create_loupe(
            assay(sce, "X"),
            clusters = clust, 
            projections = proj,
            output_name = "output",
        )
        """)