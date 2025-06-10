import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects.packages import importr


def install_package(packages: list) -> bool:
    utils = importr("utils")
    utils.chooseCRANmirror(ind=1)
    for p in packages:
        utils.install_packages(p)
    return True


def get_converter():
    return ro.default_converter + ro.pandas2ri.converter + anndata2ri.converter

def R_preload() -> bool:
    ro.r('''
    library(scry)
    library(scater)
    library(scDblFinder)
    library(BiocParallel)
    library(Seurat)
    library(dplyr)
    library(DoubletFinder)
    options(Seurat.object.assay.version = "v3")
    ''')
