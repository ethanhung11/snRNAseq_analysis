from typing import Iterable
from anndata import AnnData

import os
import liana as li
import pandas as pd
from ccc import liana_mouse_resource


def cell2cell_interactions(adata:AnnData,
                           methods:Iterable=None,
                           resource_opts:Iterable[str]=['baccin2019', 'cellcall', 'cellchatdb', 'cellinker', 'cellphonedbv5', 'cellphonedb', 'celltalkdb', 'connectomedb2020', 'consensus', 'embrace', 'guide2pharma', 'hpmr', 'icellnet', 'italk', 'kirouac2010', 'lrdb', 'mouseconsensus', 'ramilowski2015'],
                           cores:int=None,
):
    if cores is None:
        cores = os.cpu_count()
        
    # pull interaction databases
    ccc_db = liana_mouse_resource(resource_opts)

    # run all methods
    if methods is None:
        rank_aggregate_custom = li.method.rank_aggregate
    else:
        rank_aggregate_custom = li.mt.AggregateClass(li.mt.aggregate_meta, methods=methods)
    rank_aggregate_custom(
        adata, groupby="cell_type", layer='normalized', use_raw=False, key_added="ccc_aggregate",
        return_all_lrs=True, verbose=True, n_jobs=cores, resource=ccc_db[["ligand", "receptor"]]
    )

    # save reference database info
    adata.uns["ccc_aggregate"] = adata.uns["ccc_aggregate"].merge(ccc_db[['ligand', 'receptor', 'db_sources']],
                                                                  left_on=["ligand_complex", "receptor_complex"],
                                                                  right_on=['ligand','receptor'],
                                                                  how='left')
    
    # filter for quality interactions
    ccc_filters = []

    # Cell Specificity filters
    if "cellphone_pvals" in adata.uns["ccc_aggregate"].columns: # CellphoneDB
        ccc_filters.append(adata.uns["ccc_aggregate"]["cellphone_pvals"] <= 0.05) 
    if "gmean_pvals" in adata.uns["ccc_aggregate"].columns: # CellphoneDB V2
        ccc_filters.append(adata.uns["ccc_aggregate"]["gmean_pvals"] <= 0.05)
    if "cellchat_pvals" in adata.uns["ccc_aggregate"].columns: # CellChat
        ccc_filters.append(adata.uns["ccc_aggregate"]["cellchat_pvals"] <= 0.05)
    if "lr_logfc" in adata.uns["ccc_aggregate"].columns: # log2FC
        ccc_filters.append(adata.uns["ccc_aggregate"]["lr_logfc"] > 0)
    if "scaled_weight" in adata.uns["ccc_aggregate"].columns: # Connectome
        ccc_filters.append(adata.uns["ccc_aggregate"]["scaled_weight"] > adata.uns["ccc_aggregate"]["scaled_weight"].quantile(.95))
    if "spec_weight" in adata.uns["ccc_aggregate"].columns: # NATMI
        ccc_filters.append(adata.uns["ccc_aggregate"]["spec_weight"] > adata.uns["ccc_aggregate"]["spec_weight"].quantile(.95))

    # Magnitue filters
    if "lr_probs" in adata.uns["ccc_aggregate"].columns: # CellChat
        ccc_filters.append(adata.uns["ccc_aggregate"]["lr_probs"] <= 0.05)
    if "lrscore" in adata.uns["ccc_aggregate"].columns: # SingleCellSignalR
        ccc_filters.append(adata.uns["ccc_aggregate"]["lrscore"] > 0.6)
    if "expr_prod" in adata.uns["ccc_aggregate"].columns: # NATMI/Connectome
        ccc_filters.append(adata.uns["ccc_aggregate"]["expr_prod"] > adata.uns["ccc_aggregate"]["expr_prod"].quantile(.95))
    if "magnitude_rank" in adata.uns["ccc_aggregate"].columns: # Liana (aggregated score)
        ccc_filters.append(adata.uns["ccc_aggregate"]["magnitude_rank"] <= 0.05)

    df = pd.concat(ccc_filters, axis=1)
    ccc_filter_all = df.all(axis=1)
    adata.uns["ccc_aggregate_filtered"] = adata.uns["ccc_aggregate"][ccc_filter_all].copy()
    
    return adata