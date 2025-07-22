from typing import Iterable
from anndata import AnnData

import os
import scanpy as sc
import liana as li
import decoupler as dc
import pandas as pd
from .ccc import liana_mouse_resource


def cell2cell_interactions(
    adata: AnnData,
    cell_group="cell_type",
    key="ccc",
    resource_opts: Iterable[str] = [
        "baccin2019",
        "cellcall",
        "cellchatdb",
        "cellinker",
        "cellphonedbv5",
        "cellphonedb",
        "celltalkdb",
        "connectomedb2020",
        "consensus",
        "embrace",
        "guide2pharma",
        "hpmr",
        "icellnet",
        "italk",
        "kirouac2010",
        "lrdb",
        "ramilowski2015",
        "mouseconsensus",
    ],
    methods: Iterable = None,
    cores: int = None,
):
    if cores is None:
        cores = os.cpu_count()

    # pull interaction databases
    ccc_db = liana_mouse_resource(resource_opts)

    # run all methods
    if methods is None:
        rank_aggregate_custom = li.method.rank_aggregate
    else:
        rank_aggregate_custom = li.mt.AggregateClass(
            li.mt.aggregate_meta, methods=methods
        )

    rank_aggregate_custom(
        adata,
        groupby=cell_group,
        layer="normalized",
        use_raw=False,
        key_added=key,
        de_method="wilcoxon",
        return_all_lrs=True,
        verbose=True,
        n_jobs=cores,
        resource=ccc_db[["ligand", "receptor"]],
    )

    # save reference database info
    adata.uns[key] = adata.uns[key].merge(
        ccc_db[["ligand", "receptor", "db_sources"]],
        left_on=["ligand_complex", "receptor_complex"],
        right_on=["ligand", "receptor"],
        how="left",
    )

    # filter for quality interactions
    ccc_filters = []

    # Cell Specificity filters
    if "cellphone_pvals" in adata.uns[key].columns:  # CellphoneDB
        ccc_filters.append(adata.uns[key]["cellphone_pvals"] <= 0.05)
    if "gmean_pvals" in adata.uns[key].columns:  # CellphoneDB V2
        ccc_filters.append(adata.uns[key]["gmean_pvals"] <= 0.05)
    if "cellchat_pvals" in adata.uns[key].columns:  # CellChat
        ccc_filters.append(adata.uns[key]["cellchat_pvals"] <= 0.05)
    if "lr_logfc" in adata.uns[key].columns:  # log2FC
        ccc_filters.append(adata.uns[key]["lr_logfc"] > 0)
    if "scaled_weight" in adata.uns[key].columns:  # Connectome
        ccc_filters.append(
            adata.uns[key]["scaled_weight"]
            > adata.uns[key]["scaled_weight"].quantile(0.95)
        )
    if "spec_weight" in adata.uns[key].columns:  # NATMI
        ccc_filters.append(
            adata.uns[key]["spec_weight"] > adata.uns[key]["spec_weight"].quantile(0.95)
        )

    # Magnitue filters
    if "lr_probs" in adata.uns[key].columns:  # CellChat
        ccc_filters.append(adata.uns[key]["lr_probs"] <= 0.05)
    if "lrscore" in adata.uns[key].columns:  # SingleCellSignalR
        ccc_filters.append(adata.uns[key]["lrscore"] > 0.6)
    if "expr_prod" in adata.uns[key].columns:  # NATMI/Connectome
        ccc_filters.append(
            adata.uns[key]["expr_prod"] > adata.uns[key]["expr_prod"].quantile(0.95)
        )
    if "magnitude_rank" in adata.uns[key].columns:  # Liana (aggregated score)
        ccc_filters.append(adata.uns[key]["magnitude_rank"] <= 0.05)

    df = pd.concat(ccc_filters, axis=1)
    ccc_filter_all = df.all(axis=1)
    adata.uns[key + "_filtered"] = adata.uns[key][ccc_filter_all].reset_index(drop=True)

    return adata


def GO_Enrich(
    adata,
    groupby,
    key,
    sources=["GO:MF", "GO:CC", "GO:BP", "KEGG", "REAC"],
    pval_cutoff=1e-4,
    log2fc_min=2,
):
    # see here for other gprofilier args: https://biit.cs.ut.ee/gprofiler/page/apis

    GO_enrichments = {}
    for src in sources:
        GO_enrichments[src] = {}

    for category in adata.obs[groupby].unique():
        df = sc.get.rank_genes_groups_df(
            adata,
            group=category,
            key=key,
            pval_cutoff=pval_cutoff,
            log2fc_min=log2fc_min,
        )
        for src in sources:
            GO_enrichments[src][category] = sc.queries.enrich(
                df.names.to_list(),
                org="mmusculus",
                gprofiler_kwargs={"sources": [src], "all_results": True},
            )

    return GO_enrichments


def GSEA_decoupler(
    adata: AnnData,
    name: str,
    type: str = "ULM",
    geneset_dir: str = None,
    geneset=None,
):
    if geneset is None:
        assert geneset_dir is not None
        geneset = dc.pp.read_gmt(geneset_dir)

    if type == "ULM":
        dc.mt.ulm(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_ulm"] = adata.obsm["score_ulm"]
        adata.obsm[f"{name}_padj_ulm"] = adata.obsm["padj_ulm"]
        del adata.obsm["score_ulm"], adata.obsm["padj_ulm"]

    elif type == "GSEA":
        dc.mt.gsea(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_gsea"] = adata.obsm["score_gsea"]
        adata.obsm[f"{name}_padj_gsea"] = adata.obsm["padj_gsea"]
        del adata.obsm["score_gsea"], adata.obsm["padj_gsea"]

    elif type == "AUCell":
        dc.mt.aucell(data=adata, net=geneset, layer="normalized", verbose=True)
        adata.obsm[f"{name}_score_aucell"] = adata.obsm["score_aucell"]
        del adata.obsm["score_aucell"]

    return adata
