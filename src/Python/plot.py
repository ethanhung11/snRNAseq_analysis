import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from anndata import AnnData
import scanpy as sc
import glasbey

def empty_axs(axs:np.ndarray):
    for rs in axs[:]:
        for ax in rs[:]:
            ax.remove()
    return


def order_obs(adata, col, order):
    adata.obs[col] = pd.Categorical(adata.obs[col], categories=order, ordered=True)
    return

def color_gen(groups):
    return glasbey.create_palette(palette_size=len(groups.unique()))


def check_integration(adata:AnnData, category:str, f, embeddings = ['X_umap', 'LocalMAP'], 
                      nrow:int=None, ncol:int=None, mini=False,
                      ):
    print(f"Category {category} has {len(adata.obs[category].unique())} groups!")
    
    sf = f.subfigures(1,len(embeddings))
    int_colors = color_gen(adata.obs[category])

    for e,obsm in enumerate(embeddings):
        if mini is True:
            ax = sf[e].subplots(1,1)
            sc.pl.embedding(adata, basis=obsm, color=category,
                            ax=ax, show=False,
                            palette=int_colors)
                            # legend_loc='none' if e<len(embeddings)-1 else "right margin")

        else:
            assert nrow is not None and ncol is not None
            axs = sf[e].subplots(nrow*2, ncol)
            gs = axs[0,0].get_gridspec()
            empty_axs(axs)
            ax = sf[e].add_subplot(gs[:nrow,:])
            sc.pl.embedding(adata, basis=obsm, color=category,
                            ax=ax, show=False,
                            palette=int_colors)

            for n,group in enumerate(adata.obs[category].unique()):
                ax = sf[e].add_subplot(gs[nrow + n//nrow, n%nrow])
                sc.pl.embedding(adata[adata.obs[category] == group],
                                basis=obsm, color=category,
                                ax=ax, show=False,
                                palette=[int_colors[n]],
                                legend_loc='none')
                ax.set_title(group)

    return


def cluster_violinplot(adata, markers, group, f, layer='normalized'):
    axs = f.subplots(len(markers),1)
    for n,m in enumerate(markers):
        sc.pl.violin(adata, m, groupby=group, use_raw=False, layer=layer,
                    show=False, ax=axs[n],)
        if n<len(markers)-1:
            axs[n].set_xlabel('')
            axs[n].set_xticklabels(['']*len(axs[n].get_xticklabels()))
        axs[n].set_ylabel(axs[n].get_ylabel(),size=12)
    return


def cluster_stackedbarplot(adata:AnnData, col1:str, col2:str, pct:bool = False, colors:list = None, ax = None):
    if colors is None:
        colors = color_gen(adata.obs[col2])
    if ax is None:
       f, ax = plt.subplots(1, 2, figsiz=(5,5), layout="constrained")

    crosstab_counts = pd.crosstab(adata.obs[col1], adata.obs[col2])
    crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=1), axis=0) * 100

    if pct is False:
        # Counts
        crosstab_counts.plot(kind='bar', stacked=True, ax=ax, color=colors)
        ax.set_title('Cluster Composition by Group (Counts)', fontsize=14, fontweight='bold')
        ax.set_xlabel(col1, fontsize=12)
        ax.set_ylabel('Cell Count', fontsize=12)
        ax.legend(title=col2, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=45)
        
    else:
        # Percentages
        crosstab_pct.plot(kind='bar', stacked=True, ax=ax, color=colors)
        ax.set_title('Cluster Composition by Group (Percentages)', fontsize=14, fontweight='bold')
        ax.set_xlabel(col1, fontsize=12)
        ax.set_ylabel('Percentage (%)', fontsize=12)
        ax.legend(title=col2, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=45)
        ax.set_ylim(0, 100)
    
    return

