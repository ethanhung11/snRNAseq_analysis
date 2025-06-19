from typing import Iterable
from anndata import AnnData

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import glasbey

def empty_axs(axs:np.ndarray):
    for rs in axs[:]:
        for ax in rs[:]:
            ax.remove()
    return

def order_obs(adata:AnnData, col:str, order:Iterable[str]):
    adata.obs[col] = pd.Categorical(adata.obs[col], categories=order, ordered=True)
    return

def color_gen(groups:pd.Series | list | np.array, index=None):
    cs = glasbey.create_palette(palette_size=len(groups.unique()))
    if index is not None:
        return pd.Series(cs,index=index)
    else:
        return pd.Series(cs,index=groups.unique())


def check_integration(adata:AnnData, category:str, f, embeddings:Iterable[str] = ['X_umap', 'LocalMAP'], 
                      nrow:int=None, ncol:int=None, mini=False,
                      ):
    print(f"Category {category} has {len(adata.obs[category].unique())} groups!")
    
    sf = f.subfigures(1,len(embeddings))
    int_colors = color_gen(adata.obs[category],adata.obs[category].cat.categories)

    for e,obsm in enumerate(embeddings):
        if mini is True:
            ax = sf[e].subplots(1,1)
            sc.pl.embedding(adata, basis=obsm, color=category,
                            ax=ax, show=False,
                            palette=int_colors.to_list())
                            # legend_loc='none' if e<len(embeddings)-1 else "right margin")

        else:
            assert nrow is not None and ncol is not None
            axs = sf[e].subplots(nrow*2, ncol)
            gs = axs[0,0].get_gridspec()
            empty_axs(axs)
            ax = sf[e].add_subplot(gs[:nrow,:])
            sc.pl.embedding(adata, basis=obsm, color=category,
                            ax=ax, show=False,
                            palette=int_colors.to_list())

            for n,group in enumerate(adata.obs[category].unique()):
                ax = sf[e].add_subplot(gs[nrow + n//nrow, n%nrow])
                sc.pl.embedding(adata[adata.obs[category] == group],
                                basis=obsm, color=category,
                                ax=ax, show=False,
                                palette=[int_colors[group]],
                                legend_loc='none')
                ax.set_title(group)

    return


def cluster_violinplot(adata:AnnData, markers:Iterable[str], group:str, f, layer:str='normalized'):
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
        crosstab_pct = crosstab_counts.div(crosstab_counts.sum(axis=1), axis=0) * 100

        crosstab_pct.plot(kind='bar', stacked=True, ax=ax, color=colors)
        ax.set_title('Cluster Composition by Group (Percentages)', fontsize=14, fontweight='bold')
        ax.set_xlabel(col1, fontsize=12)
        ax.set_ylabel('Percentage (%)', fontsize=12)
        ax.legend(title=col2, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.tick_params(axis='x', rotation=45)
        ax.set_ylim(0, 100)
    
    return

def cluster_c2c(adata:AnnData, key:str, pval:float=0.05, top_n:int=20, 
                  sources:Iterable[str]=None, targets:Iterable[str]=None, figsize=(15, 8)):
    """Lightweight cell-cell communication dotplot with controlled dot sizes."""
    
    df = adata.uns[key]
    df = df[df['specificity_rank'] < pval]
    
    if sources: df = df[df['source'].isin(sources)]
    if targets: df = df[df['target'].isin(targets)]
    
    df = df.nsmallest(top_n, 'magnitude_rank')
    
    # Use complex columns if available, fallback to simple
    lig_col = 'ligand_complex' if 'ligand_complex' in df.columns else 'ligand'
    rec_col = 'receptor_complex' if 'receptor_complex' in df.columns else 'receptor'
    df['lr'] = df[lig_col] + ' -> ' + df[rec_col]
    
    # Efficient size mapping: 5 quantile-based bins, min size 10
    mag_vals = -np.log10(np.maximum(df['magnitude_rank'], 1e-10))
    size_bins = np.array([10, 50, 100, 150, 200])
    # Use 1 as minimum, then map higher values to larger sizes
    size_norm = plt.Normalize(vmin=1, vmax=mag_vals.max())
    df['dot_size'] = size_bins[np.searchsorted(np.percentile(mag_vals, [20, 40, 60, 80]), np.maximum(mag_vals, 1))]
    
    src_list = sorted(df['source'].unique())
    tgt_list = sorted(df['target'].unique())
    lr_list = df['lr'].drop_duplicates().tolist()
    
    fig, axes = plt.subplots(1, len(src_list), figsize=figsize, sharey=True)
    if len(src_list) == 1: axes = [axes]
    
    # Pre-compute color normalization
    spec_vals = -np.log10(np.maximum(df['specificity_rank'], 1e-10))
    norm = plt.Normalize(vmin=1, vmax=spec_vals.max())
    
    for i, src in enumerate(src_list):
        ax = axes[i]
        data = df[df['source'] == src]
        
        if len(data):
            x = [tgt_list.index(t) for t in data['target']]
            y = [lr_list.index(lr) for lr in data['lr']]
            colors = plt.cm.Greens(norm(-np.log10(np.maximum(data['specificity_rank'], 1e-10))))
            ax.scatter(x, y, c=colors, s=data['dot_size'], alpha=0.8, edgecolors=colors, lw=1)
        
        ax.set(xlim=(-0.5, len(tgt_list)-0.5), ylim=(-0.5, len(lr_list)-0.5),
               xticks=range(len(tgt_list)), title=src)
        ax.set_xticklabels(tgt_list, rotation=45, ha='right')
        ax.grid(alpha=0.3)
    
    axes[0].set(yticks=range(len(lr_list)), ylabel='Ligand -> Receptor')
    axes[0].set_yticklabels(lr_list)
    fig.suptitle('Source')
    fig.text(0.5, 0.02, 'Target', ha='center')
    
    plt.tight_layout()
    fig.subplots_adjust(right=0.75)
    
    # Colorbar
    sm = plt.cm.ScalarMappable(cmap='Greens', norm=norm)
    plt.colorbar(sm, ax=axes, label='-log10(specificity p-val)')
    
    # Create proper size legend using evenly spaced integer values
    if len(mag_vals) > 0:  # Only create legend if we have data
        # Create evenly spaced integer thresholds from 1 to max
        max_val = int(np.ceil(mag_vals.max()))
        if max_val <= 1:
            size_labels = ['1']
        else:
            # Create 5 evenly spaced integers from 1 to max
            thresholds = np.linspace(1, max_val, 5).astype(int)
            # Remove duplicates while preserving order
            thresholds = np.array(sorted(set(thresholds)))
            size_labels = [str(t) for t in thresholds]
        
        # Create legend elements for matplotlib
        from matplotlib.lines import Line2D
        legend_elements = []
        # Use only the number of sizes we actually have
        for i, label in enumerate(size_labels):
            if i < len(size_bins):
                size = size_bins[i]
                legend_elements.append(Line2D([0], [0], marker='o', color='w', 
                                            markerfacecolor='gray', markersize=np.sqrt(size/10),
                                            alpha=0.7, markeredgecolor='k', label=label))
        
        # Add legend to the figure
        fig.legend(handles=legend_elements, title='Magnitude\n(-log10 p-val)', 
                  loc='center right', bbox_to_anchor=(0.83, 0.5), frameon=True, 
                  fancybox=True, shadow=True)
    
    return