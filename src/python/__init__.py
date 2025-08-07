from .plot import (
    clear_uns,
    empty_axs,
    order_obs,
    color_gen,
    check_integration,
    checkDoublets,
    plot_violinplot,
    plot_cluster_violinplot,
    plot_cluster_stackedbarplot,
    plot_c2c,
    plot_gsea_dc,
    plot_go_enrichment,
)
from .preprocess import (
    Filter_QC,
    Filter_Doublet,
    Normalize,
    FindVariableGenes,
    Integrate,
    Visualize,
)
from .analysis import cell2cell_interactions, GSEA_decoupler, GO_Enrich
from .utils import GetSingleCellData, create_cloupe
from .R import get_converter, R_preload

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib as mpl


DEFAULT_CMAP = "Reds"
plt.rcParams['image.cmap'] = DEFAULT_CMAP
sns.set_palette(DEFAULT_CMAP) 