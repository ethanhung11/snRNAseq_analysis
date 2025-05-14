# snrnaseq-analysis

Replicating analysis by [So et al. 2025](https://elifesciences.org/articles/97981#s4-9-1)
* data made available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241987
* analysis only done on [Sample GSM7747185](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM7747185)
* note, I use [Rye](https://rye.astral.sh/) for package management.

### Step 1. Access data
* Download .tsv & .mtx data. This is a CellRanger output, as described [here](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices).
* Unzip data as necessary using `gzip -d filename.gz`, and save to a desired directory (here, `./data/GSM7747185_Chow-eWAT`). Rename corresponding files to `barcodes.tsv`, `genes.tsv`, and `matrix.mtx`.

### Step 2. CellBender (Python [CellBender, PyTables])
* Used to clean technical artifacts from sequencing data.
* This can take a while, so run it on an HPC if necessary. You'll need [Cellbender](https://cellbender.readthedocs.io/en/latest/installation/index.html), and should follow this [tutorial](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for more details.
```
cellbender remove-background --input ./data/GSM7747185_Chow-eWAT --output cleaned_output.h5
```
* There's an issue document on the [README of the GitHub page](https://github.com/broadinstitute/CellBender), regarding checkpoint/saving issues on v0.3.1. This is echoed [here](https://github.com/broadinstitute/CellBender/issues/386), and the suggested solution was to pull from this [recent commit](https://github.com/lukabor/CellBender/commit/4334e8966217c3591bf7c545f31ab979cdc6590d):
```
pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@4334e8966217c3591bf7c545f31ab979cdc6590d
```
(Hopefully, this is resolved in v0.3.2.)

* To complete analysis using Seurat in R (see [the tutorial](https://cellbender.readthedocs.io/en/latest/tutorial/index.html#open-in-seurat)) you must also run this command, which requires [PyTables](https://www.pytables.org/usersguide/installation.html):
```
ptrepack --complevel 5 data/CellBender/cleaned_output_filtered.h5:/matrix data/CellBender/cleaned_output_filtered_seurat.h5:/matrix
```

### Step 3. Data Processing (R [Seurat, DoubleFinder]).
* See the .Rmd file for details. Broadly, the steps are:
1. Process each dataset:
    1. Pull the CellBender data into Seurat.
    2. Filter samples based on UMIs (counts), genes (features), mitochrondial gene ratio, and UMIs/gene.
    3. Run [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) without ground truth, from [McGinnis et al. 2019](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30073-0).
2. Pool datasets togther and integrate.

### Step 4a. Cluster Identification
* Identify graph-based clusters, visualize (PCA, UMAP, tSNE), and annotate clusters by marker genes using Seurat.
* For improved granularity, select a subset of cells within a clusters, then cluster again and annotate using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html), from [Wu et al. 2021](https://www.sciencedirect.com/science/article/pii/S2666675821000667?via%3Dihub).

### Step 4b. Differentially Expressed Genes (R & webtools [Morpheus & Metascape])
* Group counts by sample and condition using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* Identify DEGs based on 0.05FDR and absolute 2FC.
* Cluster genes within each condition based on K-means using [Morpheus](https://software.broadinstitute.org/morpheus/), from the Broad Institute.
* Pathway analysis for each gene cluster using [Metascape](https://metascape.org/gp/index.html#/main/step1), from [Zhou et al. 2019](https://www.nature.com/articles/s41467-019-09234-6).
