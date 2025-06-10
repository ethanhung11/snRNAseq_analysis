# snRNAseq-analysis


Replicating analysis by [So et al. 2025](https://elifesciences.org/articles/97981#s4-9-1)
* snRNAseq data was made available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241987
* note, I use [Rye](https://rye.astral.sh/) for package management.
* [**Click me to see how what step I'm working on now.**](#anchor)

### Repo Structure
``` py
.
└── snRNAseq-analysis
    ├── .venv #gitignored
    ├── data #gitignored
    │   ├── cellbender
    │   │   ├── experiment_set1
    │   │   │   ├── experiment1-ANNOTATION-sample1
    │   │   │   │   ├── output_filtered.h5
    │   │   │   │   ├── output_filtered_seurat.h5
    │   │   │   │   └── #cellbender outputs...
    │   │   │   └── #...
    │   │   └── #...
    │   ├── cellranger
    │   │   ├── experiment_set1
    │   │   │   ├── experiment1
    │   │   │   │   ├── sample1
    │   │   │   │   │   └── #cellranger outputs...
    │   │   │   │   └── #...
    │   │   │   └── #...
    │   │   └── #...
    │   ├── raw
    │   │   ├── experiment1
    │   │   │   ├── sample1.fastq.gz
    │   │   │   └── #...
    │   │   └── #...
    │   └── preprocessed
    │   │   ├── object.pickle
    │   │   ├── object.Rds
    │   │   ├── object.h5ad
    │   │   ├── object.h5
    │   │   └── #...
    ├── plots                   # plots if desired
    ├── references              # ref genomes, databases, etc.
    ├── src                     # Anything starting with a # is not yet implemented
    │   ├── Python
    │   │   ├── preprocess.py
    │   │   ├── analysis.py
    │   │   └── R.py
    │   ├── R
    │   │   ├── preprocesss.R
    │   │   ├── preprocesssing.Rmd
    │   │   ├── #analysis.R
    │   │   ├── analysis.Rmd
    │   │   └── pipeline.Rmd
    │   └── scripts
    │       ├── download_data.sh
    │       ├── cellcounting.sh
    │       ├── cellbending.sh
    │       ├── #preprocess-py.sh
    │       ├── #analyze-py.sh
    │       ├── preprocess-R.sh
    │       └── #analyze-R.sh
    ├── snRNAseq_analysis.Rproj
    ├── .gitignore
    ├── pyproject.toml          # Rye project management
    ├── requirements-dev.lock   # Rye project management
    └── requirements.lock       # Rye project management
```

# Steps
In this example, I have raw data from 4 experiments `GSM7747185` through `GSM7747188` (each is it's own condition, and has 1-2 samples each). These 4 experiments are grouped into 1 experiment set called `paper_processed`.

### Step 0a. Access Data
* **PREPROCESSED DATA:**
    * Download .tsv & .mtx data. This is a CellRanger output, as described [here](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices).
    * Unzip data as necessary using `gzip -d [filename].[ext].gz`. Rename the corresponding files to `barcodes.tsv`, `genes.tsv`, and `matrix.mtx`, and save to a desired directory (using the repo structure above, I used `./data/CellRanger/[filename]/`). 
    * Skip to [Step 3](#step-3-data-processing-r-seurat-doublefinder)!

* **UNPROCESSED DATA:**
    * Pull data from European Nucleotide Archive. Search by accession for the to find this [site](https://www.ebi.ac.uk/ena/browser/view/PRJNA1010853). The downloaded script is in `./data/preprocessed/PRJNA1010853-fastq.sh`, which you can run to pull the ~170GB FASTQ dataset.
    * You'll also need [CellRanger](https://www.10xgenomics.com/support/jp/software/cell-ranger/latest/tutorials/cr-tutorial-in) to align the sequencing data.
        * You'll also need the relevant [transcriptomes](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads#reference-downloads) (authors used mm10/GRCm38, but I'm using mm39/GRCm39 as of May 2025). You can also build a custom one here following these [instructions](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads/cr-ref-build-steps).
    * Python Dependencies:
        * [Cellbender](https://cellbender.readthedocs.io/en/latest/installation/index.html); should follow this [tutorial](https://cellbender.readthedocs.io/en/latest/tutorial/index.html) for more details.
            * There's an issue document on the [README of the GitHub page](https://github.com/broadinstitute/CellBender), regarding checkpoint/saving issues on v0.3.1. A discussion on this can be found [here](https://github.com/broadinstitute/CellBender/issues/386). The suggested solution was to pull from this recent commit: [`4334e8966217c3591bf7c545f31ab979cdc6590d`](https://github.com/lukabor/CellBender/commit/4334e8966217c3591bf7c545f31ab979cdc6590d):\
            For rye:
            ```
            rye add cellbender --git=https://github.com/broadinstitute/CellBender.git@4334e8966217c3591bf7c545f31ab979cdc6590d
            ```
            or normal pip installation:
            ```
            pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@4334e8966217c3591bf7c545f31ab979cdc6590d
            ```
        * [PyTables](https://www.pytables.org/usersguide/installation.html)

### Step 0b: Prepare Analysis Dependencies
* R Dependencies:
    * anything in the scripts
        * Seurat (and other Seurat packages, see [here](https://satijalab.org/seurat/articles/install.html#seurat-v5-seurat-5-install-from-github))
        * DoubletFinder (use `chris-mcginnis-ucsf/DoubletFinder`)
        * presto (use `immunogenics/presto`)
Example on how to install remote Github package:
```
library(remotes)
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
```

### Step 1. Generate Count Matrix (CellRanger)
* `cellranger count` used to align/map FASTQ reads.
* See `cellcounting.sh`, which searches `./data/raw` for various experiments (each of which should have the fastq.gz files inside).
* Authors use CellBender to call cells, but CellRanger has this functionality. See here for a [discussion of differences](https://bioinformatics.stackexchange.com/questions/20497/how-do-cellranger-and-cellbender-call-cells-what-is-the-difference-between-them).
To run `cellcounting.sh`, see the example below:
```
./src/scripts/cellcounting.sh -i GSM7747185,GSM7747186,GSM7747187,GSM7747188 -o cellcounting.out
```

### Step 2. Call Cells (CellBender)
* `cellbender remove-background` used to clean technical artifacts from sequencing data.
* The script checks for the cellbender env. If you have one, just change the name in the script. Otherwise, the script will create one using `mamba` (change `mamba` -> `conda` if you don't have mamba)
* See `cellbending.sh`, which searches the provided directory (e.g. experiment_set1) for various experiments, which should have at least one sample that's in cellranger output format.
* To complete analysis using Seurat in R (see [the tutorial](https://cellbender.readthedocs.io/en/latest/tutorial/index.html#open-in-seurat) the command `ptrepack` must also be run for compatibility, which requires PyTables. This is implemented in the shell script.
To run `cellbending.sh`, see the example below:
```
./src/scripts/cellbending.sh -i ./data/cellranger/paper_processed -o cellbending.out
```

### Step 3. Data Pre-Processing (Seurat).
See the .Rmd file for details. Broadly, the steps are:
1. Process each dataset:
    1. Pull the CellBender data into Seurat.
    2. Filter samples based on UMIs (counts), genes (features), mitochrondial gene ratio, and UMIs/gene.
    3. Run [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) without ground truth.
    4. If using Seurat, can run cell stage scoring & regress out those genes (and mt/rb genes).
2. Pool datasets togther and integrate.
3. Identify graph-based clusters & visualize (PCA, UMAP, tSNE).


### Step 4a. Cluster Identification (Seurat)
* Annotate clusters by marker genes using Seurat (see [FindMarkers](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#finding-differentially-expressed-features-cluster-biomarkers)) or in this case, markers are provided in the paper.
* For improved granularity, select a `subset` of cells within a clusters, then cluster again.
> [!IMPORTANT]
> Haven't looked into this yet!
> * Annotate using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) if desired (per the example paper).



### Step 4b. Differentially Expressed Genes & Pathway Analysis (Seurat, Morpheus, & Metascape)
* Group counts by sample and condition using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

> [!NOTE]
> <a name="anchor"></a> <span style="color: blue;">I'm here!</span>

* Identify DEGs based on 0.05FDR and absolute 2FC.
* Cluster genes within each condition based with K-means using [Morpheus](https://software.broadinstitute.org/morpheus/), per the example paper.
* Pathway analysis for each gene cluster using [Metascape](https://metascape.org/gp/index.html#/main/step1), per the example paper.

### Step 4c. Cell-Cell Interactions (NicheNet & CellChat)
> [!IMPORTANT]
> Haven't looked into this yet!

# Bibliography & Readings
I will create a .bib at some point, but until then here are some useful reads:
* [BioConductor Intro to SComics (R, but not Seurat)](https://bioconductor.org/books/3.13/OSCA/)
* [SeuratV3 tutorial (R)](https://satijalab.org/seurat/articles/pbmc3k_tutorial) (don't use the DGEs here for anything other than annotation)
* [SCVerse Best Practices (Python)](https://www.sc-best-practices.org/preamble.html)
* [CellRanger (Lun et al. 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)
* [CellBender (Fleming et al. 2023)](https://www.nature.com/articles/s41592-023-01943-7)
* [DoubletFinder (McGinnes et al. 2019)](https://www.sciencedirect.com/science/article/pii/S2405471219300730?via%3Dihub)
* [SCTransform in Seurat (Hafemeister & Satija 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
* [clusterProfiler (Wu et al. 2021)](https://www.sciencedirect.com/science/article/pii/S2666675821000667?via%3Dihub)
* [MetaScape (Zhou et al. 2019)](https://www.nature.com/articles/s41467-019-09234-6))

# A list of other problems I ran into:
- Issue: can't add CellBender into Rye
    - Solution: Install cellbender into a conda environment. `conda_cellbender.yml` attached for convenience.
- Issue: Rye installation of `package` led to clang not found error:
    - Notes: probably occurs because Presidio doesn't have clang
    - Solved by: https://github.com/astral-sh/rye/issues/836
    - Solution: `CC=$(which gcc) CXX=$(which g++) rye add [package]`
