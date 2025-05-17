# snRNAseq-analysis

Replicating analysis by [So et al. 2025](https://elifesciences.org/articles/97981#s4-9-1)
* snRNAseq data was made available here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241987
* note, I use [Rye](https://rye.astral.sh/) for package management.
* [**Click me to see how what step I'm working on now.**](#anchor)

### Repo Structure
``` py
.
└── snRNAseq-analysis
    ├── [.venv] #gitignored
    ├── data #
    │   ├── [CellBender]
    │   │   ├── GSM7747185
    │   │   │   ├── output.h5
    │   │   │   ├── output_filtered.h5
    │   │   │   └── #...
    │   │   └── #...
    │   ├── [CellRanger]
    │   │   ├── GSM7747185
    │   │   └── #...
    │   ├── raw # large FASTQ files, and the .sh to get them
    │   │   ├── PRJNA1010853-fastq.sh
    │   │   ├── SRR25819152_1.fastq.gz
    │   │   └── #...
    │   └── [preprocessed] #provided by authors
    │   │   ├── GSM7747185_Chow-eWAT
    │   │   │   ├── barcodes.tsv
    │   │   │   ├── genes.tsv
    │   │   │   └── matrix.mtx
    │   │   └── #...
    ├── src
    │   ├── analysis.Rmd
    │   └── cellbending.sh
    ├── snRNAseq_analysis.Rproj
    ├── .gitignore
    ├── pyproject.toml          # Rye project management
    ├── requirements-dev.lock   # Rye project management
    └── requirements.lock       # Rye project management
```

### Step 0. Access data & dependencies
* Download .tsv & .mtx data. This is a CellRanger output, as described [here](https://kb.10xgenomics.com/hc/en-us/articles/115000794686-How-is-the-MEX-format-used-for-the-gene-barcode-matrices).
* Unzip data as necessary using `gzip -d [filename].[ext].gz`. Rename the corresponding files to `barcodes.tsv`, `genes.tsv`, and `matrix.mtx`, and save to a desired directory (using the repo structure above, I used `./data/CellRanger/[filename]/`). 
* [CellRanger](https://www.10xgenomics.com/support/jp/software/cell-ranger/latest/tutorials/cr-tutorial-in)
    * Also, get the relevant [transcriptomes](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads#reference-downloads) (authors used mm10/GRCm38, but I'm using mm39/GRCm39 as of May 2025). You can also build a custom one here following these [instructions](https://www.10xgenomics.com/support/jp/software/cell-ranger/downloads/cr-ref-build-steps).

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

* R Dependencies:\
Can be found in the Rmd file too
    * BioConductor:
        * Seurat
        * patchwork
        * hdf5r
    * General:
        * remote
        * here
        * dplyr
    * Remote:
        * DoubletFinder (use `install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)`)

### Step 1. CellRanger
* `cellranger count` used to align/map FASTQ reads. Authors use CellBender to call cells, but CellRanger has this functionality. See here for a [discussion of differences](https://bioinformatics.stackexchange.com/questions/20497/how-do-cellranger-and-cellbender-call-cells-what-is-the-difference-between-them).

### Step 2. CellBender (Python [CellBender, PyTables])
* Used to clean technical artifacts from sequencing data. This can take a while, so run it on an HPC/server as necessary.
* See `cellbending.sh` for proper usage. Besides input/output, the other settings are just to conserve computing resources as Presidio doesn't have a scheduler.
* To complete analysis using Seurat in R (see [the tutorial](https://cellbender.readthedocs.io/en/latest/tutorial/index.html#open-in-seurat)) you must also use the command `ptrepack` for compatabilitiy, which requires PyTables:
```
ptrepack --complevel 5 data/CellBender/[filename].h5:/matrix data/CellRanger/[filename]_seurat.h5:/matrix
```


### Step 3. Data Processing (R [Seurat, DoubleFinder]).
See the .Rmd file for details. Broadly, the steps are:
1. Process each dataset:
    1. Pull the CellBender data into Seurat.
    2. Filter samples based on UMIs (counts), genes (features), mitochrondial gene ratio, and UMIs/gene.
    3. Run [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) without ground truth.
2. Pool datasets togther and integrate.
> [!NOTE]
> <a name="anchor"></a> $${\color{lightblue}I'm \space Here!}$$\
> Awaiting Presidio access to run CellBender on additional datasets.


### Step 4a. Cluster Identification
* Identify graph-based clusters, visualize (PCA, UMAP, tSNE), and annotate clusters by marker genes using Seurat.
* For improved granularity, select a subset of cells within a clusters, then cluster again and annotate using [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html).


### Step 4b. Differentially Expressed Genes (R & webtools [Morpheus & Metascape])
* Group counts by sample and condition using [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
* Identify DEGs based on 0.05FDR and absolute 2FC.
* Cluster genes within each condition based with K-means using [Morpheus](https://software.broadinstitute.org/morpheus/).
* Pathway analysis for each gene cluster using [Metascape](https://metascape.org/gp/index.html#/main/step1).

I will create a .bib at some point, but until then here are some useful reads:
[BioConductor Intro to SComics](https://bioconductor.org/books/3.13/OSCA/)
[CellRanger (Lun et al. 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)
[CellBender (Fleming et al. 2023)](https://www.nature.com/articles/s41592-023-01943-7)
[DoubletFinder (McGinnes et al. 2019)](https://www.sciencedirect.com/science/article/pii/S2405471219300730?via%3Dihub)
[SCTransform in Seurat (Jafemeister & Satija 2019)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)
[clusterProfiler (Wu et al. 2021)](https://www.sciencedirect.com/science/article/pii/S2666675821000667?via%3Dihub)
[MetaScape (Zhou et al. 2019)](https://www.nature.com/articles/s41467-019-09234-6))