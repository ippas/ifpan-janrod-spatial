# ifpan-janrod-spatial

This repository contains scripts and analysis workflows for spatial transcriptomics data from three mouse experiments investigating pharmacological and genotype-dependent effects on gene expression in the brain.

- **Reference data used in the analyses** – genomic and transcriptomic reference datasets used during preprocessing.
- **Analyses** – experimental design, analysis workflow and publication information for each dataset.
- **Interactive spatial transcriptomics browser** – interactive exploration of spatial clusters, gene expression and statistical results.

#### Repository structure

* `analysis/` – experiment-specific analysis documentation;
* `preprocessing/` – preprocessing, integration, statistical analysis and visualization scripts;
* `data/` – processed data, metadata and reference files; not stored in the Git repository;
* `results/` – analysis outputs, tables and figures; not stored in the Git repository.


## 1. Reference data used in the analyses

Several reference datasets were used across the spatial transcriptomics analyses:

- **LTR annotation:** `data/ltr-grcm38-mm10.bed`, prepared from the UCSC RepeatMasker annotation for the GRCm38/mm10 mouse genome (Dec. 2011 assembly), including chromosome, start and end positions, strand, repeat name, repeat class and repeat family;
- **gene annotation:** `data/mart-export-v102-mm10.bed`, prepared using Ensembl BioMart v102 and containing chromosome, gene start, gene end, gene stable ID, gene name and strand;
- **striatal Nanopore cDNA reference:** `data/str-cdna-nanopore.bed`, prepared from transcripts detected by Nanopore cDNA sequencing of mouse striatal tissue.

The Nanopore sequencing data were prepared as part of the study [*μ-Opioid receptor transcriptional variants in the murine forebrain and spinal cord*](https://doi.org/10.1016/j.gene.2024.148890).

The processing of the Nanopore transcript data is described in the [ifpan-janrod-nanopore repository](https://github.com/ippas/ifpan-janrod-nanopore).

Detailed preparation of these reference datasets is documented below.

## 2. Analyses

### 2.1. L-DOPA

Spatial transcriptomics analysis of the mouse forebrain following L-DOPA treatment.

- **Experimental design:** TIF-IADATCreERT2 mouse model with progressive loss of dopaminergic neurons; forebrain tissue; 12 samples divided into three groups: mutant saline (n=5), mutant L-DOPA (n=5), and control saline (n=2). Mutant mice were treated with L-DOPA or saline.
- **Analysis:** [L-DOPA analysis workflow](analysis/analysis-ldopa.md)
- **Publication:** [Molecular Neurobiology, 2025](https://doi.org/10.1007/s12035-025-04957-8)

### 2.2. Risperidone

Spatial transcriptomics analysis of the mouse forebrain following acute risperidone treatment.

- **Experimental design:** adult male C57BL/6N mice; forebrain tissue; 12 samples divided into two groups: saline (n=6) and risperidone (n=6). Risperidone was administered acutely at 0.5 mg/kg, i.p.
- **Analysis:** [Risperidone analysis workflow](analysis/analysis-risperidone.md)
- **Publication:** [Frontiers in Molecular Neuroscience, 2026](https://doi.org/10.3389/fnmol.2026.1844705)

### 2.3. Risperidone × 3q29 deletion

Spatial transcriptomics analysis of the mouse brain in the 3q29 deletion model following risperidone treatment.

- **Experimental design:** wild-type and heterozygous 3q29 deletion mice; brain tissue; 26 samples divided into four experimental groups: wild-type saline (`salWt`), 3q29 deletion saline (`salDel`), wild-type risperidone (`risWt`), and 3q29 deletion risperidone (`risDel`).
- **Analysis:** [Risperidone × 3q29 analysis workflow](analysis/analysis-risperidone-3q29.md)
- **Status:** ongoing analysis.

## 3. Interactive spatial transcriptomics browser

Spatial transcriptomics datasets and analysis results can be explored using the interactive browser:

- **Browser:** [https://spatial.labpgx.com/](https://spatial.labpgx.com/)
- **Source code:** [`spatial-transcriptomics-drug-browser4`](preprocessing/shiny-app/spatial-transcriptomics-drug-browser4/)

The browser provides interactive visualization of spatial clusters and gene expression across individual tissue sections, together with access to statistical analysis results.

### 3.1. Local deployment

```bash
docker run --rm -p 4000:3838 \
  -v "$PWD/preprocessing/functions:/home/rstudio/preprocessing/functions" \
  -v "$PWD/data:/home/rstudio/data" \
  -v "$PWD/results:/home/rstudio/results" \
  -v "$PWD/preprocessing/shiny-app/spatial-transcriptomics-drug-browser4:/srv/shiny-server/" \
  -v "$PWD/preprocessing/shiny-app/spatial-transcriptomics-drug-browser4:/etc/shiny-server/" \
  matzieb/shiny-4.0.3-seurat4-spatial
```

<br>
<br>
<br>


## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present
- tools - downoload tools
