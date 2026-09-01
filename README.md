# ifpan-janrod-spatial

This repository contains scripts and analysis workflows for spatial transcriptomics data from three mouse experiments investigating pharmacological and genotype-dependent effects on gene expression in the brain.

- **Reference data used in the analyses** – genomic and transcriptomic reference datasets used during preprocessing.
- **Analyses** – experimental design, analysis workflow and publication information for each dataset.
- **Interactive spatial transcriptomics browser** – interactive exploration of spatial clusters, gene expression and statistical results.

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

<!-- Backup of old documentation -->
##### Backup of old documentation:
<br>





# ifpan-janrod-spatial

#### Project logline (technique, organism, tissue type)
Short description of treatment groups/subjects


## Methods
This sections should be a description of preprocessing and analysis ready to be included in the publication


## Preprocessing
#### Information about scripts
##### [spaceranger-analysis.sh](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/spaceranger-analysis.sh)
Script automates the execution of the `spacerange count` analysis for multi-sample using `samples-spatial-metadata.tsv` file contains metadata about samples.

Script arguments:
- -\-transcriptome - indicate folder with reference of transcriptome
- -\-output-dir - indicate output dir of results
- -\-samples - indicate for samples execute spaceranger count analysis. Set many samples using comma separated
- -\-metadata - flag is used to provide a metadata TSV file for the analysis
- -\-localcores - flag is used to specify the number of CPU cores to be allocated for parallel processing during the analysis.

##### [prepare-annotate-peaks.sh](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/prepare-annotate-peaks.sh) 
The script contains commands which prepare a file with annotated peaks:
 - To assess peaks detected by `macs3` was used bedtools intersect and `ltr-grcm38-mm10.bed` file. If the peak intersected with any ltr in `bed` file was assigned to ltr if not was assigned to gene ant results were saved to `bed` file. 
 - Using samtools prepare a `bam` file and index `bam.bai` them with minus ad plus strand for `merged-samples.bam`.
 - For prepared earlier files with ltr and gene peaks from macs3, extracted information about coverage of strand using `samtools bedcov` using plus mand minus `bam` files. 
 - Based on which strand was greater coverage assessment strand of peaks.
 - Sorted `mart-export-v102-mm10.bed` using `bedtools sorted`.
 - Using `bedtools insert` and `mart-export-v102-mm10.bed` assigned gene name for ltr and gene peaks in the range +/- 30000 from gene.
 - Merged two last created files and sorted them.
 
These steps are aimed at preparing a file containing information about peaks that a peak belongs to gene or ltr, information to the nearest gene, and also information about a strand of a peak, therefore, coverage was counted for each peak.

Script arguments:
  - -\-gene-bed - Path to the gene BED file
  - -\-ltr-bed - Path to the LTR BED file
  - -\-macs-dir - Path to the peaks file obtained by MACS3
  - -\-data-dir	- Path to the data directory
  - -\-number-threads - Number of threads to use (default: 1)

##### [bed2gtf-spaceranger.py](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/bed2gtf-spaceranger.py)
The script convert `bed` to `gtf` file.

Script arguments:
  - -\-input - input file
  - -\-output - output file 


#### Information about data
1. `.fastq` files are available on the [SRA database](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1143882).
2. `.json` files for risperidone and saline are availabel [here](https://github.com/ippas/ifpan-janrod-spatial/tree/master/analysis/json-files)
3. Brightfield H&E tissue images (`.tif` files) are available upon request.


## Analysis
### Prepare spaceranger

Spaceranger 1.3.1. [download from](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation)

Command to initial spaceranger:
```
export PATH=/home/mateusz/projects/ifpan-janrod-spatial/tools/spaceranger-1.3.1:$PATH
```

Referance genome for spaceranger [download from](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/installation)

#### Download information about ltr and gene. 
Prepare `ltr-grcm38-mm10.bed file` [download from](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1285659565_Hc0UpfZ6D2O5aH8QfBBNIGYqANc9&clade=mammal&org=Mouse&db=mm10&hgta_group=varRep&hgta_track=joinedRmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56%2C694%2C976-56%2C714%2C605&hgta_outputType=primaryTable&hgta_outFileName=) using parameters:
 -  clade: mammal
 - genome: Mouse
 - assembly: Dec. 2011 (GRCm38/mm10)
 - group: Variation and Repeats
 - track: RepeatMasker

Prepare `mart-export-v102-mm10.bed` file [download from](http://nov2020.archive.ensembl.org/biomart/martview/41fc32d9a3d3d980eaf9f536c5256275) with:
 - chromosome
 - gene start
 - gene end
 - gene stable ID
 - gene name
 - strand
 

### Prepare data from nanopore
The [repository](https://github.com/ippas/ifpan-janrod-nanopore) provides information on how was prepared transcripts for nanopore sequencing.

A `bed` file was  prepared from a `csv` file containing transcirpts find in striatum by nanopore sequencing experiment using command:
```
cat data/str-cdna-peaks.csv | 
  grep -v MT | \
  sed '1d ;s/,/\t/g' | \
  awk 'BEGIN{OFS="\t"}{print "chr"$2, $6-15000, $6+15000, $1, $3}' \
  > data/str-cdna-nanopore.bed
```


### Analysis spaceranger count

To execute the first analysis in spaceranger run command:
```
bash preprocessing/spaceranger-analysis.sh \
  --transcriptome raw/refdata-gex-mm10-2020-A/ \
  --output-dir data/spaceranger/
  --samples `cat data/samples-spatial-metadata.tsv | sed 1d | cut -f1 | tr "\n" "," | sed 's/,$/\n/'
```

### [Analysis ldopa](https://github.com/ippas/ifpan-janrod-spatial/blob/master/analysis/analysis-ldopa.md)

### [Analysis risperidone](https://github.com/ippas/ifpan-janrod-spatial/blob/master/analysis/analysis-risperidone.md)


### Example command to run shiny app
```
docker run --rm -p 4000:3838  \
  -v /home/mateusz/projects/ifpan-janrod-spatial/preprocessing/functions:/home/rstudio/preprocessing/functions \
  -v /home/mateusz/projects/ifpan-janrod-spatial/data:/home/rstudio/data \
  -v /home/mateusz/projects/ifpan-janrod-spatial/results:/home/rstudio/results \
  -v /home/mateusz/projects/ifpan-janrod-spatial/preprocessing/shiny-app/spatial-browser-v2.2-deploy:/srv/shiny-server/ \
  -v /home/mateusz/projects/ifpan-janrod-spatial/preprocessing/shiny-app/spatial-browser-v2.2-deploy:/etc/shiny-server/   matzieb/shiny-4.0.3-seurat4-spatial
```

## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present
- tools - downoload tools
