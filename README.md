# ifpan-janrod-spatial

#### Project logline (technique, organism, tissue type)
Short description of treatment groups/subjects


## Methods
This sections should be a description of preprocessin and analysis ready to be included in the publication


## Preprocessing
#### Information about scripts
##### [spaceranger-analysis.sh](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/spaceranger-analysis.sh)
Script automates the execution of the `spacerange count` analysis for multi-sample using `samples-spatial-metadata.tsv` file contains metadata about samples.

Script arguments:
- -\-transcriptome - indicate folder with reference of transcriptome
- -\-output-dir - indicate output dir of results
- -\-samples - indicate for samples execute spaceranger count analysis. Set many samples using comma separated

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
  - -\-gene-bed - file contain infromation about gene
  - -\-ltr-bed - file contain information about ltr
  - -\-macs-dir - indicate path to results from MACS3
  - -\-data-dir	- indicate where are needed files and output dir for `peaks-annotate-sort.bed`
  - -\-number-threads - number of threads to use 

##### [bed2gtf-spaceranger.py](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/bed2gtf-spaceranger.py)
The script convert `bed` to `gtf` file.

Script arguments:
  -  -\-input - input file
  - -\-output - output file 

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




## About this template
Directories:
- _root_ - README.md, *.Rproj, general configuration files, etc.
- raw - raw data
- preprocessing - scripts
- data - useful data, created by scripts/tools/preprocessing
- analysis - analysis source code
- results - output ready to present
- tools - downoload tools
