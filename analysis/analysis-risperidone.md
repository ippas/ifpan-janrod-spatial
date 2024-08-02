# Analysis for risperidone (12 samples)

First analysis using spaceranger using command:

```
bash preprocessing/spaceranger-analysis.sh \
  --transcriptome raw/refdata-gex-mm10-2020-A/ \
  --output-dir data/spaceranger/ \
  --metadata data/metadata-antipsychotics.tsv \
  --samples `cat data/metadata-antipsychotics.tsv | grep -P "risperidone|saline" | sed 1d | cut -f1 | tr "\n" "," | sed 's/,$/\n/'`
  --localcores 40
```

### Analysis for risperidone, opus20, (12 samples) 
Merged samples for risperidone:
```
samtools merge \
  -@ 42  \
  data/risperidone/gene-annotation/merged-samples.bam `cat data/metadata-antipsychotics.tsv  |  \
    grep -P "saline|risperidone" | \
    cut -f1 | \
    xargs -i bash -c 'find data/spaceranger/{}/SPATIAL_RNA_COUNTER_CS/ -name *.bam'` 
```

Prepared index file for merged-samples.bam using:
```
samtools index data/risperidone/gene-annotation/merged-samples.bam
```

##### Analysis [MACS3](https://github.com/macs3-project/MACS)
Command to run macs3 analysis using [matzieb/macs3](https://hub.docker.com/repository/docker/matzieb/macs3/general):
```
docker run \
--rm \
 -u 10007:10000 \
  --name macs3 \
  -v $PWD:/ifpan-janrod-spatial/ \
  matzieb/macs3 \
  macs3 callpeak \
  -t /ifpan-janrod-spatial/data/risperidone/gene-annotation/merged-samples.bam \
  --outdir /ifpan-janrod-spatial/data/risperidone/gene-annotation/macs3
```

The structure of a peak file consists of several columns:
1. Chromosoem (chr): The chromosome on which the peak is located
2. Start position (start)
3. End position (end)
4. Peak name (name): A unique identifier for each peak, ususally in the format "prefix_peak_#", where # is a consecutive number.
5. Fold enrichment (score): The ration fo the ChIP-seq signal in the peak region compared to the control (input) sample. This value indicates the strength of the ChIP signal at the peak.
6. Strand
7. Peak summint (summit): The genomic coordinate of the highest point within the peak region, represtnting the most likely binding site of the protein of interest.
8. -log10(p-value)(pVAlue): The negative log10-transformed p-value, indicating the statistical significantce of the peak. 
9. Fold enrichment vs. control (foldChange): The fold enrichment of the ChIP-Seq signal at the peak summit compared to the local backgroud signal.
10. -log10(q-value)(qValue): The negative log10-transformed q-value (or FDR, false discovery rate), which is the p-value adjusted for multiple testing. 


#### Annotate trascripts

Run [prepare-annotate-peaks.sh](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/prepare-annotate-peaks.sh) using command:
```
bash preprocessing/prepare-annotate-peaks.sh   \
  --gene-bed data/mart-export-v102-mm10.bed  \
  --ltr-bed data/ltr-grcm38-mm10.bed  \
  --macs-peaks data/risperidone/gene-annotation/macs3/risperidone_peaks.narrowPeak \
  --data-dir data/risperidone/gene-annotation/ \
  --number-threads 40
```


#### Filter peaks using results from nanopore
To identify the cDNA associated with risperidone in nanopore sequencing data, the following command was used:
```
bedtools intersect \
  -a data/risperidone/gene-annotation/peaks-annotate-sort.bed \
  -b data/str-cdna-nanopore.bed \
  -wa -wb | \
    awk '{print $0, ($16==$22?_:"NOT") "MATCH"}' | \
    grep -v NOT | \
    cut -f 1-17 | \
    sort | \
    uniq \
    > data/ldopa/gene-annotation/peaks-annotate-filt-nanopore.bed
```


Co zrobić gdy jeden pik przypisany do dwóch genów?
When one peak was describe to two genes, the genes was merged or remov

#### Reduction peak
After annotating peaks, run the [reduction-peaks.R](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/reduction-peaks.R) script, which is responsible for filtering out peaks that are classified as artificial in DNA sequencing:
```
docker run \
  -v `pwd`:`pwd` \
  -u 10007:10000 \
  --rm  matzieb/rstudio-4.1.2-bioconductor-peaks  \
    Rscript `pwd`/preprocessing/reduction-peaks.R \
      --file_info_peaks `pwd`/data/risperidone/gene-annotation/peaks-annotate-filt-nanopore.bed \
      --bam_file `pwd`/data/risperidone/gene-annotation/merged-samples.bam \
      --output_tsv `pwd`/data/risperidone/gene-annotation/peaks-annotate-reduction2.tsv \
      --number_cores 40 \
      --peak_counts 1200 \
      --amplitude_peak 400 \
      --ratio_counts_amplitude 1.4
```

It is possible execute this step by run rstudio:
```
docker run \
  -d   \
  -p 9777:8787   \
  -e USER=[your user]  \
  -e USERID=$(id | grep -o -P "uid=\d+" | cut -d '=' -f2)   \
  -e GROUPID=[your groupid]   \
  -e UMASK=002   \
  -e PASSWORD=[your password]   \
  -v `pwd`:`pwd` \
  --name [name of docker container]  matzieb/rstudio-4.1.2-bioconductor-peaks
```



#### Create gtf file for peaks

Using  `bed`  file with peaks created  `gtf`  file using [bed2gtf-spaceranger.py](https://github.com/ippas/ifpan-janrod-spatial/blob/master/preprocessing/bed2gtf-spaceranger.py):

```
python3 preprocessing/bed2gtf-spaceranger.py \
  --input data/risperidone/gene-annotation/peaks-annotate-reduction.tsv \
  --output data/risperidone/gene-annotation/peaks-annotate-reduction.gtf
```

#### Spaceranger mkref - corrected

Created a corrected transcriptome reference employing  [spaceranger mkref](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/advanced/references) using command:
```
spaceranger mkref \
  --genome=corrected-reference \
  --fasta=raw/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes=data/risperidone/gene-annotation/peaks-annotate-reduction.gtf && \
  mv corrected-reference data/risperidone
```

#### Spaceranger counts - corrected
The spaceranger count analysis was re-executed using the updated transcriptome reference for risperidoned samples, using command:
```
bash preprocessing/spaceranger-analysis.sh \
  --transcriptome data/risperidone/corrected-reference \
  --output-dir data/spaceranger/ \
  --metadata data/metadata-antipsychotics.tsv \
  --samples `cat data/metadata-antipsychotics.tsv | sed 1d | cut -f1 | tr "\n" "," | sed 's/,$/\n/'`
  --localcores 40
```
