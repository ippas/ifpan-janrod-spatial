# Analysis for L-DOPA (12 samples)

Merged samples for L-DOPA:
```
samtools merge \
  -@ 16  \
  data/ldopa/gene-annotation/merged-samples.bam `cat data/samples-spatial-metadata.tsv | \
    grep "tif-ldopa" | \
    cut -f1 | \
    xargs -i bash -c 'find data/spaceranger/{}/SPATIAL_RNA_COUNTER_CS/ -name *.bam'` 
```

Prepare index file for merged-samples.bam using:
```
samtools index data/gene-annotation/ldopa/merged-samples.bam
```

#### Analysis [MACS3](https://github.com/macs3-project/MACS)
Command to run macs3 analysis from docker:
```
docker run \
 -u 10001:10000 \
  --name macs3 \
  -v $PWD:/ipan-jrparkitna-spatial/ \
  matzieb/macs3 \
  macs3 callpeak \
  --name merged-samples \
  -t /ipan-jrparkitna-spatial/data/ldopa/gene-annotation/merged-samples.bam \
  --outdir /ipan-jrparkitna-spatial/data/ldopa/gene-annotation/macs3
```

#### Annotate trascripts
Run prepare-annotate-peaks.sh using command:
```
bash preprocessing/prepare-annotate-peaks.sh \
  --gene-bed data/mart-export-v102-mm10.bed \
  --ltr-bed data/ltr-grcm38-mm10.bed \
  --macs-dir data/ldopa/gene-annotation/macs3/ \
  --data-dir data/ldopa/gene-annotation/ \
  --number-threads 16
```

#### Filter peaks using results from nanopore

Reduction peaks (for 12 samples - L-DOPA)
After annotating peaks run `reduction-peaks.R` which is responsible for reduction peaks

Parameters to filtering peaks: (to change)
- score > 350
- counts > 800
- summit.log.p > 35
- amplitude > 400
- count/amplitude > 1.4
- strand peaks agree with strand gene
- distance from gene < 30000

#### Create gtf file for peaks
Using `bed` file with peaks created `gtf` file using py:
```
python3 preprocessing/bed2gtf-spaceranger.py \
  --input data/ldopa/gene-annotation/peaks-annotate-reduction.tsv \
  --output data/ldopa/gene-annotation/peaks-annotate-reduction.gtf
```

#### Spaceranger mkref - corrected
Created a corrected transcriptome reference employing [spaceranger mkref](https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/advanced/references) and a prepared `gtf` file using the command:
```
spaceranger mkref \
  --genome=corrected-reference \
  --fasta=raw/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes=data/ldopa/gene-annotation/peaks-annotate-reduction.gtf && \
  mv corrected-reference data/ldopa
```


#### Spaceranger counts - corrected
The `spaceranger count` analysis was performed again with the corrected transcriptome reference and samples for L-DOPA using the command:
```
bash preprocessing/spaceranger-analysis.sh \
  --transcriptome some \
  --output-dir some2 \
  --samples `cat data/samples-spatial-metadata.tsv | grep "tif\-ldopa" | cut -f1 | tr "\n" "," | sed 's/,$/\n/'`
