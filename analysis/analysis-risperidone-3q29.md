
# 1. Data preprocessing

## 1.1. Initial read alignment and gene counting

The initial read alignment and gene counting were performed using **Space Ranger v3.1.3** and the standard 10x Genomics mouse transcriptome reference (`refdata-gex-mm10-2020-A`).

**Script:** [`spaceranger-analysis-v2.sh`](../preprocessing/spaceranger-analysis-v2.sh)

```bash
bash preprocessing/spaceranger-analysis-v2.sh \
  --transcriptome raw/refdata-gex-mm10-2020-A/ \
  --output-dir data/risperidone-3q29/spaceranger-raw \
  --metadata data/risperidone-3q29/metadata-3q29-ris.tsv \
  --samples $(cat data/risperidone-3q29/metadata-3q29-ris.tsv | sed 1d | tail -4 | cut -f1 | tr "\n" "," | sed 's/,$//') \
  --input-directory raw/risperidone-3q29 \
  --localcores 20
```

The script runs spaceranger count using the specified transcriptome reference (--transcriptome), output directory (--output-dir), sample metadata (--metadata), selected sample IDs (--samples), input sequencing data directory (--input-directory), and number of processing cores (--localcores). For each sample, the metadata file provides the corresponding Visium slide ID, capture area, tissue-positioning JSON file and sample-specific information required to match the sequencing data with the spatial image. FASTQ files are retrieved from the specified input directory.

---

## 1.2. Preparation and annotation of sequenced transcript regions

The purpose of this step was to identify transcribed genomic regions not captured by the default reference annotation, determine their overlap with protein-coding genes from BioMart v102 and LTR elements, and use this information to construct a corrected transcriptome reference.

BAM files prepared by Space Ranger for all 26 samples were combined into a single alignment file and indexed using **samtools v1.19.2**. Combining all samples provided an experiment-wide collection of aligned reads and allowed a common set of sequenced genomic regions to be identified across the complete dataset.

```bash
samtools merge \
  -@ 42 \
  data/risperidone-3q29/gene-annotation/merged-26samples-3q29-ris.bam \
  `cat data/risperidone-3q29/metadata-3q29-ris.tsv | \
    grep -P "saline|risperidone" | \
    cut -f1 | \
    xargs -i bash -c 'find data/risperidone-3q29/spaceranger-raw/{}/SPATIAL_RNA_COUNTER_CS/ -name *.bam'`

samtools index \
  data/risperidone-3q29/gene-annotation/merged-26samples-3q29-ris.bam
```

**MACS3 v3.0.0a6** was then applied to the merged BAM file to identify genomic regions with concentrated sequencing signal.

```bash
docker run \
  --rm \
  -u 10007:10000 \
  --name macs3 \
  -v $PWD:/ifpan-janrod-spatial/ \
  matzieb/macs3 \
  macs3 callpeak \
  -t /ifpan-janrod-spatial/data/risperidone-3q29/gene-annotation/merged-26samples-3q29-ris.bam \
  --outdir /ifpan-janrod-spatial/data/risperidone-3q29/gene-annotation/macs3
```

The identified regions were subsequently annotated to determine their genomic context, strand and associated gene.

**Script:** [`prepare-annotate-peaks-v2.sh`](../preprocessing/prepare-annotate-peaks-v2.sh)

```bash
bash preprocessing/prepare-annotate-peaks-v2.sh \
  --gene-bed data/mart-export-v102-mm10.bed \
  --ltr-bed data/ltr-grcm38-mm10.bed \
  --macs-peaks data/risperidone-3q29/gene-annotation/macs3/ris3q29-26s_peaks.narrowPeak \
  --data-dir data/risperidone-3q29/gene-annotation/ \
  --number-threads 20 \
  --prefix ris3q29 \
  --bam data/risperidone-3q29/gene-annotation/merged-26samples-3q29-ris.bam
```

The script:

- intersects the detected regions with the LTR annotation;
- prepares strand-specific BAM files;
- calculates coverage separately for the positive and negative strands;
- assigns strand information based on sequencing coverage;
- compares detected regions with the mouse gene annotation and assigns an associated gene within a ±30 kb interval;
- combines the resulting annotations into a single BED file for subsequent filtering and reference construction.

---


## 1.3. Peak filtering and transcript annotation reduction

Nanopore cDNA sequencing data from the striatum were used as an independent transcript-level reference to retain only regions with evidence of transcription in this brain region. Genomic interval operations were performed using **bedtools v2.31.1**.

```bash
awk 'BEGIN{OFS="\t"} {$4 = gensub(/^NA_/, "", "g", $4) "_" $16; print}' \
  ris3q29-peaks-annotate-sort.bed \
  > ris3q29-peaks-annotate-sort-renamed.bed
```

```bash
bedtools intersect \
  -a data/risperidone-3q29/gene-annotation/ris3q29-peaks-annotate-sort-renamed.bed \
  -b data/str-cdna-nanopore.bed \
  -wa -wb | \
  awk '{print $0, ($16==$22?_:"NOT") "MATCH"}' | \
  grep -v NOT | \
  cut -f 1-17 | \
  sort | \
  uniq \
  > data/risperidone-3q29/gene-annotation/peaks-annotate-filt-nanopore.bed
```

The remaining regions were further reduced to remove low-coverage and weakly supported transcript regions. The filtering used a minimum peak count of **1200**, a minimum peak amplitude of **400**, and a count-to-amplitude ratio threshold of **1.4**. This step reduced the annotation to regions with sufficient sequencing support for inclusion in the corrected reference.



**Script:** [`reduction-peaks.R`](../preprocessing/reduction-peaks.R)

```bash
docker run \
  -v `pwd`:`pwd` \
  -u 10007:10000 \
  --rm \
  --memory=150g \
  matzieb/rstudio-4.1.2-bioconductor-peaks \
  Rscript `pwd`/preprocessing/reduction-peaks.R \
  --file_info_peaks `pwd`/data/risperidone-3q29/gene-annotation/peaks-annotate-filt-nanopore.bed \
  --bam_file `pwd`/data/risperidone-3q29/gene-annotation/merged-26samples-3q29-ris.bam \
  --output_tsv `pwd`/data/risperidone-3q29/gene-annotation/peaks-annotate-reduction.tsv \
  --number_cores 10 \
  --peak_counts 1200 \
  --amplitude_peak 400 \
  --ratio_counts_amplitude 1.4
```

---

## 1.4. Preparation of the customized reference, read alignment and gene counting

The filtered transcript regions were used to prepare a customized transcriptome reference for the risperidone × 3q29 dataset. The reduced annotation was first converted to GTF format and then combined with the mouse genome sequence to create the customized Space Ranger reference.

**Script:** [`bed2gtf-spaceranger.py`](../preprocessing/bed2gtf-spaceranger.py)

```bash
python3 preprocessing/bed2gtf-spaceranger.py \
  --input data/risperidone-3q29/gene-annotation/peaks-annotate-reduction.tsv \
  --output data/risperidone-3q29/gene-annotation/peaks-annotate-reduction.gtf
```

The customized reference was prepared using **Space Ranger v3.1.3**.

```bash
spaceranger mkref \
  --genome=corrected-reference \
  --fasta=raw/refdata-gex-mm10-2020-A/fasta/genome.fa \
  --genes=data/risperidone-3q29/gene-annotation/peaks-annotate-reduction.gtf && \
  mv corrected-reference data/risperidone-3q29
```

All samples were then processed with **Space Ranger v3.1.3** using the customized reference. Space Ranger aligned sequencing reads to the updated transcriptome reference, assigned reads to annotated genes, and prepared spatial gene expression count matrices used for downstream Seurat analysis.

**Script:** [`spaceranger-analysis-v2.sh`](../preprocessing/spaceranger-analysis-v2.sh)

```bash
bash preprocessing/spaceranger-analysis-v2.sh \
  --transcriptome data/risperidone-3q29/corrected-reference \
  --output-dir data/risperidone-3q29/spaceranger-corrected/ \
  --metadata data/risperidone-3q29/metadata-3q29-ris.tsv \
  --samples $(cat data/risperidone-3q29/metadata-3q29-ris.tsv | sed 1d | tail -4 | cut -f1 | tr "\n" "," | sed 's/,$//') \
  --input-directory raw/risperidone-3q29 \
  --localcores 20
```

# 2. Data integration and clustering

## 2.1. Seurat integration and clustering

Spatial gene expression counts and associated metadata prepared by Space Ranger were processed and integrated using **Seurat v4.4.0** to identify transcriptionally similar groups of spots across samples.

Sample **S13839Nr3** was excluded due to substantially lower UMI counts compared with the remaining samples and persistent separation from other samples after integration and clustering.

**Main script:** [`ris3q29-analysis-seurat-v2.R`](../preprocessing/risperidone-3q29/ris3q29-analysis-seurat-v2.R)

**Functions:** [`functions-spatial-data.R`](../preprocessing/functions/functions-spatial-data.R)

The main processing steps included:

- loading gene expression counts, sample metadata and spatial information prepared by Space Ranger;
- normalization of each sample using `LogNormalize`;
- identification of **2,000 highly variable genes** using the variance-stabilizing transformation (`vst`) method;
- selection of **2,000 integration features**;
- identification of integration anchors using `FindIntegrationAnchors` based on canonical correlation analysis (CCA);
- integration of samples using `IntegrateData` to prepared a common expression space and reduce sample-specific technical differences;
- scaling and centering of the integrated expression data using `ScaleData`;
- principal component analysis (PCA) using **30 principal components**;
- UMAP dimensionality reduction using PCA dimensions **1–30**;
- construction of nearest-neighbor and shared nearest-neighbor graphs using PCA dimensions **1–30** and **k = 20** nearest neighbors;
- graph-based clustering of spatial spots across multiple clustering resolutions.

## 2.2. Cluster marker identification and selection of the clustering resolution

Marker genes were identified for each cluster to support anatomical interpretation of the clustering results.

**Marker identification:** [`ris3q29-find-markers.R`](../preprocessing/risperidone-3q29/ris3q29-find-markers.R)

Marker genes were identified by comparing spots assigned to each cluster with all remaining spots using the **Wilcoxon rank-sum test** implemented in Seurat. Markers were required to show a minimum **log2 fold change of 0.5**, to be detected in at least **25% of spots** within the respective cluster, and to have **FDR < 0.05**.

The clustering resolution was selected based on a combination of:

- [**Average Silhouette Width (ASW)**](https://doi.org/10.1186/s13059-024-03361-0) to assess cluster separation and within-cluster similarity, with higher values indicating better-defined clusters;
- [**Percentage of Abnormal Spots (PAS)**](https://doi.org/10.1186/s13059-024-03361-0) to assess local spatial homogeneity, with lower values indicating fewer spots inconsistent with their surrounding cluster;
- [**Spatial Chaos Score (CHAOS)**](https://doi.org/10.1186/s13059-024-03361-0) to assess spatial continuity, with lower values indicating more spatially coherent clusters;
- spatial localization of the identified clusters;
- expression of cluster-enriched marker genes and their correspondence with expected anatomical structures;
- expert assessment of the anatomical organization of the brain sections.

**Spatial clustering evaluation:** [`01_evaluateClusteringSpatialContinuity_PAS_CHAOS_ASW.R`](../preprocessing/risperidone-3q29/clustering-spatial-continuity/01_evaluateClusteringSpatialContinuity_PAS_CHAOS_ASW.R)

Based on the combined quantitative and anatomical evaluation, **resolution 0.4** was selected for statistical analyses.

## 2.3. Preparation of cluster-level expression summaries

Mean gene expression was calculated separately for each sample and cluster to prepare cluster-level expression summaries. Sample-cluster combinations containing fewer than **20 spots** were excluded to avoid calculating mean expression from poorly represented clusters.

**Script:** [`ris3q29_preapreMeanExpressionPersSampleMin20spots_09.06.2026.R`](../preprocessing/risperidone-3q29/ris3q29_preapreMeanExpressionPersSampleMin20spots_09.06.2026.R)


# 3. Statistical analysis

## 3.1. Cluster-level differential expression analysis

For each sample and spatial cluster, raw gene counts from individual spots were summed to prepare pseudobulk expression profiles, with individual animals treated as biological replicates. Sample-cluster combinations containing fewer than **20 spots** were excluded. The resulting pseudobulk counts were analyzed separately for each spatial cluster using [**edgeR**](https://bioconductor.org/packages/edgeR/).

**Main script:** [`ris3q29-edgerStatistics.R`](../preprocessing/risperidone-3q29/ris3q29-edgerStatistics.R)

The statistical workflow included:

- removal of lowly expressed genes using `edgeR::filterByExpr`;
- normalization of library sizes using edgeR normalization factors;
- estimation of dispersion and fitting of negative-binomial generalized linear models using the edgeR quasi-likelihood framework;
- robust model fitting (`robust = TRUE`);
- separate statistical analysis of each spatial cluster.

The statistical model included two experimental factors:

- **mouse genotype:** wild type (`wtwt`) or 3q29 deletion (`wtdel`);
- **treatment:** saline or risperidone.

The model was used to estimate:

- the **genotype main effect**;
- the **risperidone treatment main effect**;
- the **genotype × treatment interaction**.

Marginal main effects were estimated using proportional weighting across the experimental groups.

Four pairwise contrasts were calculated in the main analysis:

- `salDel vs salWt` – genotype effect in saline-treated animals;
- `risDel vs risWt` – genotype effect in risperidone-treated animals;
- `risWt vs salWt` – risperidone effect in wild-type animals;
- `risDel vs salDel` – risperidone effect in 3q29 deletion animals.

Multiple-testing correction was performed using the Benjamini–Hochberg procedure.

## 3.2. Additional cross-condition contrasts

Two additional contrasts between experimental groups were calculated and added to the main edgeR results:

- `risDel vs salWt` – risperidone-treated 3q29 deletion animals versus saline-treated wild-type animals;
- `risWt vs salDel` – risperidone-treated wild-type animals versus saline-treated 3q29 deletion animals.

**Additional contrasts:** [`edger_cross_cell_contrasts_by_cluster_29.06.2026.R`](../preprocessing/risperidone-3q29/edgeR/edger_cross_cell_contrasts_by_cluster_29.06.2026.R)

The two additional contrasts complemented the four comparisons from the main analysis, resulting in a total of **six pairwise contrasts** between the four experimental groups.

## 3.3. Statistical output

For each cluster and tested effect or contrast, the edgeR output included:

- `logFC` – estimated log2 fold change;
- `logCPM` – average expression on the edgeR log-counts-per-million scale;
- `F` – quasi-likelihood F statistic;
- `PValue` – nominal p-value;
- `FDR` – Benjamini–Hochberg adjusted p-value.

Mean expression values for the four experimental groups (`salWt`, `salDel`, `risWt`, and `risDel`) were additionally included in the result tables as descriptive measures of gene expression level. These values were used to support biological interpretation of the edgeR results and were not used for statistical significance testing.

The resulting tables contained the genotype and treatment main effects, genotype × treatment interaction, six pairwise contrasts, and supporting group-level mean expression values for each gene and spatial cluster.

# 4. Visualization and interactive data browser

Additional scripts were used for visualization and interactive exploration of the spatial transcriptomics results.

## 4.1. Genotype × treatment interaction heatmap

Genes associated with the genotype × treatment interaction were visualized using sample-level cluster expression data.

**Heatmap script:** [`plot_heatmap_zscoreSalineType_customDist.R`](../preprocessing/risperidone-3q29/interaction-heatmap-16.10.2025/plot_heatmap_zscoreSalineType_customDist.R)

Quantile-normalized (`qn`) expression values were used, and Z-scores were calculated relative to the saline group separately for each genotype. Heatmap rows were hierarchically clustered using Pearson correlation-based distances and complete linkage.

To preserve information about the spatial origin of the expression signal, distances between genes originating from the same spatial cluster were reduced by **25%**, whereas distances between genes originating from different clusters were increased by **25%**.

The robustness of the genotype × treatment interaction was additionally evaluated using permutation-based FDR estimation. Sample group labels were repeatedly randomized and the interaction was recalculated. The frequency of significant interaction signals obtained in the permuted datasets relative to the observed data was used to estimate an empirical FDR.

**Permutation-based interaction FDR:** [`ris3q29_edgerMarginaProportional_permutationFDRInteraction.R`](../preprocessing/risperidone-3q29/edgeR/ris3q29_edgerMarginaProportional_permutationFDRInteraction.R)

## 4.2. Spatial visualization of clusters and gene expression

Cluster assignments and gene expression patterns were visualized directly on individual tissue sections.

**Visualization functions:** [`visualization-functions.R`](../preprocessing/functions/visualization-functions.R)

The script includes functions for:

- visualization of spatial cluster assignments using `spatial_cluster`;
- visualization of gene expression on tissue sections using `spatial_feature_plot` and `spatial_gene_plot`.

Additional result visualizations were prepared using:

**Script:** [`ris3q29-visualization-results-test.R`](../preprocessing/risperidone-3q29/ris3q29-visualization-results-test.R)

## 4.3. Interactive spatial transcriptomics browser

An interactive Shiny application was prepared for exploration of spatial clusters, gene expression and statistical results.

**Shiny source code:** [`spatial-transcriptomics-drug-browser4`](../preprocessing/shiny-app/spatial-transcriptomics-drug-browser4/)

**Interactive browser:** [https://spatial.labpgx.com/](https://spatial.labpgx.com/)





