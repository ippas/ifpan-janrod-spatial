#!/bin/bash

macs_dir=data/gene-annotation/macs3
data_dir=data/gene-annotation
# create tmp folder for temporary files
tmp_dir="data/gene-annotation/tmp"
if [ ! -d "$tmp_dir" ]; then
  mkdir $tmp_dir
fi


# echo "1. Preprafe files bed with ltr ang gene peaks"
# Prepare peaks which are intersect with ltr 
# bedtools intersect \
#   -u \
#   -a $macs_dir/merged-samples_peaks.narrowPeak \
#   -b $data_dir/ltr-grcm38-mm10.bed  2>/dev/null | \
#     awk '{print $0"\tltr"}' > $tmp_dir/tmp-peaks-ltr.bed

# # Prepare peaks which are not intersect with ltr 
# bedtools intersect \
#   -v \
#   -a $macs_dir/merged-samples_peaks.narrowPeak \
#   -b $data_dir/ltr-grcm38-mm10.bed  2>/dev/null | \
#     awk '{print $0"\tgene"}' > $tmp_dir/tmp-peaks-gene.bed


# echo "2. Prepare file bam with minus and plus strand"
# # prepare bam file with minus strand
# samtools view -b \
#   -f 16 $data_dir/merged-samples.bam > $tmp_dir/tmp-merged-samples-minus.bam

# # prepare bam file with plus strand
# samtools view -b \
#   -F 16 $data_dir/merged-samples.bam > $tmp_dir/tmp-merged-samples-plus.bam


# echo "3. indexing minus and plus bam files"
# # indexing minus and plus created bam files
# samtools index $tmp_dir/tmp-merged-samples-minus.bam
# samtools index $tmp_dir/tmp-merged-samples-plus.bam

# echo "4. Prepare file with minus and plus strand coverage for ltr"
# # prepare file with minus strand coverage for ltr 
# samtools bedcov $tmp_dir/tmp-peaks-ltr.bed $tmp_dir/tmp-merged-samples-minus.bam > $tmp_dir/tmp-peaks-ltr-coverage-minus.bed
# # prepare file with plus strand coverage for ltr
# samtools bedcov $tmp_dir/tmp-peaks-ltr.bed $tmp_dir/tmp-merged-samples-plus.bam > $tmp_dir/tmp-peaks-ltr-coverage-plus.bed


# echo "5. Prepare files with minus and plus strand coverage for genes"
# # prepare file with minus strand coverage for genes
# samtools bedcov $tmp_dir/tmp-peaks-gene.bed $tmp_dir/tmp-merged-samples-minus.bam > $tmp_dir/tmp-peaks-gene-coverage-minus.bed
# # prepare file with plus strand coverage for genes
# samtools bedcov $tmp_dir/tmp-peaks-gene.bed $tmp_dir/tmp-merged-samples-plus.bam > $tmp_dir/tmp-peaks-gene-coverage-plus.bed


# echo "6. Strand assesment for ltr"
# # strand assesment for peaks ltr based on minus and plus coverage
# # paste plus and minus coverage together for ltr
# paste $tmp_dir/tmp-peaks-ltr-coverage-plus.bed $tmp_dir/tmp-peaks-ltr-coverage-minus.bed |
#     awk '{print $0"\t+"$12-$24}' |
#     awk 'BEGIN{FS=OFS"\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
#     awk -F"\t" '{OFS=FS}{ $6=$25 ; print   }' |
#     cut -f1-11 > $tmp_dir/tmp-peaks-ltr-strand.bed


# echo "7. Strand assesment for genes"
# # strand assesment fo peaks gene based on minus and plus coverage
# # paste plus and minus coverage together for genes
# paste $tmp_dir/tmp-peaks-gene-coverage-plus.bed $tmp_dir/tmp-peaks-gene-coverage-minus.bed |
#     awk '{print $0"\t+"$12-$24}' |
#     awk 'BEGIN{FS=OFS"\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
#     awk -F"\t" '{OFS=FS}{ $6=$25 ; print   }' |
#     cut -f1-11  > $tmp_dir/tmp-peaks-gene-strand.bed


echo "8. Prepare bed file with ltr ang genes which will convert to gtf format."
# sort mart-export-v10bed2-mm10.bed
bedtools sort -i $data_dir/mart-export-v102-mm10.bed \
  | uniq >  $tmp_dir/tmp-mart-export-v102-mm10-sorted.bed

# using peaks which are defined to genes chose the intersect peak to gene in +/-30000 range from gene 
# peaks are describtion to nearest gene -> assesment strand with the same way as for ltr
bedtools intersect \
  -a $tmp_dir/tmp-peaks-gene-strand.bed \
  -b <( cat $tmp_dir/tmp-mart-export-v102-mm10-sorted.bed | \
    awk 'BEGIN{OFS = "\t"} {print $1, $2-30000, $3+30000, $4, $5, $6, $2, $3}' | \
    awk 'BEGIN{OFS = "\t"} {if ($6=="-1") {$6="-"} else {$6="+"}}1') \
  -wa -wb -s  2>/dev/null | \
  awk 'BEGIN{OFS = "\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $19, $15, $16, $17}' \
  > $tmp_dir/tmp-peaks2gtf-gene.bed


# bedtools closest \
#   -a $tmp_dir/tmp-peaks-gene-strand.bed \
#   -b $tmp_dir/tmp-mart-export-v102-mm10-sorted.bed \
#   -t first 2>/dev/null \
#   > $tmp_dir/tmp-peaks2gtf-gene.bed

# using peaks which are defined to ltr chose the intersect peak to gene in +/-30000 range from gene
# ltr are describtion to nearest gene
bedtools intersect \
  -a $tmp_dir/tmp-peaks-ltr-strand.bed \
  -b <( cat $tmp_dir/tmp-mart-export-v102-mm10-sorted.bed | \
    awk 'BEGIN{OFS = "\t"} {print $1, $2-30000, $3+30000, $4, $5, $6, $2, $3}' | \
    awk 'BEGIN{OFS = "\t"} {if ($6=="-1") {$6="-"} else {$6="+"}}1') \
  -wa -wb -s  2>/dev/null | \
  awk 'BEGIN{OFS = "\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $19, $15, $16, $17}' \
  > $tmp_dir/tmp-peaks2gtf-ltr.bed

# bedtools closest \
#   -a $tmp_dir/tmp-peaks-ltr-strand.bed \
#   -b $tmp_dir/tmp-mart-export-v102-mm10-sorted.bed \
#   -t first 2>/dev/null \
#   > $tmp_dir/tmp-peaks2gtf-ltr.bed

echo "9. Conect file bed for ltr and gene"
# conect together peaks for ltr and genes
cat $tmp_dir/tmp-peaks2gtf-gene.bed $tmp_dir/tmp-peaks2gtf-ltr.bed | \
  # remove row with GL4 and JH5
  grep -Pv "GL4|JH5" \
  > $tmp_dir/tmp-peaks-annotate.bed


echo "10. sort file bed"
# sort tmp-peaks-annotate.bed
bedtools sort -i $tmp_dir/tmp-peaks-annotate.bed > $data_dir/peaks-annotate-sort.bed

# remove tmp folder
# rm -r $tmp_dir

