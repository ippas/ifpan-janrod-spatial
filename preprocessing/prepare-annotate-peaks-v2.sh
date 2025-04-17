#!/bin/bash

# Default values
prefix="tmp"
bam_file=""

# Parse command-line arguments
while test $# -gt 0; do
    case "$1" in
        --gene-bed)
            shift
            gene_bed=$1
            shift
            ;;
        --ltr-bed)
            shift
            ltr_bed=$1
            shift
            ;;
        --macs-peaks)
            shift
            macs_peaks=$1
            shift
            ;;
        --data-dir)
            shift
            data_dir=$1
            shift
            ;; 
        --number-threads)
            shift
            number_threads=$1
            shift
            ;;
        --prefix)
            shift
            prefix=$1
            shift
            ;;
        --bam)
            shift
            bam_file=$1
            shift
            ;;
        *)
            echo "$1 is not a recognized flag!"                 
            break;
            ;;
    esac
done  

# Set default BAM file if not provided
if [ -z "$bam_file" ]; then
    bam_file="$data_dir/merged-samples.bam"
fi

# Create temporary directory with prefix
tmp_dir=$data_dir/${prefix}_tmp

# Clean and recreate temp directory
if [ -d "$tmp_dir" ]; then
    echo "Removing old temporary directory: $tmp_dir"
    rm -rf "$tmp_dir"
fi
mkdir -p "$tmp_dir"

echo "1. Prepare BED files for LTR and gene peaks"
bedtools intersect -u -a $macs_peaks -b $ltr_bed 2>/dev/null | \
  awk '{print $0"\tltr"}' > $tmp_dir/${prefix}-peaks-ltr.bed

bedtools intersect -v -a $macs_peaks -b $ltr_bed 2>/dev/null | \
  awk '{print $0"\tgene"}' > $tmp_dir/${prefix}-peaks-gene.bed

echo "2. Extract reads by strand from BAM"
samtools view -@ $number_threads -b -f 16 $bam_file > $tmp_dir/${prefix}-merged-samples-minus.bam
samtools view -@ $number_threads -b -F 16 $bam_file > $tmp_dir/${prefix}-merged-samples-plus.bam

echo "3. Index the strand-specific BAM files"
samtools index -@ $number_threads $tmp_dir/${prefix}-merged-samples-minus.bam
samtools index -@ $number_threads $tmp_dir/${prefix}-merged-samples-plus.bam

echo "4. Compute coverage for LTR peaks"
samtools bedcov $tmp_dir/${prefix}-peaks-ltr.bed $tmp_dir/${prefix}-merged-samples-minus.bam > $tmp_dir/${prefix}-peaks-ltr-coverage-minus.bed
samtools bedcov $tmp_dir/${prefix}-peaks-ltr.bed $tmp_dir/${prefix}-merged-samples-plus.bam > $tmp_dir/${prefix}-peaks-ltr-coverage-plus.bed

echo "5. Compute coverage for gene peaks"
samtools bedcov $tmp_dir/${prefix}-peaks-gene.bed $tmp_dir/${prefix}-merged-samples-minus.bam > $tmp_dir/${prefix}-peaks-gene-coverage-minus.bed
samtools bedcov $tmp_dir/${prefix}-peaks-gene.bed $tmp_dir/${prefix}-merged-samples-plus.bam > $tmp_dir/${prefix}-peaks-gene-coverage-plus.bed

echo "6. Assign strand to LTR peaks"
paste $tmp_dir/${prefix}-peaks-ltr-coverage-plus.bed $tmp_dir/${prefix}-peaks-ltr-coverage-minus.bed |
    awk '{print $0"\t+"$12-$24}' |
    awk 'BEGIN{FS=OFS="\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
    awk -F"\t" '{OFS=FS}{ $6=$25 ; print }' |
    cut -f1-11 > $tmp_dir/${prefix}-peaks-ltr-strand.bed

echo "7. Assign strand to gene peaks"
paste $tmp_dir/${prefix}-peaks-gene-coverage-plus.bed $tmp_dir/${prefix}-peaks-gene-coverage-minus.bed |
    awk '{print $0"\t+"$12-$24}' |
    awk 'BEGIN{FS=OFS="\t"} {gsub(/+-[0-9]*/, "-" $3)} 1 {gsub(/+[0-9]*/, "+" $3)} 1' |
    awk -F"\t" '{OFS=FS}{ $6=$25 ; print }' |
    cut -f1-11  > $tmp_dir/${prefix}-peaks-gene-strand.bed

echo "8. Annotate peaks with gene information"
bedtools sort -i $gene_bed | uniq >  $tmp_dir/${prefix}-mart-export-sorted.bed

bedtools intersect \
  -a $tmp_dir/${prefix}-peaks-gene-strand.bed \
  -b <(cat $tmp_dir/${prefix}-mart-export-sorted.bed | \
    awk 'BEGIN{OFS="\t"} {print $1, $2-30000, $3+30000, $4, $5, $6, $2, $3}' | \
    awk 'BEGIN{OFS="\t"} {if ($6=="-1") {$6="-"} else {$6="+"}}1') \
  -wa -wb -s 2>/dev/null | \
  awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $19, $15, $16, $17}' \
  > $tmp_dir/${prefix}-peaks2gtf-gene.bed

bedtools intersect \
  -a $tmp_dir/${prefix}-peaks-ltr-strand.bed \
  -b <(cat $tmp_dir/${prefix}-mart-export-sorted.bed | \
    awk 'BEGIN{OFS="\t"} {print $1, $2-30000, $3+30000, $4, $5, $6, $2, $3}' | \
    awk 'BEGIN{OFS="\t"} {if ($6=="-1") {$6="-"} else {$6="+"}}1') \
  -wa -wb -s 2>/dev/null | \
  awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $18, $19, $15, $16, $17}' \
  > $tmp_dir/${prefix}-peaks2gtf-ltr.bed

echo "9. Merge annotated LTR and gene peaks"
cat $tmp_dir/${prefix}-peaks2gtf-gene.bed $tmp_dir/${prefix}-peaks2gtf-ltr.bed | \
  grep -Pv "GL4|JH5" \
  > $tmp_dir/${prefix}-peaks-annotate.bed

echo "10. Sort and save final annotation file"
bedtools sort -i $tmp_dir/${prefix}-peaks-annotate.bed > $data_dir/${prefix}-peaks-annotate-sort.bed

# Optional: remove temporary files
# rm -r "$tmp_dir"
