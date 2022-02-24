#!/bin/bash

# assigning arguments from flags to varibles
while test $# -gt 0; do
    case "$1" in
        --transcriptome)
            shift
                transcriptome=$1
                shift
                ;;
            --output-dir)
                shift
                output_dir=$1
                shift
                ;;
            *)
                echo "$1 is not a recognized flag!"                 
                break;
                ;;
    esac
done  


# output_dir="data/spaceranger/"

# check that folder for output exist 
# if not create folder
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

# loop for run spaceranger
cat data/samples-spatial-metadata.tsv | sed 1d | \
while read line; do

  # create variable needed for spaceranger
  sample=$(echo $line | awk '{print $1}')
  slide=$(echo $line | awk '{print $2}')
  area=$(echo $line | awk '{print $3}')
  tif_file=(raw/tif/$sample.tif)
  json_file=(raw/json/$sample"-"$slide"-"$area.json)
  localcores=16

#   echo ""
#   echo $sample
#   echo $slide
#   echo $area
#   echo $tif_file
#   echo $json_file
  
  # run spaceranger for samples in /samples-spatial-metadata.tsv
  spaceranger count --id=$sample \
                    --transcriptome=$transcriptome \
                    --fastqs=raw/fastq \
                    --image=$tif_file \
                    --slide=$slide \
                    --area=$area \
                    --loupe-alignment=$json_file \
                    --localcores=$localcores \
                    --sample=$sample
  
  # move folder with output spaceranger to output_dir
  mv $sample $output_dir

done 