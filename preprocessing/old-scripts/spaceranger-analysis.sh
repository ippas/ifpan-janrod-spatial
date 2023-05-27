#!/bin/bash

# initialize the variable with the default value
localcores=$(nproc)

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
            --metadata)
                shift
                metadata=$1
                shift
                ;;
            --samples)
                shift
                samples=$1
                shift
                ;;
            --localcores)
                shift
                localcores=$1
                shift
                ;;
            *)
                echo "$1 is not a recognized flag!"                 
                break;
                ;;
    esac
done  

# create variable to grep many samples by grep
grep_samples=$(echo $samples | sed 's/,/\|/g; s/^/"|/; s/$/|"/')

# output_dir="data/spaceranger/"

# check that folder for output exist 
# if not create folder
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi


# cat data/samples-spatial-metadata.tsv | sed 1d | grep 

# loop for run spaceranger
cat $metadata | sed 1d | grep -P $grep_samples | \
while read line; do

  # create variable needed for spaceranger
  sample=$(echo $line | awk '{print $1}')
  slide=$(echo $line | awk '{print $2}')
  area=$(echo $line | awk '{print $3}')
  tif_file=(raw/$sample/*.tif)
  json_file=(raw/$sample/*.json)


  echo ""
  echo $sample
  echo $slide
  echo $area
  echo $tif_file
  echo $json_file
  echo $localcores
  echo ""


  # run spaceranger for samples in /samples-spatial-metadata.tsv
  spaceranger count --id=$sample \
                    --transcriptome=$transcriptome \
                    --fastqs=raw/$sample \
                    --image=$tif_file \
                    --slide=$slide \
                    --area=$area \
                    --loupe-alignment=$json_file \
                    --localcores=$localcores \
                    --sample=$sample
  
  # move folder with output spaceranger to output_dir
  mv $sample $output_dir

done 
