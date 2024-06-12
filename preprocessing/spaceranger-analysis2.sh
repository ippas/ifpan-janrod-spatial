#!/bin/bash

# initialize the variable with the default value
localcores=$(nproc)
json_sufix=""

echo "hello"

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
            --json-sufix)
                shift
                json_sufix=$1
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

# if json_sufix is not empty add .json to json_sufix
if [ ! -z "$json_sufix" ]; then
  json_sufix="-$json_sufix"
fi

# create variable to grep many samples by grep
grep_samples=$(echo $samples | sed 's/,/\\t|/g; s/^/"|/; s/$/\\t|"/')

# output_dir="data/spaceranger/"

# check that folder for output exist 
# if not create folder
if [ ! -d "$output_dir" ]; then
  mkdir $output_dir
fi

echo "hello"

# cat data/samples-spatial-metadata.tsv | sed 1d | grep 

# loop for run spaceranger
cat $metadata | sed 1d | grep -P $grep_samples | \
while read line; do

  # create variable needed for spaceranger
  sample=$(echo $line | awk '{print $1}')
  slide=$(echo $line | awk '{print $2}')
  area=$(echo $line | awk '{print $3}')
  tif_file=(raw/$sample/*.tif)
  json_file=(raw/$sample/$sample-$slide-$area$json_sufix.json)
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
