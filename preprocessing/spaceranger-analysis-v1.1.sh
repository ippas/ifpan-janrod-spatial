#!/bin/bash

# initialize variables with default values
localcores=$(nproc)
localmem=""
json_sufix=""
input_directory="raw"

# assigning arguments from flags to variables
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
        --localmem)
            shift
                localmem=$1
                shift
                ;;
        --input-directory)
            shift
                input_directory=$1
                shift
                ;;
        *)
            echo "$1 is not a recognized flag!"
            exit 1
            ;;
    esac
done

# check required arguments
if [ -z "$transcriptome" ]; then
  echo "Error: --transcriptome is required"
  exit 1
fi

if [ -z "$output_dir" ]; then
  echo "Error: --output-dir is required"
  exit 1
fi

if [ -z "$metadata" ]; then
  echo "Error: --metadata is required"
  exit 1
fi

if [ -z "$samples" ]; then
  echo "Error: --samples is required"
  exit 1
fi

# if json_sufix is not empty add "-" before suffix
if [ ! -z "$json_sufix" ]; then
  json_sufix="-$json_sufix"
fi

# create variable to grep many samples by grep
grep_samples=$(echo "$samples" | sed 's/,/\\t|/g; s/^/"|/; s/$/\\t|"/')

# check that folder for output exist, if not create folder
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
fi

# loop for run spaceranger
cat "$metadata" | sed 1d | grep -P "$grep_samples" | \
while read line; do

  # create variables needed for spaceranger
  sample=$(echo "$line" | awk '{print $1}')
  slide=$(echo "$line" | awk '{print $2}')
  area=$(echo "$line" | awk '{print $3}')

  tif_file=($input_directory/$sample/*.tif)
  json_file=($input_directory/$sample/$sample-$slide-$area$json_sufix.json)

  echo ""
  echo "$sample"
  echo "$slide"
  echo "$area"
  echo "$tif_file"
  echo "$json_file"
  echo "$localcores"
  echo "$localmem"
  echo ""

  # check input files
  if [ ! -f "$json_file" ]; then
    echo "Error: JSON file not found: $json_file"
    exit 1
  fi

  if [ ! -f "$tif_file" ]; then
    echo "Error: TIFF file not found: $tif_file"
    exit 1
  fi

  # run spaceranger for samples in metadata
  if [ ! -z "$localmem" ]; then
    spaceranger count --id="$sample" \
                      --transcriptome="$transcriptome" \
                      --fastqs="$input_directory/$sample" \
                      --image="$tif_file" \
                      --slide="$slide" \
                      --area="$area" \
                      --loupe-alignment="$json_file" \
                      --localcores="$localcores" \
                      --localmem="$localmem" \
                      --sample="$sample"
  else
    spaceranger count --id="$sample" \
                      --transcriptome="$transcriptome" \
                      --fastqs="$input_directory/$sample" \
                      --image="$tif_file" \
                      --slide="$slide" \
                      --area="$area" \
                      --loupe-alignment="$json_file" \
                      --localcores="$localcores" \
                      --sample="$sample"
  fi

  # move folder with output spaceranger to output_dir only if it exists
  if [ -d "$sample" ]; then
    mv "$sample" "$output_dir"
  else
    echo "Warning: output directory '$sample' was not created, skipping mv"
    exit 1
  fi

done
