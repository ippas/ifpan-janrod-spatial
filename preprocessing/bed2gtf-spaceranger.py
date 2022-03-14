# Script write in python 3.8.2

##########################
# import needed packages #
##########################
import csv
import os
import pandas as pd
import sys
import argparse


def create_gtf_gene(data_frame):
    list_gtf_gene = []

    for number in range(0, len(data_frame)):
        gtf_gene = [data_frame.iloc[number, 0], 'MACS3', "gene", data_frame.iloc[number, 1],
                    data_frame.iloc[number, 2], ".", data_frame.iloc[number, 5], ".",
                    f'gene_id "{data_frame.iloc[number, 3]}"; gene_version "1"; gene_type "protein_coding"; '
                    # f'gene_name "{data_frame.iloc[number, 14]}";']
                    f'gene_name "{data_frame.iloc[number, 3]}";']

        gtf_transcript = [data_frame.iloc[number, 0], 'MACS3', "transcript", data_frame.iloc[number, 1],
                          data_frame.iloc[number, 2], ".", data_frame.iloc[number, 5], ".",
                          f'gene_id "{data_frame.iloc[number, 3]}"; '
                          f'transcript_id "{data_frame.iloc[number, 3]}"; gene_version '
                          f'"1"; gene_type "protein_coding"; '
                          # f'gene_name "{data_frame.iloc[0, 14]}";']
                          f'gene_name "{data_frame.iloc[number, 3]}";']

        gtf_exon = [data_frame.iloc[number, 0], 'MACS3', "exon", data_frame.iloc[number, 1],
                    data_frame.iloc[number, 2], ".", data_frame.iloc[number, 5], ".",
                    f'gene_id "{data_frame.iloc[number, 3]}"; '
                    f'transcript_id "{data_frame.iloc[number, 3]}"; gene_version '
                    f'"1"; gene_type "protein_coding"; '
                    # f'gene_name "{data_frame.iloc[0, 14]}";']
                    f'gene_name "{data_frame.iloc[number, 3]}";']

        list_gtf_gene.append(gtf_gene)
        list_gtf_gene.append(gtf_transcript)
        list_gtf_gene.append(gtf_exon)

    return list_gtf_gene

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest = "inputfile")
parser.add_argument("-o", "--output", dest = "outputfile")

args = parser.parse_args()

input_file = args.inputfile
output_file = args.outputfile

bed_to_gtf = pd.read_csv(input_file, sep='\t')

# print(bed_to_gtf[bed_to_gtf.duplicated(subset=['peak_id', 'gene_name'], keep=False)])
df_gtf_all = pd.DataFrame()
list_gene_names = bed_to_gtf["gene_name"].unique().tolist()

for gene in list_gene_names:
    tmp_gene_bed2gtf = bed_to_gtf[bed_to_gtf.gene_name == gene]
    df_gtf_all = pd.concat([df_gtf_all, pd.DataFrame(create_gtf_gene(data_frame=tmp_gene_bed2gtf))], ignore_index=True)

gtf_head = """##description: evidence-based annotation of the mouse genome (GRCm38), version M23 (Ensembl 98)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2019-09-06
"""

file_gtf = open("tmp_file1.txt", "w")
file_gtf.write(gtf_head)
file_gtf.close()

pd.DataFrame(df_gtf_all).to_csv("tmp_file2.tsv", index=False, header=False, sep="\t", quoting=csv.QUOTE_NONE)

data = data2 = ""

with open('tmp_file1.txt') as fp:
    data = fp.read()

with open('tmp_file2.tsv') as fp:
    data2 = fp.read()

data += data2

with open(output_file, 'w') as fp:
    fp.write(data)

os.remove("tmp_file1.txt")
os.remove("tmp_file2.tsv")
