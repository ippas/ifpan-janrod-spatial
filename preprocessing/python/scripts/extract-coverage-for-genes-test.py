# %%
import pysam
import pandas as pd
from pybiomart import Dataset, Server
import os
import sys
import numpy as np
import re
import seaborn as sns
import matplotlib.pyplot as plt


os.chdir("/home/mateusz/projects/ifpan-janrod-spatial/")

sys.path.append("preprocessing/python/")


from src.mouse_annotations_downloader import get_mouse_annotations
from src.workflow_utils import find_all_bam_files
from src.bam_coverage_utils import get_per_base_coverage 
from src.bam_coverage_utils import extract_max_coverage_for_genes
from src.plotting_functions import plot_boxplot_with_jitter

# print current working directory
print(os.getcwd())

# read metadata from data/risperidone-3q29/metadata-3q29-ris.tsv 
metadata_3q29_df = pd.read_csv(
    "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone-3q29/metadata-3q29-ris.tsv",
    sep="\t",
    header=0
)

metadata_3q29_df.columns = metadata_3q29_df.columns.str.lower()

metadata_3q29_df["label"] = (
    metadata_3q29_df["sample_id"].astype(str) + "-" +
    metadata_3q29_df["mouse_genotype"].astype(str) + "-" +
    metadata_3q29_df["treatment"].astype(str)
)

metadata_3q29_df.head()


mouse_genes_3q29 = [
    "Bdh1",
    "Dlg1",
    "Meltf",
    "Pigz",
    "Ncbp2",
    "Senp5",
    "Pak2",
    "Pigx",
    "Cep19",
    "Nrros",
    "Fbxo45",
    "Wdr53",
    "Rnf168",
    "Ubxn7",
    "Tm4sf19",
    "Tctex1d2",
    "Pcyt1a",
    "Slc51a",
    "Zdhhc19",
    "Tfrc",
    "Smco1"
]


# genes_mouse_michkor_indicate_T_S1_Cx43DGE = [
mouse_genes_3q29 = [
    "Slc1a2", "Kcng1", "Glul", "Gse1", "Prrx1", "Sdc3", "Ptn", "Acsbg1", "Chrdl1", "Kif13a", "Msantd4", "C1orf61",
    "Nkain4", "Slc13a5", "Il17rb", "Gnao1", "Dmd", "Sox10", "Ophn1", "Shank2", "Itm2c", "Slco1c1", "Slc4a4",
    "Pdgfrb", "Slc25a18", "Megf11", "Nt5e", "Gprc5b", "Elovl2", "Mmd2", "C3orf70", "Tom1l2", "Ncan", "Agbl4",
    "Atp1b2", "Epn2", "Ac009879.3", "Grin2c", "Bcan", "Etnppl", "Smtn", "Prss35", "F3", "Ankrd18b", "Oaf", "Epas1",
    "Cables1", "Dio2", "Mpdz", "Atp13a4", "Lrp4", "Mfge8", "Top3a", "Scd5", "Lrrc8a", "Ednrb", "Lcnl1", "Kat2b",
    "Kiaa1841", "Aldoc", "Hif3a", "Caskin2", "Ptgds", "Atp13a5", "Rhob", "Slc25a26", "Plpp3"
]

mouse_genes_3q29 = [
    "Fkbp5", "Htra1", "Mrc1", "Aldoc", "Fmod", "Gpd1", "Fmo2", "Map3k6", "Ndrg2", "Cyp1b1", "Hif3a", "Hspb1", "Prg4",
    "Ankdd1b", "Hopx", "Slc26a2", "Tnfrsf11a", "Flnc", "Apod", "Pabpc1l", "Cryab", "Igf2", "Pdk4", "S100a11", "Zbtb16",
    "Mt3", "Il1r1", "Wfdc1", "Atp1a2", "Slc26a7", "Cd163", "Omg", "Ogn", "Lmod1", "Rtkn", "Aebp1", "F13a1", "Tgfa",
    "Ca14", "Sult1a2", "Tpm1", "Ctxn3", "Etnppl", "Fstl1", "Bmp5", "Anln", "Prnp", "Lama1", "Gpatch4", "Lgals1"
]

mouse_gene_annotations_v102 = get_mouse_annotations(biomart_version="v102")

# filter gene for Mfi2
mouse_gene_annotations_v102[
    mouse_gene_annotations_v102['gene_name'].isin(mouse_genes_3q29)
]

bam_dir = "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone-3q29/spaceranger-raw/tmp/"
bam_file_list = find_all_bam_files(bam_dir)

max_coverage_3q29del_genes_df = extract_max_coverage_for_genes(
    genes=mouse_genes_3q29,
    files=bam_file_list,
    gene_annotation=mouse_gene_annotations_v102
)

mouse_genes_ref = [
    "Gapdh",
    "Pgk1",
    # "Actb",
    # "Ppia",
    # "Rplp0",
    "Tbp",
    "Atp5pb"
]

max_coverage_ref_genes_df = extract_max_coverage_for_genes(
    genes=mouse_genes_ref,
    files=bam_file_list,
    gene_annotation=mouse_gene_annotations_v102
)

max_coverage_ref_genes_df.head()


preprocessing_3q29del_genes_coverage_df = pd.merge(
    max_coverage_3q29del_genes_df,
    metadata_3q29_df,
    on="sample_id",
    how="left"
)

wide_max_coverage_ref_genes_df = max_coverage_ref_genes_df.pivot(
    index="sample_id",
    columns="gene_symbol",
    values="max_coverage"
).reset_index()

wide_max_coverage_ref_genes_df.columns.name = None

preprocessing_3q29del_genes_coverage_df = pd.merge(
    preprocessing_3q29del_genes_coverage_df,
    wide_max_coverage_ref_genes_df,
    on="sample_id",
    how="left"  # lub "inner", jeśli chcesz tylko dopasowane rekordy
)

for ref in mouse_genes_ref:
    col_name = f"norm_by_{ref.lower()}"
    preprocessing_3q29del_genes_coverage_df[col_name] = (
        preprocessing_3q29del_genes_coverage_df["max_coverage"] /
        preprocessing_3q29del_genes_coverage_df[ref]
    )


preprocessing_3q29del_genes_coverage_df.columns



# mouse_genes_ref = [
#     "Gapdh",
#     "Pgk1",
#     "Actb",
#     "Ppia",
#     "Rplp0"
# ]

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(12, 28))
sort_by_options = ["sample", "mean", "median", "IQR", "max"]

for i, sort_method in enumerate(sort_by_options):
    plot_boxplot_with_jitter(
        data=preprocessing_3q29del_genes_coverage_df,
        x="label",
        y="norm_by_gapdh",
        jitter_strength=0.2,
        point_size=3,
        sort_by=sort_method,
        top_n_labels=5,
        show_top_genes=True,
        ax=axes[i]
    )

# fig.suptitle("Comparison of sort methods for 3q29del gene coverage", fontsize=16, y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

print("end")


