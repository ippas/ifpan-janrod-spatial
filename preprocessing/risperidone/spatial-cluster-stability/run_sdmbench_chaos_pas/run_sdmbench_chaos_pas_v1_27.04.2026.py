#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import sys

SDMBENCH_PATH = "/home/mateusz/projects/ifpan-janrod-spatial/tools/SDMBench/SDMBench"

if SDMBENCH_PATH not in sys.path:
    sys.path.insert(0, SDMBENCH_PATH)

import pandas as pd
import numpy as np
import anndata as ad
from SDMBench import sdmbench


# ============================================================
# Paths
# ============================================================

output_dir = "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone/spatial-cluster-stability"

coords_path = f"{output_dir}/risperidone_half_spatial_coords.tsv"
clusters_path = f"{output_dir}/risperidone_half_clusters.tsv"

out_path = f"{output_dir}/risperidone_half_SDMbench_CHAOS_PAS.tsv"


# ============================================================
# Load input data
# ============================================================

print("\n[1/6] Loading input files...")

coords_df = pd.read_csv(coords_path, sep="\t")
clusters_df = pd.read_csv(clusters_path, sep="\t", index_col=0)

print(f"  Spatial coords file: {coords_path}")
print(f"  Cluster file:        {clusters_path}")
print(f"  Spatial coords:      {coords_df.shape[0]} spots x {coords_df.shape[1]} columns")
print(f"  Cluster matrix:      {clusters_df.shape[0]} spots x {clusters_df.shape[1]} resolutions")


# ============================================================
# Check required columns
# ============================================================

print("\n[2/6] Checking required columns...")

required_cols = {"barcode", "imagerow", "imagecol"}
missing_cols = required_cols - set(coords_df.columns)

if missing_cols:
    raise ValueError(f"Missing required columns in coords file: {missing_cols}")

print("  Required columns found: barcode, imagerow, imagecol")


# ============================================================
# Match coords and cluster labels
# ============================================================

print("\n[3/6] Matching spatial coords with cluster labels by barcode...")

df = coords_df.merge(
    clusters_df,
    left_on="barcode",
    right_index=True,
    how="inner"
)

print(f"  Matched spots: {df.shape[0]}")

if df.shape[0] == 0:
    raise ValueError("No matching barcodes between coords and cluster files.")

unmatched_coords = coords_df.shape[0] - df.shape[0]
unmatched_clusters = clusters_df.shape[0] - df.shape[0]

if unmatched_coords > 0:
    print(f"  WARNING: {unmatched_coords} spots from coords file were not matched.")

if unmatched_clusters > 0:
    print(f"  WARNING: {unmatched_clusters} spots from cluster file were not matched.")


# ============================================================
# Create minimal AnnData object for SDMBench
# ============================================================

print("\n[4/6] Creating minimal AnnData object...")

adata = ad.AnnData(X=np.zeros((df.shape[0], 1)))
adata.obs_names = df["barcode"].astype(str).values

# IMPORTANT:
# CHAOS/PAS use spatial tissue coordinates, not PCA/UMAP.
adata.obsm["spatial"] = df[["imagerow", "imagecol"]].to_numpy()

cluster_cols = list(clusters_df.columns)

for col in cluster_cols:
    adata.obs[col] = df[col].astype(str).values

print(f"  AnnData object created with {adata.n_obs} spots")
print("  Spatial coordinates stored in: adata.obsm['spatial']")
print(f"  Added {len(cluster_cols)} cluster columns")


# ============================================================
# Compute CHAOS and PAS
# ============================================================

print("\n[5/6] Computing CHAOS and PAS...")

results = []

for i, col in enumerate(cluster_cols, start=1):
    resolution = col.replace("cluster_resolution_", "")
    n_clusters = adata.obs[col].nunique()

    print(
        f"  [{i}/{len(cluster_cols)}] "
        f"resolution={resolution} | "
        f"column={col} | "
        f"n_clusters={n_clusters}"
    )

    chaos = sdmbench.compute_CHAOS(
        adata,
        pred_key=col,
        spatial_key="spatial"
    )

    pas = sdmbench.compute_PAS(
        adata,
        pred_key=col,
        spatial_key="spatial"
    )

    results.append({
        "resolution": resolution,
        "cluster_column": col,
        "n_spots": adata.n_obs,
        "n_clusters": n_clusters,
        "CHAOS": chaos,
        "PAS": pas
    })


# ============================================================
# Save results
# ============================================================

print("\n[6/6] Saving results...")

results_df = pd.DataFrame(results)

results_df.to_csv(
    out_path,
    sep="\t",
    index=False
)

print("\nDONE")
print(results_df)
print(f"\nSaved to:\n{out_path}")
