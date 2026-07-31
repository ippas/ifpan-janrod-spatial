#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import sys
import numpy as np
import pandas as pd
import anndata as ad

from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


# ============================================================
# Add SDMBench path
# ============================================================

SDMBENCH_PATH = "/home/mateusz/projects/ifpan-janrod-spatial/tools/SDMBench/SDMBench"

if SDMBENCH_PATH not in sys.path:
    sys.path.insert(0, SDMBENCH_PATH)

from SDMBench import sdmbench


# ============================================================
# Paths
# ============================================================

output_dir = "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone/spatial-cluster-stability"

coords_path = f"{output_dir}/risperidone_half_spatial_coords.tsv"
clusters_path = f"{output_dir}/risperidone_half_clusters.tsv"

out_path = f"{output_dir}/risperidone_half_SDMbench_CHAOS_PAS.tsv"


# ============================================================
# Patch broken SDMBench private functions
# ============================================================

def _compute_CHAOS_patch(pred, spatial):
    spatial = StandardScaler().fit_transform(np.asarray(spatial))
    pred = np.asarray(pred).astype(str)

    chaos_values = []

    for cluster in np.unique(pred):
        cluster_spatial = spatial[pred == cluster]

        if cluster_spatial.shape[0] < 2:
            continue

        nbrs = NearestNeighbors(n_neighbors=2)
        nbrs.fit(cluster_spatial)

        distances, _ = nbrs.kneighbors(cluster_spatial)

        # distance to nearest neighbor from the same cluster
        chaos_values.append(np.mean(distances[:, 1]))

    if len(chaos_values) == 0:
        return np.nan

    return np.mean(chaos_values)


def _compute_PAS_patch(pred, spatial, k=10):
    spatial = StandardScaler().fit_transform(np.asarray(spatial))
    pred = np.asarray(pred).astype(str)

    if spatial.shape[0] <= k:
        return np.nan

    nbrs = NearestNeighbors(n_neighbors=k + 1)
    nbrs.fit(spatial)

    _, indices = nbrs.kneighbors(spatial)

    abnormal = []

    for i in range(spatial.shape[0]):
        neighbor_idx = indices[i, 1:]  # remove self-neighbor
        neighbor_labels = pred[neighbor_idx]

        values, counts = np.unique(neighbor_labels, return_counts=True)
        majority_label = values[np.argmax(counts)]

        abnormal.append(pred[i] != majority_label)

    return np.mean(abnormal)


# IMPORTANT:
# SDMBench compute_CHAOS() searches for _compute_CHAOS
# in the function global namespace.
sdmbench.compute_CHAOS.__globals__["_compute_CHAOS"] = _compute_CHAOS_patch
sdmbench.compute_PAS.__globals__["_compute_PAS"] = _compute_PAS_patch


# ============================================================
# Load data
# ============================================================

print("\n[1/5] Loading input files...")

coords_df = pd.read_csv(coords_path, sep="\t")
clusters_df = pd.read_csv(clusters_path, sep="\t", index_col=0)

print(f"  Spatial coords file: {coords_path}")
print(f"  Cluster file:        {clusters_path}")
print(f"  Spatial coords:      {coords_df.shape[0]} spots x {coords_df.shape[1]} columns")
print(f"  Cluster matrix:      {clusters_df.shape[0]} spots x {clusters_df.shape[1]} resolutions")


# ============================================================
# Check columns
# ============================================================

print("\n[2/5] Checking required columns...")

required_cols = {"barcode", "imagerow", "imagecol"}
missing_cols = required_cols - set(coords_df.columns)

if missing_cols:
    raise ValueError(f"Missing required columns in coords file: {missing_cols}")

print("  Required columns found: barcode, imagerow, imagecol")


# ============================================================
# Match coords with clusters
# ============================================================

print("\n[3/5] Matching coords with cluster labels by barcode...")

df = coords_df.merge(
    clusters_df,
    left_on="barcode",
    right_index=True,
    how="inner"
)

print(f"  Matched spots: {df.shape[0]}")

if df.shape[0] == 0:
    raise ValueError("No matching barcodes between coords and cluster files.")

if df.shape[0] < coords_df.shape[0]:
    print(f"  WARNING: {coords_df.shape[0] - df.shape[0]} coords spots were not matched.")

if df.shape[0] < clusters_df.shape[0]:
    print(f"  WARNING: {clusters_df.shape[0] - df.shape[0]} cluster spots were not matched.")


# ============================================================
# Create AnnData for SDMBench
# ============================================================

print("\n[4/5] Creating AnnData and computing SDMBench metrics...")

adata = ad.AnnData(X=np.zeros((df.shape[0], 1)))
adata.obs_names = df["barcode"].astype(str).values

# CHAOS/PAS use tissue coordinates, not PCA or UMAP
adata.obsm["spatial"] = df[["imagerow", "imagecol"]].to_numpy()

cluster_cols = list(clusters_df.columns)

for col in cluster_cols:
    adata.obs[col] = df[col].astype(str).values

print(f"  AnnData:          {adata.n_obs} spots")
print(f"  Spatial key:      adata.obsm['spatial']")
print(f"  Cluster columns:  {len(cluster_cols)}")


results = []

for i, col in enumerate(cluster_cols, start=1):
    resolution = col.replace("cluster_resolution_", "")
    n_clusters = adata.obs[col].nunique()

    print(
        f"  [{i}/{len(cluster_cols)}] "
        f"resolution={resolution} | "
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
        "SDMBench_CHAOS": chaos,
        "SDMBench_PAS": pas
    })


# ============================================================
# Save results
# ============================================================

print("\n[5/5] Saving results...")

results_df = pd.DataFrame(results)

results_df.to_csv(
    out_path,
    sep="\t",
    index=False
)

print("\nDONE")
print(results_df)
print(f"\nSaved to:\n{out_path}")
