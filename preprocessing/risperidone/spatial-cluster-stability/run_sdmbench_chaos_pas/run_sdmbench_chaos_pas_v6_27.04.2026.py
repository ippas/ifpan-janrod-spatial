#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import sys
import numpy as np
import pandas as pd
import anndata as ad

from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score


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

out_path = f"{output_dir}/risperidone_half_SDMbench_original_metrics.tsv"


# ============================================================
# Original SDMBench helper functions
# ============================================================

def _compute_CHAOS(clusterlabel, location):
    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = StandardScaler().fit_transform(location)

    clusterlabel_unique = np.unique(clusterlabel)
    dist_val = np.zeros(len(clusterlabel_unique))
    count = 0

    for k in clusterlabel_unique:
        location_cluster = matched_location[clusterlabel == k, :]

        if len(location_cluster) <= 2:
            continue

        n_location_cluster = len(location_cluster)
        results = [fx_1NN(i, location_cluster) for i in range(n_location_cluster)]
        dist_val[count] = np.sum(results)
        count = count + 1

    return np.sum(dist_val) / len(clusterlabel)


def fx_1NN(i, location_in):
    location_in = np.array(location_in)
    dist_array = distance_matrix(location_in[i, :][None, :], location_in)[0, :]
    dist_array[i] = np.inf
    return np.min(dist_array)


def fx_kNN(i, location_in, k, cluster_in):
    location_in = np.array(location_in)
    cluster_in = np.array(cluster_in)

    dist_array = distance_matrix(location_in[i, :][None, :], location_in)[0, :]
    dist_array[i] = np.inf
    ind = np.argsort(dist_array)[:k]

    cluster_use = np.array(cluster_in)

    if np.sum(cluster_use[ind] != cluster_in[i]) > (k / 2):
        return 1
    else:
        return 0


def _compute_PAS(clusterlabel, location):
    clusterlabel = np.array(clusterlabel)
    location = np.array(location)
    matched_location = location

    results = [
        fx_kNN(i, matched_location, k=10, cluster_in=clusterlabel)
        for i in range(matched_location.shape[0])
    ]

    return np.sum(results) / len(clusterlabel)


# ============================================================
# Patch SDMBench globals
# ============================================================

sdmbench.compute_CHAOS.__globals__["_compute_CHAOS"] = _compute_CHAOS
sdmbench.compute_CHAOS.__globals__["fx_1NN"] = fx_1NN

sdmbench.compute_PAS.__globals__["_compute_PAS"] = _compute_PAS
sdmbench.compute_PAS.__globals__["fx_kNN"] = fx_kNN


# ============================================================
# Load data
# ============================================================

print("\n[1/5] Loading input files...")

coords_df = pd.read_csv(coords_path, sep="\t")
clusters_df = pd.read_csv(clusters_path, sep="\t", index_col=0)

print(f"  Spatial coords: {coords_df.shape[0]} spots x {coords_df.shape[1]} columns")
print(f"  Clusters:       {clusters_df.shape[0]} spots x {clusters_df.shape[1]} resolutions")


# ============================================================
# Check data
# ============================================================

print("\n[2/5] Checking required columns...")

required_cols = {"barcode", "imagerow", "imagecol"}
missing_cols = required_cols - set(coords_df.columns)

if missing_cols:
    raise ValueError(f"Missing required columns in coords file: {missing_cols}")

print("  Required columns found: barcode, imagerow, imagecol")


# ============================================================
# Match coords and clusters
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
    raise ValueError("No matching barcodes between coords and clusters.")


# ============================================================
# Create AnnData
# ============================================================

print("\n[4/5] Creating AnnData and computing SDMBench metrics...")

adata = ad.AnnData(X=np.zeros((df.shape[0], 1)))
adata.obs_names = df["barcode"].astype(str).values
adata.obsm["spatial"] = df[["imagerow", "imagecol"]].to_numpy()

cluster_cols = list(clusters_df.columns)

for col in cluster_cols:
    adata.obs[col] = df[col].astype(str).values

print(f"  AnnData spots:   {adata.n_obs}")
print(f"  Spatial key:     adata.obsm['spatial']")
print(f"  Cluster columns: {len(cluster_cols)}")


results = []

for i, col in enumerate(cluster_cols, start=1):
    resolution = col.replace("cluster_resolution_", "")
    n_clusters = adata.obs[col].nunique()

    print(f"  [{i}/{len(cluster_cols)}] resolution={resolution} | n_clusters={n_clusters}")

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

    asw = sdmbench.compute_ASW(
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
        "SDMBench_PAS": pas,
        "SDMBench_ASW_spatial": asw
    })


# ============================================================
# Save results
# ============================================================

print("\n[5/5] Saving results...")

results_df = pd.DataFrame(results)
results_df.to_csv(out_path, sep="\t", index=False)

print("\nDONE")
print(results_df)
print(f"\nSaved to:\n{out_path}")
