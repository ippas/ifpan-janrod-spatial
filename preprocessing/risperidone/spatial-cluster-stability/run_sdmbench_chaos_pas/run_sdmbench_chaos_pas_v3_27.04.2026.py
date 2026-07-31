#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import sys
import numpy as np
import pandas as pd

SDMBENCH_PATH = "/home/mateusz/projects/ifpan-janrod-spatial/tools/SDMBench/SDMBench"

if SDMBENCH_PATH not in sys.path:
    sys.path.insert(0, SDMBENCH_PATH)

from SDMBench import sdmbench
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


# ============================================================
# Paths
# ============================================================

output_dir = "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone/spatial-cluster-stability"

coords_path = f"{output_dir}/risperidone_half_spatial_coords.tsv"
clusters_path = f"{output_dir}/risperidone_half_clusters.tsv"
pca_path = f"{output_dir}/risperidone_half_pca_embeddings_1_30.tsv"

out_path = f"{output_dir}/risperidone_half_SDMbench_metrics.tsv"


# ============================================================
# Monkey patch for broken SDMBench private functions
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
        chaos_values.append(np.mean(distances[:, 1]))

    return np.mean(chaos_values) if len(chaos_values) > 0 else np.nan


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
        neighbor_idx = indices[i, 1:]
        neighbor_labels = pred[neighbor_idx]

        values, counts = np.unique(neighbor_labels, return_counts=True)
        majority_label = values[np.argmax(counts)]

        abnormal.append(pred[i] != majority_label)

    return np.mean(abnormal)


# SDMBench ma bug: compute_CHAOS() woła _compute_CHAOS, którego nie widzi.
# Dlatego dokładamy tę funkcję do namespace pakietu.
sdmbench._compute_CHAOS = _compute_CHAOS_patch

# Na wszelki wypadek także PAS, jeśli u Ciebie też byłby popsuty.
sdmbench._compute_PAS = _compute_PAS_patch


# ============================================================
# Load data
# ============================================================

print("\n[1/5] Loading input files...")

coords_df = pd.read_csv(coords_path, sep="\t")
clusters_df = pd.read_csv(clusters_path, sep="\t", index_col=0)

print(f"  Spatial coords: {coords_df.shape}")
print(f"  Clusters:       {clusters_df.shape}")


# ============================================================
# Match data
# ============================================================

print("\n[2/5] Matching coords with clusters...")

df = coords_df.merge(
    clusters_df,
    left_on="barcode",
    right_index=True,
    how="inner"
)

print(f"  Matched spots: {df.shape[0]}")

if df.shape[0] == 0:
    raise ValueError("No matching barcodes.")


# ============================================================
# Build minimal AnnData for SDMBench
# ============================================================

print("\n[3/5] Creating AnnData for SDMBench...")

import anndata as ad

adata = ad.AnnData(X=np.zeros((df.shape[0], 1)))
adata.obs_names = df["barcode"].astype(str).values
adata.obsm["spatial"] = df[["imagerow", "imagecol"]].to_numpy()

cluster_cols = list(clusters_df.columns)

for col in cluster_cols:
    adata.obs[col] = df[col].astype(str).values

print(f"  AnnData: {adata.n_obs} spots")
print(f"  Cluster columns: {len(cluster_cols)}")


# ============================================================
# Compute SDMBench metrics
# ============================================================

print("\n[4/5] Computing SDMBench metrics...")

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
results_df.to_csv(out_path, sep="\t", index=False)

print("\nDONE")
print(results_df)
print(f"\nSaved to:\n{out_path}")
