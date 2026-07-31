#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


# ============================================================
# Paths
# ============================================================

output_dir = "/home/mateusz/projects/ifpan-janrod-spatial/data/risperidone/spatial-cluster-stability"

coords_path = f"{output_dir}/risperidone_half_spatial_coords.tsv"
clusters_path = f"{output_dir}/risperidone_half_clusters.tsv"

out_path = f"{output_dir}/risperidone_half_local_CHAOS_PAS.tsv"


# ============================================================
# Local CHAOS / PAS functions
# ============================================================

def compute_CHAOS_local(coords, labels):
    """
    CHAOS:
    mean nearest-neighbor distance within the same spatial cluster.
    Lower CHAOS = more spatially compact clusters.
    """
    coords = StandardScaler().fit_transform(coords)
    labels = np.asarray(labels).astype(str)

    chaos_values = []

    for cluster in np.unique(labels):
        cluster_coords = coords[labels == cluster]

        if cluster_coords.shape[0] < 2:
            continue

        nbrs = NearestNeighbors(n_neighbors=2)
        nbrs.fit(cluster_coords)

        distances, _ = nbrs.kneighbors(cluster_coords)

        # distances[:, 0] = distance to itself
        # distances[:, 1] = nearest neighbor from the same cluster
        chaos_values.append(np.mean(distances[:, 1]))

    if len(chaos_values) == 0:
        return np.nan

    return np.mean(chaos_values)


def compute_PAS_local(coords, labels, k=10):
    """
    PAS:
    fraction of spots whose cluster label differs from
    the majority label among their spatial nearest neighbors.
    Lower PAS = better local spatial coherence.
    """
    coords = StandardScaler().fit_transform(coords)
    labels = np.asarray(labels).astype(str)

    n_spots = coords.shape[0]

    if n_spots <= k:
        return np.nan

    nbrs = NearestNeighbors(n_neighbors=k + 1)
    nbrs.fit(coords)

    _, indices = nbrs.kneighbors(coords)

    abnormal = []

    for i in range(n_spots):
        # remove self-neighbor
        neighbor_idx = indices[i, 1:]
        neighbor_labels = labels[neighbor_idx]

        values, counts = np.unique(neighbor_labels, return_counts=True)
        majority_label = values[np.argmax(counts)]

        abnormal.append(labels[i] != majority_label)

    return np.mean(abnormal)


# ============================================================
# Load data
# ============================================================

print("\n[1/5] Loading input files...")

coords_df = pd.read_csv(coords_path, sep="\t")
clusters_df = pd.read_csv(clusters_path, sep="\t", index_col=0)

print(f"  Spatial coords: {coords_df.shape[0]} spots x {coords_df.shape[1]} columns")
print(f"  Cluster matrix: {clusters_df.shape[0]} spots x {clusters_df.shape[1]} resolutions")


# ============================================================
# Check input data
# ============================================================

print("\n[2/5] Checking required columns...")

required_cols = {"barcode", "imagerow", "imagecol"}
missing_cols = required_cols - set(coords_df.columns)

if missing_cols:
    raise ValueError(f"Missing required columns in coords file: {missing_cols}")

print("  Required columns found: barcode, imagerow, imagecol")


# ============================================================
# Match coords and clusters by barcode
# ============================================================

print("\n[3/5] Matching spatial coords with cluster labels by barcode...")

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
# Prepare spatial coordinates
# ============================================================

coords = df[["imagerow", "imagecol"]].to_numpy()
cluster_cols = list(clusters_df.columns)

print("\n[4/5] Computing CHAOS and PAS...")
print("  Using tissue coordinates: imagerow, imagecol")
print("  Lower CHAOS = more spatially compact clusters")
print("  Lower PAS   = better local spatial coherence")

results = []

for i, col in enumerate(cluster_cols, start=1):
    resolution = col.replace("cluster_resolution_", "")
    labels = df[col].astype(str).to_numpy()
    n_clusters = len(np.unique(labels))

    print(
        f"  [{i}/{len(cluster_cols)}] "
        f"resolution={resolution} | "
        f"n_clusters={n_clusters}"
    )

    chaos = compute_CHAOS_local(coords, labels)
    pas = compute_PAS_local(coords, labels, k=10)

    results.append({
        "resolution": resolution,
        "cluster_column": col,
        "n_spots": df.shape[0],
        "n_clusters": n_clusters,
        "CHAOS": chaos,
        "PAS": pas
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
