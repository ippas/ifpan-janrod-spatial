#!/usr/bin/env python3

import os
import resource

# =========================
# 0. Limit RAM i wątki
# =========================

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

limit_gb = 50
limit_bytes = limit_gb * 1024**3
resource.setrlimit(resource.RLIMIT_AS, (limit_bytes, limit_bytes))

print(f"[DEBUG] Ustawiono limit RAM: {limit_gb} GB")


# =========================
# 1. Importy
# =========================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    adjusted_rand_score,
    normalized_mutual_info_score
)

from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

import scanpy as sc
import anndata as ad


# =========================
# 2. Parametry
# =========================

input_dir = "data/risperidone/silhouette_analysis"

pca_file = os.path.join(
    input_dir,
    "risperidone_half_pca_embeddings_1_30.tsv"
)

clusters_file = os.path.join(
    input_dir,
    "risperidone_half_clusters.tsv"
)

output_file = os.path.join(
    input_dir,
    "resolution_compare_0.4_vs_0.85_FAST.tsv"
)

plot_file = os.path.join(
    input_dir,
    "resolution_compare_0.4_vs_0.85_FAST.png"
)

dendrogram_dir = os.path.join(
    input_dir,
    "dendrograms_0.4_vs_0.85_FAST"
)

os.makedirs(dendrogram_dir, exist_ok=True)

target_resolutions = [0.4, 0.85]

# bardzo luźne/szybkie parametry
sample_size = 2000
n_stability_runs = 3
subsample_fraction = 0.5
n_neighbors = 10

min_cluster_size = 20
random_state = 123


# =========================
# 3. Funkcje pomocnicze
# =========================

def get_resolution_from_col(colname):
    return float(colname.replace("cluster_resolution_", ""))


def run_leiden(X_array, resolution, seed):
    """
    Szybka ponowna klasteryzacja Leiden na PCA embeddings.
    """

    adata = ad.AnnData(X_array)

    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        use_rep="X",
        random_state=seed
    )

    sc.tl.leiden(
        adata,
        resolution=resolution,
        random_state=seed,
        key_added="leiden",
        flavor="igraph",
        n_iterations=2,
        directed=False
    )

    return adata.obs["leiden"].astype(str).values


def save_dendrogram(X_df, labels, resolution, output_dir):
    """
    Dendrogram centroidów klastrów.
    """

    tmp = X_df.copy()
    tmp["cluster"] = labels

    centroids = tmp.groupby("cluster").mean()

    if centroids.shape[0] < 2:
        return np.nan

    Z = linkage(centroids.values, method="ward")

    plt.figure(figsize=(8, 5))
    dendrogram(
        Z,
        labels=centroids.index.astype(str).tolist(),
        leaf_rotation=90
    )

    plt.title(f"Cluster dendrogram, resolution = {resolution}")
    plt.xlabel("Cluster")
    plt.ylabel("Ward distance")
    plt.tight_layout()

    out_file = os.path.join(
        output_dir,
        f"dendrogram_resolution_{resolution}.png"
    )

    plt.savefig(out_file, dpi=300)
    plt.close()

    centroid_distances = pdist(centroids.values, metric="euclidean")
    return float(np.mean(centroid_distances))


# =========================
# 4. Wczytanie danych
# =========================

print("[1/7] Wczytuję PCA embeddings...")
X = pd.read_csv(pca_file, sep="\t", index_col=0)

print("[2/7] Wczytuję klastry...")
clusters_df = pd.read_csv(clusters_file, sep="\t", index_col=0)


# =========================
# 5. Dopasowanie ID
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"[3/7] Dopasowano spoty: {X.shape[0]}")
print(f"Liczba PC: {X.shape[1]}")


# =========================
# 6. Wybór resolution 0.4 i 0.85
# =========================

print("[4/7] Wybieram resolution 0.4 i 0.85...")

selected_cols = []

for col in clusters_df.columns:
    res = get_resolution_from_col(col)

    if any(abs(res - target) < 1e-6 for target in target_resolutions):
        selected_cols.append(col)

if len(selected_cols) == 0:
    raise ValueError("Nie znaleziono kolumn dla resolution 0.4 lub 0.85")

clusters_df = clusters_df[selected_cols]

print("Wybrane kolumny:")
for col in selected_cols:
    print(" -", col)


# =========================
# 7. Analiza
# =========================

print("[5/7] Liczę metryki i stability...")

rng = np.random.default_rng(random_state)

results = []

for col in clusters_df.columns:

    resolution = get_resolution_from_col(col)

    print(f"\nResolution = {resolution}")

    labels = clusters_df[col].astype(str)

    cluster_sizes = labels.value_counts()

    n_clusters = int(cluster_sizes.shape[0])
    min_size = int(cluster_sizes.min())
    median_size = float(cluster_sizes.median())
    max_size = int(cluster_sizes.max())

    n_small_clusters = int((cluster_sizes < min_cluster_size).sum())
    n_spots_in_small_clusters = int(
        cluster_sizes[cluster_sizes < min_cluster_size].sum()
    )

    percent_small_clusters = 100 * n_small_clusters / n_clusters
    percent_spots_in_small_clusters = (
        100 * n_spots_in_small_clusters / labels.shape[0]
    )

    # -------------------------
    # sample do metryk
    # -------------------------

    n_sample = min(sample_size, X.shape[0])

    sample_idx = rng.choice(
        X.shape[0],
        size=n_sample,
        replace=False
    )

    X_sample = X.iloc[sample_idx].values
    labels_sample = labels.iloc[sample_idx].values

    if len(np.unique(labels_sample)) < 2:
        silhouette = np.nan
        calinski_harabasz = np.nan
        davies_bouldin = np.nan
    else:
        silhouette = silhouette_score(
            X_sample,
            labels_sample,
            metric="euclidean"
        )

        calinski_harabasz = calinski_harabasz_score(
            X_sample,
            labels_sample
        )

        davies_bouldin = davies_bouldin_score(
            X_sample,
            labels_sample
        )

    # -------------------------
    # dendrogram centroidów
    # -------------------------

    mean_centroid_distance = save_dendrogram(
        X_df=X,
        labels=labels.values,
        resolution=resolution,
        output_dir=dendrogram_dir
    )

    # -------------------------
    # szybka stability analysis
    # -------------------------

    ari_scores = []
    nmi_scores = []

    for run in range(n_stability_runs):

        seed = random_state + run
        rng_run = np.random.default_rng(seed)

        n_subsample = int(X.shape[0] * subsample_fraction)

        idx_sub = rng_run.choice(
            X.shape[0],
            size=n_subsample,
            replace=False
        )

        X_sub = X.iloc[idx_sub].values

        labels_sub_new = run_leiden(
            X_array=X_sub,
            resolution=resolution,
            seed=seed
        )

        labels_sub_original = labels.iloc[idx_sub].values

        ari = adjusted_rand_score(
            labels_sub_original,
            labels_sub_new
        )

        nmi = normalized_mutual_info_score(
            labels_sub_original,
            labels_sub_new
        )

        ari_scores.append(ari)
        nmi_scores.append(nmi)

    results.append({
        "resolution": resolution,

        "n_clusters": n_clusters,
        "min_cluster_size": min_size,
        "median_cluster_size": median_size,
        "max_cluster_size": max_size,

        "n_small_clusters_lt20": n_small_clusters,
        "percent_small_clusters_lt20": percent_small_clusters,
        "n_spots_in_small_clusters_lt20": n_spots_in_small_clusters,
        "percent_spots_in_small_clusters_lt20": percent_spots_in_small_clusters,

        "silhouette_score": silhouette,
        "calinski_harabasz_score": calinski_harabasz,
        "davies_bouldin_score": davies_bouldin,

        "mean_centroid_distance": mean_centroid_distance,

        "mean_stability_ARI": float(np.mean(ari_scores)),
        "sd_stability_ARI": float(np.std(ari_scores)),
        "mean_stability_NMI": float(np.mean(nmi_scores)),
        "sd_stability_NMI": float(np.std(nmi_scores)),

        "sample_size_metrics": n_sample,
        "n_stability_runs": n_stability_runs,
        "subsample_fraction": subsample_fraction,
        "n_neighbors": n_neighbors
    })


results_df = pd.DataFrame(results).sort_values("resolution")

results_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


# =========================
# 8. Wykres
# =========================

print("[6/7] Tworzę wykres...")

x_labels = results_df["resolution"].astype(str)

fig, axes = plt.subplots(
    nrows=7,
    ncols=1,
    figsize=(7, 14),
    sharex=True
)

axes[0].bar(x_labels, results_df["silhouette_score"])
axes[0].set_ylabel("Silhouette\nhigher = better")
axes[0].set_title("Resolution comparison: 0.4 vs 0.85")

axes[1].bar(x_labels, results_df["calinski_harabasz_score"])
axes[1].set_ylabel("Calinski-\nHarabasz")

axes[2].bar(x_labels, results_df["davies_bouldin_score"])
axes[2].set_ylabel("Davies-\nBouldin")

axes[3].bar(x_labels, results_df["n_clusters"])
axes[3].set_ylabel("Number of\nclusters")

axes[4].bar(x_labels, results_df["min_cluster_size"])
axes[4].axhline(
    min_cluster_size,
    color="black",
    linestyle=":",
    linewidth=1
)
axes[4].set_ylabel("Minimum\ncluster size")

axes[5].bar(x_labels, results_df["percent_spots_in_small_clusters_lt20"])
axes[5].set_ylabel("% spots in\nclusters <20")

axes[6].bar(x_labels, results_df["mean_stability_ARI"])
axes[6].errorbar(
    x=np.arange(len(results_df)),
    y=results_df["mean_stability_ARI"],
    yerr=results_df["sd_stability_ARI"],
    fmt="none",
    capsize=4
)
axes[6].set_ylabel("Stability\nARI")
axes[6].set_xlabel("Resolution")

plt.tight_layout()
plt.savefig(plot_file, dpi=300)


# =========================
# 9. Koniec
# =========================

print("[7/7] Gotowe ✅")
print("TSV:", output_file)
print("PNG:", plot_file)
print("Dendrogramy:", dendrogram_dir)
