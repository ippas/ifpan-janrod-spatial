#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.neighbors import NearestNeighbors
from sklearn.metrics import adjusted_rand_score

import igraph as ig
import leidenalg


# =========================
# 1. Parametry
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
    "risperidone_half_cluster_stability_ARI.tsv"
)

plot_file = os.path.join(
    input_dir,
    "risperidone_half_cluster_stability_ARI.png"
)

subsample_fraction = 0.8
n_iterations = 30
k_neighbors = 20
random_state = 123

chosen_resolution = 0.4


# =========================
# 2. Wczytanie danych
# =========================

print("[1/5] Wczytuję PCA embeddings...")

X = pd.read_csv(
    pca_file,
    sep="\t",
    index_col=0
)

print("[2/5] Wczytuję pełne klastry...")

clusters_df = pd.read_csv(
    clusters_file,
    sep="\t",
    index_col=0
)


# =========================
# 3. Dopasowanie ID
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"Dopasowano spoty: {X.shape[0]}")
print(f"Liczba PC: {X.shape[1]}")


# =========================
# 4. Funkcje
# =========================

def build_knn_graph(X_sub, k=20):
    """
    Buduje nieskierowany graf kNN na podstawie PCA embeddings.
    """

    nn = NearestNeighbors(
        n_neighbors=k + 1,
        metric="euclidean"
    )

    nn.fit(X_sub)

    neighbors = nn.kneighbors(
        X_sub,
        return_distance=False
    )

    edges = []

    for i in range(neighbors.shape[0]):
        for j in neighbors[i, 1:]:
            edges.append((i, j))

    graph = ig.Graph(
        n=X_sub.shape[0],
        edges=edges,
        directed=False
    )

    graph.simplify()

    return graph


def run_leiden(graph, resolution):
    """
    Klasteryzacja Leiden z parametrem resolution.
    """

    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        seed=random_state
    )

    return np.array(partition.membership)


# =========================
# 5. Stability ARI
# =========================

print("[3/5] Liczę stability ARI...")

rng = np.random.default_rng(random_state)

cluster_cols = list(clusters_df.columns)
results = []

n_total = X.shape[0]
n_subsample = int(n_total * subsample_fraction)

X_values = X.values

for cluster_col in cluster_cols:

    resolution = float(
        cluster_col.replace("cluster_resolution_", "")
    )

    full_labels = clusters_df[cluster_col].values

    print(f"Resolution = {resolution}")

    for iteration in range(1, n_iterations + 1):

        sampled_idx = rng.choice(
            n_total,
            size=n_subsample,
            replace=False
        )

        X_sub = X_values[sampled_idx, :]
        full_labels_sub = full_labels[sampled_idx]

        graph = build_knn_graph(
            X_sub,
            k=k_neighbors
        )

        subsample_labels = run_leiden(
            graph,
            resolution=resolution
        )

        ari = adjusted_rand_score(
            full_labels_sub,
            subsample_labels
        )

        results.append({
            "resolution": resolution,
            "iteration": iteration,
            "ARI": ari,
            "n_full_clusters": len(np.unique(full_labels_sub)),
            "n_subsample_clusters": len(np.unique(subsample_labels)),
            "subsample_fraction": subsample_fraction,
            "n_spots": n_subsample,
            "k_neighbors": k_neighbors
        })


results_df = pd.DataFrame(results)

summary_df = (
    results_df
    .groupby("resolution")
    .agg(
        mean_ARI=("ARI", "mean"),
        sd_ARI=("ARI", "std"),
        median_ARI=("ARI", "median"),
        mean_n_subsample_clusters=("n_subsample_clusters", "mean"),
        mean_n_full_clusters=("n_full_clusters", "mean")
    )
    .reset_index()
    .sort_values("resolution")
)


# =========================
# 6. Zapis wyników
# =========================

print("[4/5] Zapisuję wyniki...")

summary_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


# =========================
# 7. Wykres
# =========================

print("[5/5] Tworzę wykres...")

plt.figure(figsize=(8, 5))

plt.errorbar(
    summary_df["resolution"],
    summary_df["mean_ARI"],
    yerr=summary_df["sd_ARI"],
    marker="o",
    linewidth=1.5,
    capsize=3,
    label="Mean ARI ± SD"
)

plt.axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label=f"Chosen resolution = {chosen_resolution}"
)

plt.xlabel("Resolution")
plt.ylabel("ARI stability")
plt.title("Cluster stability across subsampling")
plt.ylim(0, 1)

plt.legend()
plt.tight_layout()
plt.savefig(plot_file, dpi=300)

print("Gotowe ✅")
print("TSV:", output_file)
print("PNG:", plot_file)
