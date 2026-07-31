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
    davies_bouldin_score
)


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
    "risperidone_half_cluster_quality_metrics.tsv"
)

plot_file = os.path.join(
    input_dir,
    "risperidone_half_cluster_quality_metrics.png"
)

chosen_resolution = 0.4
sample_size = 5000
random_state = 123


# =========================
# 3. Wczytanie danych
# =========================

print("[1/6] Wczytuję PCA embeddings...")
X = pd.read_csv(pca_file, sep="\t", index_col=0)

print("[2/6] Wczytuję klastry...")
clusters_df = pd.read_csv(clusters_file, sep="\t", index_col=0)


# =========================
# 4. Dopasowanie ID
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"[3/6] Dopasowano spoty: {len(common_ids)}")


# =========================
# 5. Losowy sample spotów
# =========================

print("[4/6] Przygotowuję sample...")

rng = np.random.default_rng(random_state)

n_total = X.shape[0]
n_sample = min(sample_size, n_total)

sample_idx = rng.choice(
    n_total,
    size=n_sample,
    replace=False
)

X_sample = X.iloc[sample_idx, :]
clusters_sample_df = clusters_df.iloc[sample_idx, :]

print(f"Używam sample: {n_sample} / {n_total} spotów")


# =========================
# 6. Liczenie metryk
# =========================

results = []
cluster_cols = list(clusters_df.columns)

print("[5/6] Liczę metryki jakości klasteryzacji...")

for i, cluster_col in enumerate(cluster_cols, start=1):

    print(f"  [{i}/{len(cluster_cols)}] {cluster_col}")

    labels = clusters_sample_df[cluster_col].values
    n_clusters = len(np.unique(labels))

    resolution = float(cluster_col.replace("cluster_resolution_", ""))

    if n_clusters < 2:
        silhouette = np.nan
        calinski_harabasz = np.nan
        davies_bouldin = np.nan
    else:
        silhouette = silhouette_score(
            X_sample.values,
            labels,
            metric="euclidean"
        )

        calinski_harabasz = calinski_harabasz_score(
            X_sample.values,
            labels
        )

        davies_bouldin = davies_bouldin_score(
            X_sample.values,
            labels
        )

    cluster_sizes = pd.Series(labels).value_counts()

    results.append({
        "resolution": resolution,
        "n_clusters": n_clusters,
        "min_cluster_size": cluster_sizes.min(),
        "median_cluster_size": cluster_sizes.median(),
        "max_cluster_size": cluster_sizes.max(),
        "silhouette_score": silhouette,
        "calinski_harabasz_score": calinski_harabasz,
        "davies_bouldin_score": davies_bouldin
    })


results_df = pd.DataFrame(results).sort_values("resolution")


# =========================
# 7. Zapis wyników
# =========================

print("[6/6] Zapisuję wyniki...")

results_df.to_csv(output_file, sep="\t", index=False)


# =========================
# 8. Wykres
# =========================

fig, axes = plt.subplots(
    nrows=4,
    ncols=1,
    figsize=(8, 10),
    sharex=True
)

# Silhouette
axes[0].plot(
    results_df["resolution"],
    results_df["silhouette_score"],
    marker="o"
)
axes[0].axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5
)
axes[0].set_ylabel("Silhouette\nhigher = better")
axes[0].set_title("Cluster quality metrics (risperidone half)")

# Calinski-Harabasz
axes[1].plot(
    results_df["resolution"],
    results_df["calinski_harabasz_score"],
    marker="o"
)
axes[1].axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5
)
axes[1].set_ylabel("Calinski-Harabasz\nhigher = better")

# Davies-Bouldin
axes[2].plot(
    results_df["resolution"],
    results_df["davies_bouldin_score"],
    marker="o"
)
axes[2].axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5
)
axes[2].set_ylabel("Davies-Bouldin\nlower = better")

# Liczba klastrów
axes[3].plot(
    results_df["resolution"],
    results_df["n_clusters"],
    marker="o"
)
axes[3].axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label=f"Chosen resolution = {chosen_resolution}"
)
axes[3].set_ylabel("Number of clusters")
axes[3].set_xlabel("Resolution")
axes[3].legend()

plt.tight_layout()
plt.savefig(plot_file, dpi=300)

print("Gotowe ✅")
print("TSV:", output_file)
print("PNG:", plot_file)
