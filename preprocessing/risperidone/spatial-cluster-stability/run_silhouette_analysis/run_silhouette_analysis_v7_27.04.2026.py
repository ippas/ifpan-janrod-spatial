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

from sklearn.metrics import silhouette_score


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
    "risperidone_half_silhouette_multi_distance_metrics.tsv"
)

plot_file = os.path.join(
    input_dir,
    "risperidone_half_silhouette_multi_distance_metrics.png"
)

chosen_resolution = 0.4

sample_size = 5000
random_state = 123

distance_metrics = [
    "euclidean",
    "cosine",
    "correlation",
    "manhattan",
    "chebyshev",
    "braycurtis",
    "canberra"
]


# =========================
# 3. Wczytanie danych
# =========================

print("[1/6] Wczytuję PCA embeddings...")

X = pd.read_csv(
    pca_file,
    sep="\t",
    index_col=0
)

print("[2/6] Wczytuję klastry...")

clusters_df = pd.read_csv(
    clusters_file,
    sep="\t",
    index_col=0
)


# =========================
# 4. Dopasowanie ID
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"[3/6] Dopasowano spoty: {X.shape[0]}")
print(f"Liczba PC: {X.shape[1]}")


# =========================
# 5. Liczenie silhouette dla wielu metryk
# =========================

print("[4/6] Liczę silhouette dla różnych metryk odległości...")

results = []

cluster_cols = list(clusters_df.columns)

for i, cluster_col in enumerate(cluster_cols, start=1):

    resolution = float(
        cluster_col.replace("cluster_resolution_", "")
    )

    labels = clusters_df[cluster_col].values
    n_clusters = clusters_df[cluster_col].nunique()

    print(f"  [{i}/{len(cluster_cols)}] resolution = {resolution}")

    for metric in distance_metrics:

        print(f"      metric = {metric}")

        if n_clusters < 2:
            score = np.nan
        else:
            try:
                score = silhouette_score(
                    X.values,
                    labels,
                    metric=metric,
                    sample_size=min(sample_size, X.shape[0]),
                    random_state=random_state
                )

            except Exception as e:
                print(f"      [WARNING] Nie udało się policzyć {metric}: {e}")
                score = np.nan

        results.append({
            "resolution": resolution,
            "distance_metric": metric,
            "silhouette_score": score,
            "n_clusters": n_clusters,
            "sample_size": min(sample_size, X.shape[0])
        })


results_df = pd.DataFrame(results).sort_values(
    ["distance_metric", "resolution"]
)


# =========================
# 6. Zapis TSV
# =========================

print("[5/6] Zapisuję wyniki...")

results_df.to_csv(
    output_file,
    sep="\t",
    index=False
)


# =========================
# 7. Wykres zbiorczy
# =========================

print("[6/6] Tworzę wykres...")

plt.figure(figsize=(9, 5))

for metric in distance_metrics:

    tmp = results_df[
        results_df["distance_metric"] == metric
    ].sort_values("resolution")

    plt.plot(
        tmp["resolution"],
        tmp["silhouette_score"],
        marker="o",
        linewidth=1.5,
        label=metric
    )

plt.axvline(
    x=chosen_resolution,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label=f"Chosen resolution = {chosen_resolution}"
)

plt.xlabel("Resolution")
plt.ylabel("Silhouette score")
plt.title("Silhouette analysis using multiple distance metrics")

plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc="upper left"
)

plt.tight_layout()
plt.savefig(plot_file, dpi=300)

print("Gotowe ✅")
print("TSV:", output_file)
print("PNG:", plot_file)
