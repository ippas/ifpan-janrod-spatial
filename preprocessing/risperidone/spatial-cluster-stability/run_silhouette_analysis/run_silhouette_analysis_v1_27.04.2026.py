#!/usr/bin/env python3

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score


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
    "risperidone_half_silhouette_scores.tsv"
)

plot_file = os.path.join(
    input_dir,
    "risperidone_half_silhouette_scores.png"
)

sample_size = 5000
random_state = 123


# =========================
# 2. Wczytanie danych
# =========================

print("[1/5] Wczytuję PCA embeddings...")
X = pd.read_csv(pca_file, sep="\t", index_col=0)

print("[2/5] Wczytuję klastry...")
clusters_df = pd.read_csv(clusters_file, sep="\t", index_col=0)


# =========================
# 3. Dopasowanie spotów
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"[3/5] Dopasowano spoty: {len(common_ids)}")


# =========================
# 4. Silhouette
# =========================

results = []

cluster_cols = list(clusters_df.columns)

print("[4/5] Liczę silhouette...")

for i, cluster_col in enumerate(cluster_cols, start=1):

    print(f"  [{i}/{len(cluster_cols)}] {cluster_col}")

    clusters = clusters_df[cluster_col].values
    n_clusters = clusters_df[cluster_col].nunique()

    if n_clusters < 2:
        sil_score = None
    else:
        sil_score = silhouette_score(
            X.values,
            clusters,
            metric="euclidean",
            sample_size=min(sample_size, X.shape[0]),
            random_state=random_state
        )

    resolution = float(cluster_col.replace("cluster_resolution_", ""))

    results.append({
        "cluster_col": cluster_col,
        "resolution": resolution,
        "n_clusters": n_clusters,
        "n_spots": X.shape[0],
        "sample_size": min(sample_size, X.shape[0]),
        "silhouette_score": sil_score
    })


results_df = pd.DataFrame(results).sort_values("resolution")


# =========================
# 5. Zapis
# =========================

print("[5/5] Zapisuję wyniki...")

results_df.to_csv(output_file, sep="\t", index=False)

plt.figure(figsize=(7, 4))
plt.plot(
    results_df["resolution"],
    results_df["silhouette_score"],
    marker="o"
)
plt.xlabel("Resolution")
plt.ylabel("Silhouette score")
plt.title("Silhouette analysis (risperidone half)")
plt.tight_layout()
plt.savefig(plot_file, dpi=300)

print("Gotowe ✅")
print("TSV:", output_file)
print("PNG:", plot_file)
