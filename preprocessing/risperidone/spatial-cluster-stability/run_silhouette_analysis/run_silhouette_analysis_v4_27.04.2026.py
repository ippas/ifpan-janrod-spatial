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
    "risperidone_half_silhouette_scores_full.tsv"
)

plot_file = os.path.join(
    input_dir,
    "risperidone_half_silhouette_scores_full_y0_1.png"
)

sample_size = None  # None = wszystkie spoty
random_state = 123


# =========================
# 3. Wczytanie danych
# =========================

print("[1/5] Wczytuję PCA embeddings...")
X = pd.read_csv(pca_file, sep="\t", index_col=0)

print("[2/5] Wczytuję klastry...")
clusters_df = pd.read_csv(clusters_file, sep="\t", index_col=0)


# =========================
# 4. Dopasowanie spotów
# =========================

common_ids = X.index.intersection(clusters_df.index)

X = X.loc[common_ids]
clusters_df = clusters_df.loc[common_ids]

print(f"[3/5] Dopasowano spoty: {len(common_ids)}")
print(f"      PCA shape: {X.shape}")
print(f"      clusters shape: {clusters_df.shape}")


# =========================
# 5. Silhouette dla wszystkich resolution
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
            sample_size=sample_size,
            random_state=random_state
        )

    resolution = float(cluster_col.replace("cluster_resolution_", ""))

    results.append({
        "cluster_col": cluster_col,
        "resolution": resolution,
        "n_clusters": n_clusters,
        "n_spots": X.shape[0],
        "sample_size": X.shape[0] if sample_size is None else min(sample_size, X.shape[0]),
        "silhouette_score": sil_score
    })


results_df = pd.DataFrame(results).sort_values("resolution")


# =========================
# 6. Zapis wyników
# =========================

print("[5/5] Zapisuję wyniki...")

results_df.to_csv(output_file, sep="\t", index=False)


# =========================
# 7. Wykres z osią Y 0-1
# =========================

plt.figure(figsize=(7, 4))

plt.plot(
    results_df["resolution"],
    results_df["silhouette_score"],
    marker="o"
)

plt.xlabel("Resolution")
plt.ylabel("Silhouette score")
plt.title("Silhouette analysis - risperidone half")

plt.ylim(0, 1)

plt.tight_layout()
plt.savefig(plot_file, dpi=300)

print("Gotowe.")
print("TSV:", output_file)
print("PNG:", plot_file)
