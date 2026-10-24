
"""
Exploratory data analysis for StateDrugResponse.

Input:
    data/processed/GSE140440_normalized.h5ad

Outputs:
    results/eda/ - QC comparisons and PCA plots

This script does not modify the input dataset.
"""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


# --------------------------------------------------
# 1. Define paths
# --------------------------------------------------

INPUT_FILE = Path("data/processed/GSE140440_normalized.h5ad")
OUTPUT_DIR = Path("results/eda")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# 2. Load data and verify metadata
# --------------------------------------------------

print("Loading normalized data...")
adata = ad.read_h5ad(INPUT_FILE)

required_columns = [
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
]

missing_columns = [
    col for col in required_columns if col not in adata.obs.columns
]

if missing_columns:
    raise ValueError(f"Missing metadata columns: {missing_columns}")

print(f"Loaded {adata.n_obs} cells and {adata.n_vars} genes.")


# --------------------------------------------------
# 3. Create a combined group label
# --------------------------------------------------

# Each cell belongs to one cell-line/response group.
adata.obs["group"] = (
    adata.obs["cell_line"].astype(str)
    + " - "
    + adata.obs["docetaxel_response"].astype(str)
)

group_order = [
    "DU145 - Sensitive",
    "DU145 - Resistant",
    "PC3 - Sensitive",
    "PC3 - Resistant",
]


# --------------------------------------------------
# 4. Save descriptive statistics
# --------------------------------------------------

metrics = [
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
]

summary = (
    adata.obs.groupby(
        ["cell_line", "docetaxel_response"], observed=True
    )[metrics]
    .agg(["count", "mean", "median", "std", "min", "max"])
)

summary.to_csv(OUTPUT_DIR / "qc_descriptive_summary.csv")

group_counts = (
    adata.obs.groupby(
        ["cell_line", "docetaxel_response"], observed=True
    )
    .size()
    .rename("cell_count")
    .reset_index()
)

group_counts.to_csv(
    OUTPUT_DIR / "group_cell_counts.csv", index=False
)

print("\nGroup cell counts:")
print(group_counts.to_string(index=False))


# --------------------------------------------------
# 5. Plot QC metrics by group
# --------------------------------------------------

plot_metrics = [
    ("total_counts", "Total counts per cell",
     "total_counts_by_group.png", True),

    ("n_genes_by_counts", "Detected genes per cell",
     "genes_by_group.png", False),

    ("pct_counts_mitochondrial", "Mitochondrial counts (%)",
     "mitochondrial_percentage_by_group.png", False),
]

for column, ylabel, filename, use_log in plot_metrics:
    fig, ax = plt.subplots(figsize=(10, 6))

    values = []
    labels = []

    for group in group_order:
        group_values = adata.obs.loc[
            adata.obs["group"] == group, column
        ].dropna()

        values.append(group_values.to_numpy())
        labels.append(group)

    ax.boxplot(values, tick_labels=labels, showfliers=True)

    if use_log:
        ax.set_yscale("log")

    ax.set_ylabel(ylabel)
    ax.set_title(f"{ylabel} across cell-line/response groups")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / filename, dpi=300)
    plt.close(fig)

    print(f"Saved {filename}")


# --------------------------------------------------
# 6. Select highly variable genes for PCA
# --------------------------------------------------

# Highly variable genes capture informative expression variation.
# This selection is for exploratory visualization only.

print("\nSelecting highly variable genes...")

# Work on a separate copy so the saved input stays unchanged.
pca_data = adata.copy()

sc.pp.highly_variable_genes(
    pca_data,
    n_top_genes=min(2000, pca_data.n_vars),
    flavor="seurat",
)

hvg_mask = pca_data.var["highly_variable"].to_numpy()

if hvg_mask.sum() < 2:
    raise ValueError("Too few highly variable genes for PCA.")

pca_data = pca_data[:, hvg_mask].copy()

print(f"Selected {pca_data.n_vars} highly variable genes.")


# --------------------------------------------------
# 7. Scale data and run PCA
# --------------------------------------------------

print("Scaling expression values and running PCA...")

# Scaling limits the influence of differences in gene variance.
# max_value clips extreme scaled values.

sc.pp.scale(pca_data, max_value=10)

sc.tl.pca(
    pca_data,
    n_comps=min(30, pca_data.n_obs - 1, pca_data.n_vars - 1),
    random_state=42,
)

# Save variance explained by each principal component.
variance = pca_data.uns["pca"]["variance_ratio"]

pd.DataFrame({
    "principal_component": [
        f"PC{i + 1}" for i in range(len(variance))
    ],
    "variance_explained": variance,
}).to_csv(OUTPUT_DIR / "pca_variance.csv", index=False)


# --------------------------------------------------
# 8. Plot the first two principal components
# --------------------------------------------------

coordinates = pca_data.obsm["X_pca"]

fig, ax = plt.subplots(figsize=(9, 7))

for group in group_order:
    mask = adata.obs["group"].to_numpy() == group

    ax.scatter(
        coordinates[mask, 0],
        coordinates[mask, 1],
        label=group,
        alpha=0.75,
        s=35,
    )

ax.set_xlabel(
    f"PC1 ({variance[0] * 100:.1f}% variance explained)"
)
ax.set_ylabel(
    f"PC2 ({variance[1] * 100:.1f}% variance explained)"
)
ax.set_title("PCA of normalized single-cell expression")
ax.legend(frameon=False)
ax.grid(alpha=0.2)

fig.tight_layout()
fig.savefig(OUTPUT_DIR / "pca_by_cell_line_and_response.png", dpi=300)
plt.close(fig)


# --------------------------------------------------
# 9. Finish
# --------------------------------------------------

print("\nEDA completed successfully.")
print(f"Results saved in: {OUTPUT_DIR.resolve()}")
print("The normalized input file was not modified.")
