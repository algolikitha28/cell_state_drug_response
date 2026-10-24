
"""Run separate exploratory PCA analyses for DU145 and PC3."""

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


INPUT = Path("data/processed/GSE140440_normalized.h5ad")
OUTPUT = Path("results/eda/pca_by_cell_line")
OUTPUT.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(INPUT)

for cell_line in ["DU145", "PC3"]:
    print(f"\nAnalyzing {cell_line}...")

    # Subset the data without changing the saved input file.
    subset = adata[adata.obs["cell_line"] == cell_line].copy()

    # Select highly variable genes within this cell line.
    sc.pp.highly_variable_genes(
        subset,
        n_top_genes=min(2000, subset.n_vars),
        flavor="seurat",
    )

    subset = subset[:, subset.var["highly_variable"]].copy()

    # Scale the selected expression matrix for PCA.
    sc.pp.scale(subset, max_value=10)

    # Calculate principal components.
    n_components = min(30, subset.n_obs - 1, subset.n_vars - 1)

    sc.tl.pca(
        subset,
        n_comps=n_components,
        random_state=42,
    )

    # Save the variance explained.
    variance = subset.uns["pca"]["variance_ratio"]

    pd.DataFrame({
        "principal_component": [
            f"PC{i + 1}" for i in range(len(variance))
        ],
        "variance_explained": variance,
    }).to_csv(
        OUTPUT / f"pca_variance_{cell_line}.csv",
        index=False,
    )

    # Plot sensitive and resistant cells.
    coordinates = subset.obsm["X_pca"]
    response = subset.obs["docetaxel_response"].astype(str)

    fig, ax = plt.subplots(figsize=(8, 6))

    for group in ["Sensitive", "Resistant"]:
        mask = response.to_numpy() == group

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
    ax.set_title(f"{cell_line}: PCA by docetaxel response")
    ax.legend(frameon=False)
    ax.grid(alpha=0.2)

    fig.tight_layout()
    fig.savefig(OUTPUT / f"pca_{cell_line}.png", dpi=300)
    plt.close(fig)

    print(f"Cells: {subset.n_obs}")
    print(f"Highly variable genes: {subset.n_vars}")
    print(f"PC1 variance: {variance[0] * 100:.2f}%")
    print(f"PC2 variance: {variance[1] * 100:.2f}%")

print("\nSeparate PCA analyses completed.")
print(f"Results saved in: {OUTPUT.resolve()}")
