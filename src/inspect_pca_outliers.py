
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

# Load the normalized dataset
data = ad.read_h5ad(
    r"data/processed/GSE140440_normalized.h5ad"
)

output_dir = Path(r"results/eda/pca_by_cell_line")
output_dir.mkdir(parents=True, exist_ok=True)

# Analyze DU145 and PC3 separately
for line in ["DU145", "PC3"]:

    print(f"\n{'=' * 60}")
    print(f"{line}: MOST DISTANT CELLS IN PCA")
    print("=" * 60)

    # Select cells from one cell line
    sub = data[data.obs["cell_line"] == line].copy()

    # Select the 2,000 most highly variable genes
    sc.pp.highly_variable_genes(
        sub,
        n_top_genes=min(2000, sub.n_vars),
        flavor="seurat"
    )

    sub = sub[:, sub.var["highly_variable"]].copy()

    # Scale gene expression and calculate PCA
    sc.pp.scale(sub, max_value=10)

    sc.tl.pca(
        sub,
        n_comps=min(30, sub.n_obs - 1, sub.n_vars - 1),
        random_state=42
    )

    # Extract the first two principal components
    scores = pd.DataFrame(
        sub.obsm["X_pca"][:, :2],
        columns=["PC1", "PC2"],
        index=sub.obs_names
    )

    scores["geo_accession"] = (
        sub.obs["geo_accession"].astype(str).values
    )
    scores["response"] = (
        sub.obs["docetaxel_response"].astype(str).values
    )

    # Calculate distance from the PCA plot's origin
    scores["distance"] = (
        scores["PC1"] ** 2 + scores["PC2"] ** 2
    ) ** 0.5

    # Display the 10 most distant cells
    top10 = scores.nlargest(10, "distance")

    print(
        top10[
            ["geo_accession", "response", "PC1", "PC2", "distance"]
        ].round(2).to_string(index=False)
    )
