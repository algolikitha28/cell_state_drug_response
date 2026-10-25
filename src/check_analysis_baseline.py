from pathlib import Path
import scanpy as sc
import pandas as pd

files = [
    Path("data/processed/GSE140440_counts.h5ad"),
    Path("data/processed/GSE140440_qc_filtered.h5ad"),
    Path("data/processed/GSE140440_gene_filtered.h5ad"),
]

for path in files:
    print(f"\n{'=' * 60}")
    print(f"FILE: {path}")

    if not path.exists():
        print("STATUS: FILE NOT FOUND")
        continue

    adata = sc.read_h5ad(path)
    print("Shape (samples/cells, genes):", adata.shape)
    print("Raw counts layer available:", "counts" in adata.layers)

    if "cell_line" in adata.obs:
        print("\nCell-line counts:")
        print(adata.obs["cell_line"].value_counts().to_string())

    if "docetaxel_response" in adata.obs:
        print("\nResponse counts:")
        print(adata.obs["docetaxel_response"].value_counts().to_string())

    if {"cell_line", "docetaxel_response"}.issubset(adata.obs.columns):
        print("\nCell-line and response groups:")
        print(
            adata.obs.groupby(
                ["cell_line", "docetaxel_response"],
                observed=True
            ).size().to_string()
        )
