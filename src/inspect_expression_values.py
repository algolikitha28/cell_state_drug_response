from pathlib import Path
import scanpy as sc
import numpy as np

files = [
    "data/processed/GSE140440_counts.h5ad",
    "data/processed/GSE140440_qc_filtered.h5ad",
    "data/processed/GSE140440_gene_filtered.h5ad",
]

for filename in files:
    print("\n" + "=" * 60)
    print("FILE:", filename)

    adata = sc.read_h5ad(Path(filename))
    X = adata.X

    if hasattr(X, "data"):
        values = X.data
    else:
        values = np.asarray(X).ravel()

    print("Shape:", adata.shape)
    print("Minimum stored value:", values.min() if values.size else "empty")
    print("Maximum stored value:", values.max() if values.size else "empty")
    print("All values non-negative:", bool(np.all(values >= 0)))
    print("All values approximately integers:",
          bool(np.allclose(values, np.round(values))))
    print("Counts layer:", "counts" in adata.layers)

    if "counts" in adata.layers:
        layer = adata.layers["counts"]
        layer_values = layer.data if hasattr(layer, "data") else np.asarray(layer).ravel()
        print("Counts layer maximum:",
              layer_values.max() if layer_values.size else "empty")
        print("Counts layer approximately integers:",
              bool(np.allclose(layer_values, np.round(layer_values))))
