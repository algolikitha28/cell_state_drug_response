
"""
Normalize the QC-filtered single-cell RNA-seq dataset.

Input:
    data/processed/GSE140440_qc_filtered.h5ad

Output:
    data/processed/GSE140440_normalized.h5ad

The original raw-count files are never overwritten.
"""

from pathlib import Path

import anndata as ad
import scanpy as sc


# --------------------------------------------------
# 1. Define input and output paths
# --------------------------------------------------

INPUT_FILE = Path("data/processed/GSE140440_qc_filtered.h5ad")

OUTPUT_FILE = Path("data/processed/GSE140440_normalized.h5ad")


# --------------------------------------------------
# 2. Load the QC-filtered dataset
# --------------------------------------------------

print("Loading QC-filtered dataset...")

adata = ad.read_h5ad(INPUT_FILE)

print(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")


# --------------------------------------------------
# 3. Preserve the unnormalized counts
# --------------------------------------------------

# Store a copy of the count matrix before normalization.
# This allows us to access the original counts later.

adata.layers["counts"] = adata.X.copy()


# --------------------------------------------------
# 4. Normalize counts per cell
# --------------------------------------------------

# Scale each cell so its total expression equals 10,000.
# This helps account for differences in sequencing depth.

print("Normalizing counts to 10,000 per cell...")

sc.pp.normalize_total(adata, target_sum=10000)


# --------------------------------------------------
# 5. Log-transform normalized expression
# --------------------------------------------------

# log1p means log(1 + expression value).
# It reduces the influence of very large expression values.

print("Applying log1p transformation...")

sc.pp.log1p(adata)


# --------------------------------------------------
# 6. Record the processing information
# --------------------------------------------------

adata.uns["normalization"] = {
    "method": "Scanpy normalize_total followed by log1p",
    "target_sum": 10000,
    "raw_counts_layer": "counts",
    "input_file": str(INPUT_FILE),
}


# --------------------------------------------------
# 7. Save the normalized dataset
# --------------------------------------------------

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

adata.write_h5ad(OUTPUT_FILE, compression="gzip")

print("\nNormalization completed successfully.")
print(f"Output file: {OUTPUT_FILE}")
print(f"Output shape: {adata.n_obs} cells x {adata.n_vars} genes")
print("Raw counts preserved in adata.layers['counts']")
