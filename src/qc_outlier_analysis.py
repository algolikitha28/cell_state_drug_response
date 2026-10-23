
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

# --------------------------------------------------
# 1. Project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "data" / "processed" / "GSE140440_counts.h5ad"
QC_DIR = ROOT / "results" / "qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT = QC_DIR / "qc_outlier_review.csv"
SUMMARY = QC_DIR / "qc_outlier_summary.csv"

# --------------------------------------------------
# 2. Load the original expression matrix
# --------------------------------------------------

print("Loading dataset...")
adata = ad.read_h5ad(INPUT)

print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")


# --------------------------------------------------
# 3. Calculate QC metrics when needed
# --------------------------------------------------

# The original H5AD file contains raw counts.
# Calculate QC metrics in memory; do not overwrite the input file.

if "mitochondrial" not in adata.var.columns:
    raise ValueError(
        "Missing mitochondrial gene annotation in adata.var['mitochondrial']."
    )

import scanpy as sc

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mitochondrial"],
    percent_top=None,
    log1p=True,
    inplace=True,
)

required = [
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
]

missing = [col for col in required if col not in adata.obs.columns]

if missing:
    raise ValueError(f"Missing required metadata/QC columns: {missing}")

print("QC metrics calculated successfully.")

# --------------------------------------------------
# 4. Build a table for review
# --------------------------------------------------

qc = adata.obs[required].copy()

qc["log10_total_counts"] = np.log10(
    qc["total_counts"].clip(lower=1)
)

qc["log10_n_genes"] = np.log10(
    qc["n_genes_by_counts"].clip(lower=1)
)


# --------------------------------------------------
# 5. Detect unusually low values within each group
# --------------------------------------------------

# Median absolute deviation (MAD) is a robust measure
# of spread. This method flags statistical outliers;
# it does not prove that a cell is poor quality.

def robust_low_outlier(values, cutoff=3.5):
    values = pd.to_numeric(values, errors="coerce")
    median = values.median()
    mad = (values - median).abs().median()

    # If there is no spread, avoid dividing by zero.
    if pd.isna(mad) or mad == 0:
        return pd.Series(False, index=values.index)

    robust_z = 0.6745 * (values - median) / mad
    return robust_z < -cutoff


qc["low_count_outlier"] = False
qc["low_gene_outlier"] = False

group_columns = ["cell_line", "docetaxel_response"]

for _, group in qc.groupby(group_columns, observed=True):
    qc.loc[group.index, "low_count_outlier"] = (
        robust_low_outlier(group["log10_total_counts"])
    )

    qc.loc[group.index, "low_gene_outlier"] = (
        robust_low_outlier(group["log10_n_genes"])
    )

# --------------------------------------------------
# 6. Add transparent biological QC review flags
# --------------------------------------------------

# These are configurable review thresholds, not
# automatic filtering rules or universal standards.
qc["high_mito_review"] = (
    qc["pct_counts_mitochondrial"] >= 20
)

qc["low_counts_review"] = (
    qc["total_counts"] < 200000
)

qc["low_genes_review"] = (
    qc["n_genes_by_counts"] < 5000
)

flag_columns = [
    "low_count_outlier",
    "low_gene_outlier",
    "high_mito_review",
    "low_counts_review",
    "low_genes_review",
]

qc["number_of_flags"] = qc[flag_columns].sum(axis=1)

qc["review_required"] = qc["number_of_flags"] > 0

# Keep every cell in the output.
qc.to_csv(OUTPUT, index=True, index_label="sample_id")

# --------------------------------------------------
# 7. Summarize flagged cells in each group
# --------------------------------------------------

summary = (
    qc.groupby(group_columns, observed=True)
    .agg(
        original_cells=("total_counts", "size"),
        cells_flagged=("review_required", "sum"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_mito_percent=("pct_counts_mitochondrial", "median"),
    )
    .reset_index()
)

summary.to_csv(SUMMARY, index=False)

print("\n=== QC OUTLIER REVIEW ===")
print(summary.to_string(index=False))

print(f"\nCell-level review: {OUTPUT}")
print(f"Group summary: {SUMMARY}")
print("\nNo cells were removed or filtered.")
