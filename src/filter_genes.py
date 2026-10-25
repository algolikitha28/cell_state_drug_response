
"""
Auditable gene filtering for StateDrugResponse.

Input:
    data/processed/GSE140440_normalized.h5ad

Outputs:
    data/processed/GSE140440_gene_filtered.h5ad
    results/gene_filtering/gene_filter_audit.csv
    results/gene_filtering/gene_filter_summary.csv

Important:
    The lowest-15% ranking statistic is a documented approximation.
    It is not confirmed to reproduce the Fluidigm Singular algorithm.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


# --------------------------------------------------
# 1. Paths and filtering settings
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "data" / "processed" / "GSE140440_normalized.h5ad"

OUTPUT = ROOT / "data" / "processed" / "GSE140440_gene_filtered.h5ad"

RESULTS_DIR = ROOT / "results" / "gene_filtering"

AUDIT_OUTPUT = RESULTS_DIR / "gene_filter_audit.csv"
SUMMARY_OUTPUT = RESULTS_DIR / "gene_filter_summary.csv"

MIN_CELLS = 10
LOWEST_FRACTION = 0.15

GROUPS = [
    ("DU145", "Sensitive"),
    ("DU145", "Resistant"),
    ("PC3", "Sensitive"),
    ("PC3", "Resistant"),
]


# --------------------------------------------------
# 2. Load and validate the input
# --------------------------------------------------

if not INPUT.exists():
    raise FileNotFoundError(f"Input file not found: {INPUT}")

adata = ad.read_h5ad(INPUT)

required_columns = {"cell_line", "docetaxel_response"}
missing = required_columns - set(adata.obs.columns)

if missing:
    raise ValueError(f"Missing metadata columns: {sorted(missing)}")

if not adata.var_names.is_unique:
    raise ValueError("Gene identifiers are not unique.")

if not adata.obs_names.is_unique:
    raise ValueError("Cell identifiers are not unique.")

if "counts" not in adata.layers:
    raise ValueError("Raw counts layer 'counts' is missing.")

counts = adata.layers["counts"]

if counts.shape != adata.shape:
    raise ValueError("Raw counts layer shape does not match the dataset.")

if not np.all(np.isfinite(counts.data if hasattr(counts, "data") else counts)):
    raise ValueError("Raw counts contain non-finite values.")

if hasattr(counts, "data"):
    count_values = counts.data
else:
    count_values = np.asarray(counts).ravel()

if np.any(count_values < 0):
    raise ValueError("Raw counts contain negative values.")

if not np.allclose(count_values, np.round(count_values)):
    raise ValueError("Counts layer does not contain integer-like raw counts.")

cell_line = adata.obs["cell_line"].astype(str).to_numpy()
response = adata.obs["docetaxel_response"].astype(str).to_numpy()

group_masks = {}

for line, status in GROUPS:
    key = f"{line}_{status}"
    mask = (cell_line == line) & (response == status)

    if not mask.any():
        raise ValueError(f"No cells found for required group: {key}")

    group_masks[key] = mask

if not np.all(
    np.logical_or.reduce(list(group_masks.values()))
):
    raise ValueError("Some cells do not belong to the four expected groups.")

print("=== STATEDRUGRESPONSE GENE FILTERING ===")
print(f"Input: {INPUT}")
print(f"Input shape: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Raw counts layer: {counts.dtype}")
print("All four required groups are present.")


# --------------------------------------------------
# 3. Calculate gene detection and expression metrics
# --------------------------------------------------

def detected_cells(matrix):
    """Count cells with raw counts greater than zero per gene."""
    return np.asarray((matrix > 0).sum(axis=0)).ravel()


n_genes = adata.n_vars
gene_ids = adata.var_names.astype(str)

du145_mask = cell_line == "DU145"
pc3_mask = cell_line == "PC3"

detected_du145 = detected_cells(counts[du145_mask, :])
detected_pc3 = detected_cells(counts[pc3_mask, :])

detected_by_group = {}

for key, mask in group_masks.items():
    detected_by_group[key] = detected_cells(counts[mask, :])

# Approximation for the lowest-15% step:
# average raw counts per cell across the entire dataset.
# This is not claimed to reproduce Fluidigm Singular.
mean_raw_counts = np.asarray(counts.mean(axis=0)).ravel()

# Use a DataFrame to keep one auditable row per original gene.
audit = pd.DataFrame({
    "gene_id": gene_ids,
    "detected_cells_DU145": detected_du145,
    "detected_cells_PC3": detected_pc3,
    "mean_raw_counts_all_cells": mean_raw_counts,
})

for key, values in detected_by_group.items():
    audit[f"detected_cells_{key}"] = values


# --------------------------------------------------
# 4. Step 1: detection in both cell lines
# --------------------------------------------------

step1 = (
    (detected_du145 >= MIN_CELLS)
    & (detected_pc3 >= MIN_CELLS)
)

audit["pass_step1_both_cell_lines"] = step1

print("\nStep 1: Genes detected in at least 10 cells in both cell lines")
print(f"Genes remaining: {int(step1.sum())}")


# --------------------------------------------------
# 5. Step 2: remove the lowest 15% by mean raw counts
# --------------------------------------------------

eligible_indices = np.flatnonzero(step1)
n_eligible = len(eligible_indices)

# Number of genes to remove, using the ceiling so that
# at least 15% of eligible genes are removed.
n_remove = int(np.ceil(LOWEST_FRACTION * n_eligible))

# Sort eligible genes by mean raw count, ascending.
# Gene ID provides deterministic tie-breaking.
rank_table = audit.loc[
    eligible_indices, ["gene_id", "mean_raw_counts_all_cells"]
].sort_values(
    ["mean_raw_counts_all_cells", "gene_id"],
    ascending=[True, True],
    kind="mergesort",
)

remove_indices = rank_table.index[:n_remove]

step2 = step1.copy()
step2[remove_indices] = False

audit["pass_step2_after_lowest15pct"] = step2
audit["removed_by_lowest15pct"] = step1 & ~step2

print("\nStep 2: Remove the lowest 15% by mean raw counts")
print(f"Eligible before step: {n_eligible}")
print(f"Genes removed: {n_remove}")
print(f"Genes remaining: {int(step2.sum())}")


# --------------------------------------------------
# 6. Step 3: detection in all four groups
# --------------------------------------------------

step3 = np.logical_and.reduce([
    detected_by_group[key] >= MIN_CELLS
    for key in group_masks
])

# Apply the final group criterion only to genes that
# survived steps 1 and 2.
final_keep = step2 & step3

audit["pass_step3_all_four_groups"] = step3
audit["final_retained"] = final_keep

print("\nStep 3: Genes detected in at least 10 cells in all four groups")
print(f"Genes meeting group criterion overall: {int(step3.sum())}")
print(f"Final genes retained: {int(final_keep.sum())}")


# --------------------------------------------------
# 7. Save audit and summary tables
# --------------------------------------------------

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

summary = pd.DataFrame([
    {
        "step": 0,
        "criterion": "All genes in normalized input",
        "genes_remaining": n_genes,
    },
    {
        "step": 1,
        "criterion": "Detected in >=10 cells in both cell lines",
        "genes_remaining": int(step1.sum()),
    },
    {
        "step": 2,
        "criterion": "Remove lowest 15% by mean raw counts",
        "genes_remaining": int(step2.sum()),
    },
    {
        "step": 3,
        "criterion": "Also detected in >=10 cells in all four groups",
        "genes_remaining": int(final_keep.sum()),
    },
])

audit.to_csv(AUDIT_OUTPUT, index=False)
summary.to_csv(SUMMARY_OUTPUT, index=False)


# --------------------------------------------------
# 8. Save a separate gene-filtered dataset
# --------------------------------------------------

filtered = adata[:, final_keep].copy()

filtered.uns["gene_filtering"] = {
    "method": "Paper-inspired, auditable approximation",
    "input_file": INPUT.name,
    "minimum_detected_cells": MIN_CELLS,
    "detection_definition": "Raw count > 0",
    "step1": "At least 10 detected cells in both DU145 and PC3",
    "step2": (
        "Remove lowest 15% of step-1 genes by mean raw counts "
        "across all input cells; approximation, not verified "
        "Fluidigm Singular reproduction"
    ),
    "step3": "At least 10 detected cells in each of four groups",
    "paper_reported_final_gene_count": 12862,
    "raw_counts_layer_preserved": True,
}

filtered.write_h5ad(OUTPUT, compression="gzip")


# --------------------------------------------------
# 9. Report comparison with the paper
# --------------------------------------------------

paper_count = 12862
final_count = int(final_keep.sum())

print("\n=== FILTERING SUMMARY ===")
print(summary.to_string(index=False))

print("\n=== COMPARISON WITH PAPER ===")
print(f"Paper-reported genes: {paper_count}")
print(f"Our final genes: {final_count}")
print(f"Difference (ours - paper): {final_count - paper_count:+d}")

print("\n=== OUTPUTS ===")
print(f"Filtered dataset: {OUTPUT}")
print(f"Per-gene audit: {AUDIT_OUTPUT}")
print(f"Step summary: {SUMMARY_OUTPUT}")

print("\nExisting input datasets were not overwritten.")
