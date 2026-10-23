
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


# --------------------------------------------------
# 1. Define project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT = (
    ROOT
    / "data"
    / "processed"
    / "GSE140440_counts.h5ad"
)

OUTPUT = (
    ROOT
    / "data"
    / "processed"
    / "GSE140440_qc_filtered.h5ad"
)

QC_DIR = ROOT / "results" / "qc"

AUDIT_OUTPUT = QC_DIR / "qc_filter_audit.csv"
SUMMARY_OUTPUT = QC_DIR / "qc_filter_group_summary.csv"


# --------------------------------------------------
# 2. Define provisional QC policy
# --------------------------------------------------

# These two cells are proposed for exclusion because
# their QC metrics show multiple warning signs.
PROPOSED_EXCLUSIONS = {
    "GSM4161522",
    "GSM4161689",
}

# This is a documented, dataset-specific decision.
# It is not a universal rule for other datasets.


# --------------------------------------------------
# 3. Load the original count matrix
# --------------------------------------------------

if not INPUT.exists():
    raise FileNotFoundError(
        f"Input dataset not found: {INPUT}"
    )

adata = ad.read_h5ad(INPUT)

print("=== STATE DRUG RESPONSE: QC FILTERING ===")
print(f"Input: {INPUT}")
print(f"Original cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")


# --------------------------------------------------
# 4. Validate required metadata
# --------------------------------------------------

required_obs = {
    "cell_line",
    "docetaxel_response",
}

missing_obs = required_obs - set(adata.obs.columns)

if missing_obs:
    raise ValueError(
        f"Missing cell metadata: {sorted(missing_obs)}"
    )

if not adata.obs_names.is_unique:
    raise ValueError("Cell IDs are not unique.")

if "mitochondrial" not in adata.var.columns:
    raise ValueError(
        "Mitochondrial gene annotation is missing."
    )

if not PROPOSED_EXCLUSIONS.issubset(
    set(adata.obs_names.astype(str))
):
    missing_ids = (
        PROPOSED_EXCLUSIONS
        - set(adata.obs_names.astype(str))
    )
    raise ValueError(
        f"Proposed exclusion IDs not found: {sorted(missing_ids)}"
    )


# --------------------------------------------------
# 5. Calculate QC metrics
# --------------------------------------------------

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mitochondrial"],
    percent_top=None,
    log1p=True,
    inplace=True,
)

print("QC metrics calculated.")


# --------------------------------------------------
# 6. Create a cell-level audit table
# --------------------------------------------------

audit = adata.obs[
    ["cell_line", "docetaxel_response"]
].copy()

audit.index.name = "sample_id"

metrics = [
    "total_counts",
    "n_genes_by_counts",
    "total_counts_mitochondrial",
    "pct_counts_mitochondrial",
]

for column in metrics:
    audit[column] = adata.obs[column]

audit["proposed_exclusion"] = (
    audit.index.astype(str).isin(PROPOSED_EXCLUSIONS)
)

audit["decision"] = audit["proposed_exclusion"].map(
    {
        True: "Exclude: provisional QC decision",
        False: "Retain",
    }
)

audit["exclusion_reason"] = ""

audit.loc[
    audit.index.astype(str) == "GSM4161522",
    "exclusion_reason",
] = (
    "Very low total counts and detected genes, "
    "with elevated mitochondrial percentage"
)

audit.loc[
    audit.index.astype(str) == "GSM4161689",
    "exclusion_reason",
] = (
    "Extremely high mitochondrial percentage "
    "and low detected gene count"
)


# --------------------------------------------------
# 7. Create a separate filtered AnnData object
# --------------------------------------------------

keep_mask = ~audit["proposed_exclusion"].to_numpy()

filtered = adata[keep_mask, :].copy()

# Record the QC policy in the output metadata.
filtered.uns["qc_policy"] = {
    "policy_name": "GSE140440 provisional QC exclusions",
    "input_file": INPUT.name,
    "excluded_sample_ids": sorted(PROPOSED_EXCLUSIONS),
    "excluded_count": len(PROPOSED_EXCLUSIONS),
    "criteria": (
        "Manual review of multiple QC metrics; "
        "not a universal threshold"
    ),
    "raw_counts_preserved": True,
}


# --------------------------------------------------
# 8. Summarize retained and excluded cells by group
# --------------------------------------------------

group_rows = []

for (cell_line, response), group in audit.groupby(
    ["cell_line", "docetaxel_response"],
    sort=True,
):
    original_count = len(group)
    excluded_count = int(group["proposed_exclusion"].sum())
    retained_count = original_count - excluded_count

    group_rows.append(
        {
            "cell_line": cell_line,
            "docetaxel_response": response,
            "original_cells": original_count,
            "retained_cells": retained_count,
            "excluded_cells": excluded_count,
            "retained_percent": round(
                100 * retained_count / original_count,
                2,
            ),
        }
    )

summary = pd.DataFrame(group_rows)


# --------------------------------------------------
# 9. Save output files
# --------------------------------------------------

QC_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)

audit.to_csv(AUDIT_OUTPUT)
summary.to_csv(SUMMARY_OUTPUT, index=False)

filtered.write_h5ad(OUTPUT, compression="gzip")


# --------------------------------------------------
# 10. Validate and report
# --------------------------------------------------

print("\n=== QC FILTER SUMMARY ===")
print(summary.to_string(index=False))

print(f"\nOriginal cells: {adata.n_obs}")
print(f"Excluded cells: {adata.n_obs - filtered.n_obs}")
print(f"Retained cells: {filtered.n_obs}")
print(f"Genes retained: {filtered.n_vars}")

print(f"\nFiltered dataset: {OUTPUT}")
print(f"Cell audit: {AUDIT_OUTPUT}")
print(f"Group summary: {SUMMARY_OUTPUT}")

print("\nOriginal input file was not overwritten.")
