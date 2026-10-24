
from pathlib import Path
import pandas as pd

# --------------------------------------------------
# 1. Define project paths
# --------------------------------------------------

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "results" / "qc" / "cell_qc_metrics.csv"
OUTPUT = ROOT / "results" / "qc" / "qc_threshold_comparison.csv"

# Load QC metrics. This script does not modify the raw matrix.
df = pd.read_csv(INPUT)

# Check that the required columns are present.
required_columns = {
    "sample_id",
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
}

missing = required_columns - set(df.columns)

if missing:
    raise ValueError(f"Missing required columns: {sorted(missing)}")

if df["sample_id"].duplicated().any():
    raise ValueError("Duplicate sample IDs found in QC metrics.")

# --------------------------------------------------
# 2. Define provisional exclusion scenarios
# --------------------------------------------------

# Stronger candidates based on multiple QC concerns.
proposed_exclusions = {
    "GSM4161522",
    "GSM4161689",
}

# Expanded scenario for sensitivity analysis.
expanded_exclusions = proposed_exclusions | {
    "GSM4161523",
    "GSM4161659",
    "GSM4161466",
}

# --------------------------------------------------
# 3. Define filtering rules
# --------------------------------------------------

rules = {
    "No filtering": lambda x: pd.Series(True, index=x.index),

    "Counts >= 200000": lambda x: (
        x["total_counts"] >= 200000
    ),

    "Genes >= 5000": lambda x: (
        x["n_genes_by_counts"] >= 5000
    ),

    "Mitochondrial < 20%": lambda x: (
        x["pct_counts_mitochondrial"] < 20
    ),

    "Counts >= 200000 AND genes >= 5000": lambda x: (
        (x["total_counts"] >= 200000)
        & (x["n_genes_by_counts"] >= 5000)
    ),

    "All three thresholds": lambda x: (
        (x["total_counts"] >= 200000)
        & (x["n_genes_by_counts"] >= 5000)
        & (x["pct_counts_mitochondrial"] < 20)
    ),

    "Proposed exclusions (2 cells)": lambda x: (
        ~x["sample_id"].isin(proposed_exclusions)
    ),

    "Expanded proposed exclusions (5 cells)": lambda x: (
        ~x["sample_id"].isin(expanded_exclusions)
    ),
}

# --------------------------------------------------
# 4. Compare retained cells in each group
# --------------------------------------------------

rows = []

for rule_name, rule in rules.items():
    keep = rule(df)

    for (cell_line, response), group in df.groupby(
        ["cell_line", "docetaxel_response"]
    ):
        retained = int(keep.loc[group.index].sum())
        original = len(group)

        rows.append({
            "rule": rule_name,
            "cell_line": cell_line,
            "response": response,
            "original_cells": original,
            "retained_cells": retained,
            "removed_cells": original - retained,
            "retained_percent": round(
                100 * retained / original, 2
            ),
        })

result = pd.DataFrame(rows)

# --------------------------------------------------
# 5. Save results
# --------------------------------------------------

result.to_csv(OUTPUT, index=False)

print("=== QC THRESHOLD SENSITIVITY ANALYSIS ===")

for rule_name in rules:
    print(f"\n--- {rule_name} ---")

    subset = result[result["rule"] == rule_name]

    print(
        subset[
            [
                "cell_line",
                "response",
                "original_cells",
                "retained_cells",
                "removed_cells",
                "retained_percent",
            ]
        ].to_string(index=False)
    )

    print(
        "Total retained:",
        int(subset["retained_cells"].sum()),
        "/",
        int(subset["original_cells"].sum()),
    )

print(f"\nSaved comparison to: {OUTPUT}")
print("No cells were removed from the original dataset.")
