from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "results" / "qc" / "cell_qc_metrics.csv"
OUTPUT = ROOT / "results" / "qc" / "qc_threshold_comparison.csv"

df = pd.read_csv(INPUT)

# Alternative thresholds for sensitivity analysis.
rules = {
    "No filtering": lambda x: pd.Series(True, index=x.index),
    "Counts >= 200000": lambda x: x["total_counts"] >= 200000,
    "Genes >= 5000": lambda x: x["n_genes_by_counts"] >= 5000,
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
}

rows = []

for rule_name, rule in rules.items():
    keep = rule(df)

    for (cell_line, response), group in df.groupby(
        ["cell_line", "docetaxel_response"]
    ):
        retained = int(keep.loc[group.index].sum())

        rows.append({
            "rule": rule_name,
            "cell_line": cell_line,
            "response": response,
            "original_cells": len(group),
            "retained_cells": retained,
            "removed_cells": len(group) - retained,
            "retained_percent": round(
                100 * retained / len(group), 2
            ),
        })

result = pd.DataFrame(rows)
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

print(f"\nSaved comparison to: {OUTPUT}")
print("This analysis does not remove any cells.")
