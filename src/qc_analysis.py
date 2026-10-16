from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent.parent

INPUT = ROOT / "data" / "processed" / "GSE140440_counts.h5ad"
QC_DIR = ROOT / "results" / "qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

print("Loading expression matrix...")
adata = ad.read_h5ad(INPUT)

print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# Confirm that mitochondrial gene annotations are available.
if "mitochondrial" not in adata.var.columns:
    raise ValueError("Mitochondrial gene annotations are missing.")

print(
    "Mitochondrial genes:",
    int(adata.var["mitochondrial"].sum())
)

# Calculate QC metrics without modifying or filtering the raw counts.
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mitochondrial"],
    percent_top=None,
    log1p=True,
    inplace=True,
)

# Save per-cell QC metrics together with sample metadata.
metric_columns = [
    "total_counts",
    "n_genes_by_counts",
    "total_counts_mitochondrial",
    "pct_counts_mitochondrial",
    "log1p_total_counts",
    "log1p_n_genes_by_counts",
]

available_columns = [
    col for col in metric_columns if col in adata.obs.columns
]

qc_table = adata.obs[
    ["cell_line", "docetaxel_response"] + available_columns
].copy()

qc_table.to_csv(QC_DIR / "cell_qc_metrics.csv")

# Overall descriptive statistics.
summary = qc_table[available_columns].describe().T
summary.to_csv(QC_DIR / "qc_summary.csv")

print("\n=== OVERALL QC SUMMARY ===")
print(summary.to_string())

# Summaries within each cell line and response group.
group_summary = (
    qc_table
    .groupby(["cell_line", "docetaxel_response"])[
        ["total_counts", "n_genes_by_counts",
         "pct_counts_mitochondrial"]
    ]
    .agg(["median", "mean", "min", "max"])
)

group_summary.to_csv(QC_DIR / "qc_group_summary.csv")

print("\n=== QC BY CELL LINE AND RESPONSE ===")
print(group_summary.to_string())

# Plot distributions so filtering decisions can be evidence-based.
plot_metrics = [
    ("total_counts", "Total counts per cell"),
    ("n_genes_by_counts", "Detected genes per cell"),
    ("pct_counts_mitochondrial", "Mitochondrial counts (%)"),
]

for column, title in plot_metrics:
    fig, ax = plt.subplots(figsize=(8, 5))

    for group, frame in qc_table.groupby(
        ["cell_line", "docetaxel_response"]
    ):
        label = f"{group[0]} - {group[1]}"
        ax.hist(
            frame[column].dropna(),
            bins=30,
            alpha=0.45,
            label=label,
        )

    ax.set_title(title)
    ax.set_xlabel(title)
    ax.set_ylabel("Number of cells")
    ax.legend(fontsize=8)
    fig.tight_layout()

    filename = column + "_distribution.png"
    fig.savefig(QC_DIR / filename, dpi=150)
    plt.close(fig)

print("\n=== QC COMPLETE ===")
print(f"Per-cell metrics: {QC_DIR / 'cell_qc_metrics.csv'}")
print(f"Overall summary: {QC_DIR / 'qc_summary.csv'}")
print(f"Group summary: {QC_DIR / 'qc_group_summary.csv'}")
print(f"Plots saved in: {QC_DIR}")
print("\nNo cells were filtered; raw counts remain unchanged.")
