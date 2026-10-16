from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
INPUT = ROOT / "results" / "qc" / "cell_qc_metrics.csv"

df = pd.read_csv(INPUT)

columns = [
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
]

print("=== TOP 15 MITOCHONDRIAL PERCENTAGES ===")
print(
    df.nlargest(15, "pct_counts_mitochondrial")[columns]
    .to_string(index=False)
)

print("\n=== LOWEST 15 TOTAL COUNTS ===")
print(
    df.nsmallest(15, "total_counts")[columns]
    .to_string(index=False)
)

print("\n=== LOWEST 15 DETECTED-GENE COUNTS ===")
print(
    df.nsmallest(15, "n_genes_by_counts")[columns]
    .to_string(index=False)
)

print("\n=== CELLS WITH MITOCHONDRIAL PERCENTAGE >= 20% ===")
high_mito = df[df["pct_counts_mitochondrial"] >= 20]

print(f"Number of cells: {len(high_mito)}")

if not high_mito.empty:
    print(high_mito[columns].sort_values(
        "pct_counts_mitochondrial", ascending=False
    ).to_string(index=False))

print("\n=== CELLS WITH FEWER THAN 5,000 DETECTED GENES ===")
low_genes = df[df["n_genes_by_counts"] < 5000]

print(f"Number of cells: {len(low_genes)}")

if not low_genes.empty:
    print(low_genes[columns].sort_values(
        "n_genes_by_counts"
    ).to_string(index=False))

print("\nNo cells have been removed.")
