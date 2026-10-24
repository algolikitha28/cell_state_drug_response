
import pandas as pd
from pathlib import Path

# Load QC metrics
qc = pd.read_csv("results/qc/cell_qc_metrics.csv")

# PCA outliers identified in the previous analysis
outlier_ids = [
    "GSM4161453", "GSM4161500", "GSM4161523",
    "GSM4161496", "GSM4161478", "GSM4161499",
    "GSM4161512", "GSM4161466", "GSM4161480",
    "GSM4161594",
    "GSM4161681", "GSM4161714", "GSM4161632",
    "GSM4161756", "GSM4161732", "GSM4161753",
    "GSM4161755", "GSM4161691", "GSM4161768",
    "GSM4161774"
]

# Keep only the cells identified as PCA outliers
cols = [
    "sample_id",
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial"
]

outliers = qc[qc["sample_id"].astype(str).isin(outlier_ids)][cols].copy()

# Sort by cell line and response
outliers = outliers.sort_values(
    ["cell_line", "docetaxel_response", "total_counts"]
)

print("\n=== PCA OUTLIERS: QC METRICS ===")
print(outliers.to_string(index=False))

# Save results for later review
output = Path("results/eda/pca_outliers_qc_comparison.csv")
output.parent.mkdir(parents=True, exist_ok=True)
outliers.to_csv(output, index=False)

print(f"\nSaved comparison to: {output}")
print(f"Matched {len(outliers)} of {len(outlier_ids)} requested cells.")
