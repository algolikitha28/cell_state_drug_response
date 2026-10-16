from pathlib import Path

import anndata as ad
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

MATRIX = ROOT / "data" / "processed" / "GSE140440_counts.h5ad"
QC_FILE = ROOT / "results" / "qc" / "cell_qc_metrics.csv"
OUTPUT = ROOT / "results" / "qc" / "qc_review_candidates.csv"

adata = ad.read_h5ad(MATRIX)
qc = pd.read_csv(QC_FILE, index_col=0)

# Use the original GSM accession stored in the AnnData index.
qc.index = qc.index.astype(str)
qc.index.name = "sample_id"

# Add mitochondrial gene annotation and the sample's QC measurements.
qc["cell_line"] = adata.obs["cell_line"].astype(str).to_numpy()
qc["docetaxel_response"] = (
    adata.obs["docetaxel_response"].astype(str).to_numpy()
)

# Flag cells for review, not automatic exclusion.
qc["high_mito_flag"] = qc["pct_counts_mitochondrial"] >= 20
qc["low_gene_flag"] = qc["n_genes_by_counts"] < 5000
qc["low_count_flag"] = qc["total_counts"] < 200000

candidates = qc[
    qc["high_mito_flag"]
    | qc["low_gene_flag"]
    | qc["low_count_flag"]
].copy()

columns = [
    "cell_line",
    "docetaxel_response",
    "total_counts",
    "n_genes_by_counts",
    "pct_counts_mitochondrial",
    "high_mito_flag",
    "low_gene_flag",
    "low_count_flag",
]

candidates[columns].to_csv(OUTPUT)

print("=== QC CANDIDATES: REVIEW ONLY ===")
print(candidates[columns].to_string())
print(f"\nCandidates flagged: {len(candidates)}")
print(f"Saved to: {OUTPUT}")

print("\n=== MITOCHONDRIAL GENE ANNOTATION ===")
mito = adata.var.loc[adata.var["mitochondrial"]]
print(f"Mitochondrial genes: {len(mito)}")
print(mito[["gene_id", "gene_symbol"]].to_string(index=False))

print("\nNo cells have been removed.")
