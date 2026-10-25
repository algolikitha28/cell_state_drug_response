from pathlib import Path
import pandas as pd

INPUT_FILE = Path("results/evidence/shared_direction_candidates.csv")
OUTPUT_DIR = Path("results/pathways")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(INPUT_FILE)

# Check that the required columns are available.
required = ["gene", "gene_symbol_DU145", "gene_symbol_PC3", "shared_direction"]
missing = [column for column in required if column not in df.columns]

if missing:
    raise ValueError(f"Missing required columns: {missing}")

# Both cell lines should have the same symbol for each Ensembl gene.
symbol_match = df["gene_symbol_DU145"].eq(df["gene_symbol_PC3"])
if not symbol_match.all():
    raise ValueError("Some gene symbols disagree between cell lines.")

# Use one symbol column for the shared-direction analysis.
df["gene_symbol"] = df["gene_symbol_DU145"].astype(str).str.strip()

# Exclude blank symbols and common missing-value strings.
invalid = df["gene_symbol"].isin(["", "nan", "None", "NA"])
df = df.loc[~invalid].copy()

# Save one list per direction.
groups = {
    "resistant": "Higher in resistant in both",
    "sensitive": "Higher in sensitive in both",
}

summary = []

for name, direction in groups.items():
    subset = df.loc[df["shared_direction"] == direction].copy()

    # Retain one record per symbol. This avoids submitting the same
    # symbol repeatedly to an enrichment tool.
    unique = subset.drop_duplicates(subset="gene_symbol", keep="first")
    genes = sorted(unique["gene_symbol"].unique())

    output_file = OUTPUT_DIR / f"{name}_shared_direction_genes.txt"
    output_file.write_text("\n".join(genes) + "\n", encoding="utf-8")

    summary.append({
        "group": name,
        "input_rows": len(subset),
        "unique_symbols": len(genes),
        "duplicate_symbol_rows_removed": len(subset) - len(genes),
        "gene_list_file": str(output_file),
    })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(OUTPUT_DIR / "gene_list_preparation_summary.csv", index=False)

print("=== GENE LIST PREPARATION ===")
print(summary_df.to_string(index=False))
print("\nSaved files:")
print(OUTPUT_DIR / "resistant_shared_direction_genes.txt")
print(OUTPUT_DIR / "sensitive_shared_direction_genes.txt")
print(OUTPUT_DIR / "gene_list_preparation_summary.csv")
