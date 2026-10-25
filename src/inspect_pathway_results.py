from pathlib import Path
import pandas as pd

folder = Path("results/pathways")

files = [
    "resistant_GO_Biological_Process_2025_enrichment.csv",
    "resistant_Reactome_Pathways_2024_enrichment.csv",
    "sensitive_GO_Biological_Process_2025_enrichment.csv",
    "sensitive_Reactome_Pathways_2024_enrichment.csv",
]

for filename in files:
    path = folder / filename
    df = pd.read_csv(path)

    print("\n" + "=" * 90)
    print(filename)
    print("Columns:", df.columns.tolist())

    if "Adjusted P-value" not in df.columns:
        print("Adjusted p-value column not found.")
        continue

    df["Adjusted P-value"] = pd.to_numeric(
        df["Adjusted P-value"], errors="coerce"
    )
    significant = df[df["Adjusted P-value"] < 0.05].copy()
    significant = significant.sort_values("Adjusted P-value")

    if significant.empty:
        print("No terms have adjusted p-value < 0.05.")
        continue

    columns = [
        c for c in [
            "Term", "Overlap", "P-value",
            "Adjusted P-value", "Genes"
        ] if c in significant.columns
    ]

    print("Significant terms:", len(significant))
    print(significant[columns].head(15).to_string(index=False))
