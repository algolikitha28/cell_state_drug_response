from pathlib import Path
import pandas as pd

folder = Path("results/pathways")
files = [
    "resistant_GO_Biological_Process_2025_enrichment.csv",
    "resistant_Reactome_Pathways_2024_enrichment.csv",
    "sensitive_GO_Biological_Process_2025_enrichment.csv",
    "sensitive_Reactome_Pathways_2024_enrichment.csv",
]

reports = []

for filename in files:
    df = pd.read_csv(folder / filename)

    df["Adjusted P-value"] = pd.to_numeric(
        df["Adjusted P-value"], errors="coerce"
    )
    df["P-value"] = pd.to_numeric(df["P-value"], errors="coerce")

    sig = df.loc[df["Adjusted P-value"] < 0.05].copy()
    if sig.empty:
        print(f"\n{filename}: no significant terms.")
        continue

    sig["overlap_genes"] = sig["Genes"].fillna("").apply(
        lambda x: [g for g in x.split(";") if g]
    )
    sig["overlap_count"] = sig["overlap_genes"].apply(len)
    sig["example_genes"] = sig["overlap_genes"].apply(
        lambda genes: ";".join(genes[:8])
    )
    sig["group"] = "resistant" if filename.startswith("resistant_") else "sensitive"
    sig["library"] = (
        "GO Biological Process"
        if "_GO_" in filename else "Reactome"
    )

    keep = [
        "group", "library", "Term", "P-value", "Adjusted P-value",
        "overlap_count", "example_genes"
    ]
    reports.append(sig[keep])

    print(f"\n=== {filename} ===")
    print(
        sig.sort_values("Adjusted P-value")[keep[2:]].head(10)
        .to_string(index=False)
    )

if reports:
    report = pd.concat(reports, ignore_index=True)
    output = folder / "significant_pathways_concise.csv"
    report.to_csv(output, index=False)
    print("\nSaved concise report:", output)
    print("Total significant term rows:", len(report))
