from pathlib import Path
import pandas as pd
import gseapy as gp
import scanpy as sc

BASE = Path(".")
RESULTS = BASE / "results" / "pathways"
RESULTS.mkdir(parents=True, exist_ok=True)

# Input gene lists created in the previous step.
gene_lists = {
    "resistant": RESULTS / "resistant_shared_direction_genes.txt",
    "sensitive": RESULTS / "sensitive_shared_direction_genes.txt",
}

# Use the genes retained by the current gene-filtering workflow as background.
adata = sc.read_h5ad(
    BASE / "data" / "processed" / "GSE140440_gene_filtered.h5ad",
    backed="r",
)
background = (
    adata.var["gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
)
background = sorted(set(x for x in background if x and x.lower() != "nan"))
adata.file.close()

libraries = [
    "GO_Biological_Process_2025",
    "Reactome_Pathways_2024",
]

summary_rows = []

print("Background genes:", len(background))
print("Starting pathway enrichment...")

for group, list_path in gene_lists.items():
    genes = sorted(set(
        line.strip()
        for line in list_path.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    ))

    # Restrict candidate genes to the tested background universe.
    tested_genes = sorted(set(genes) & set(background))

    print(f"\n{group.upper()}: {len(genes)} input symbols; "
          f"{len(tested_genes)} overlap with background")

    for library in libraries:
        print(f"Running {library} for {group}...")

        try:
            result = gp.enrich(
                gene_list=tested_genes,
                gene_sets=library,
                background=background,
                outdir=None,
                cutoff=1.0,
                no_plot=True,
            )

            table = result.results.copy()

            out_file = RESULTS / f"{group}_{library}_enrichment.csv"
            table.to_csv(out_file, index=False)

            if not table.empty and "Adjusted P-value" in table.columns:
                significant = (
                    pd.to_numeric(table["Adjusted P-value"], errors="coerce") < 0.05
                ).sum()
            else:
                significant = 0

            summary_rows.append({
                "group": group,
                "library": library,
                "input_genes": len(genes),
                "background_overlap": len(tested_genes),
                "terms_returned": len(table),
                "terms_adjusted_p_below_0.05": int(significant),
                "result_file": str(out_file),
                "status": "OK",
            })

            print(f"  Terms returned: {len(table)}")
            print(f"  Adjusted p-value < 0.05: {significant}")
            print(f"  Saved: {out_file}")

        except Exception as exc:
            summary_rows.append({
                "group": group,
                "library": library,
                "input_genes": len(genes),
                "background_overlap": len(tested_genes),
                "terms_returned": 0,
                "terms_adjusted_p_below_0.05": 0,
                "result_file": "",
                "status": f"ERROR: {exc}",
            })
            print(f"  ERROR: {exc}")

summary = pd.DataFrame(summary_rows)
summary.to_csv(RESULTS / "pathway_enrichment_summary.csv", index=False)

print("\n=== ENRICHMENT SUMMARY ===")
print(summary.to_string(index=False))
print("\nSummary saved to:", RESULTS / "pathway_enrichment_summary.csv")
