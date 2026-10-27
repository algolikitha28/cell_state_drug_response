from pathlib import Path
import pandas as pd
import scanpy as sc
import gseapy as gp

BASE = Path(".")
GENE_DIR = BASE / "results" / "pathways" / "per_cell_line"
OUTPUT_DIR = GENE_DIR / "enrichment"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(
    BASE / "data" / "processed" / "GSE140440_gene_filtered.h5ad",
    backed="r",
)

background = set(
    adata.var["gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
)
background = {
    gene for gene in background
    if gene and gene.lower() not in {"nan", "none", "na"}
}
adata.file.close()

libraries = [
    "GO_Biological_Process_2025",
    "Reactome_Pathways_2024",
]

summary_rows = []

for gene_file in sorted(GENE_DIR.glob("*_genes.txt")):
    genes = {
        line.strip()
        for line in gene_file.read_text(encoding="utf-8-sig").splitlines()
        if line.strip()
    }
    tested_genes = sorted(genes & background)

    print(f"\n=== {gene_file.stem} ===")
    print(f"Input genes: {len(genes)}")
    print(f"Genes in background: {len(tested_genes)}")

    for library in libraries:
        print(f"Running {library}...")

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
            output_file = OUTPUT_DIR / (
                f"{gene_file.stem}_{library}_enrichment.csv"
            )
            table.to_csv(output_file, index=False)

            if (
                not table.empty
                and "Adjusted P-value" in table.columns
            ):
                adjusted_p = pd.to_numeric(
                    table["Adjusted P-value"],
                    errors="coerce",
                )
                significant = int((adjusted_p < 0.05).sum())
            else:
                significant = 0

            summary_rows.append({
                "gene_list": gene_file.stem,
                "library": library,
                "input_genes": len(genes),
                "background_overlap": len(tested_genes),
                "terms_returned": len(table),
                "significant_terms": significant,
                "status": "OK",
                "result_file": str(output_file),
            })

            print(f"Terms returned: {len(table)}")
            print(f"Adjusted p < 0.05: {significant}")
            print(f"Saved: {output_file}")

        except Exception as exc:
            summary_rows.append({
                "gene_list": gene_file.stem,
                "library": library,
                "input_genes": len(genes),
                "background_overlap": len(tested_genes),
                "terms_returned": 0,
                "significant_terms": 0,
                "status": f"ERROR: {exc}",
                "result_file": "",
            })
            print(f"ERROR: {exc}")

summary = pd.DataFrame(summary_rows)
summary_file = OUTPUT_DIR / "per_cell_line_enrichment_summary.csv"
summary.to_csv(summary_file, index=False)

print("\n=== ENRICHMENT SUMMARY ===")
print(summary.to_string(index=False))
print("\nSaved:", summary_file)
