from pathlib import Path
import pandas as pd

BASE = Path(".")
INPUT_DIR = BASE / "results" / "differential_expression"
OUTPUT_DIR = BASE / "results" / "pathways" / "per_cell_line"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

cell_lines = {
    "DU145": INPUT_DIR / "DU145_resistant_vs_sensitive_descriptive.csv",
    "PC3": INPUT_DIR / "PC3_resistant_vs_sensitive_descriptive.csv",
}

required = {"gene", "gene_symbol", "direction"}
summary = []

for cell_line, input_file in cell_lines.items():
    df = pd.read_csv(input_file)

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{input_file}: missing columns {sorted(missing)}")

    if df["gene"].duplicated().any():
        raise ValueError(f"{input_file}: duplicate gene IDs found")

    for group, direction in [
        ("resistant", "Higher in resistant"),
        ("sensitive", "Higher in sensitive"),
    ]:
        subset = df.loc[df["direction"] == direction].copy()

        symbols = (
            subset["gene_symbol"]
            .dropna()
            .astype(str)
            .str.strip()
        )

        invalid = {"", "nan", "none", "na"}
        genes = sorted({
            symbol for symbol in symbols
            if symbol.lower() not in invalid
        })

        output_file = OUTPUT_DIR / f"{cell_line}_{group}_genes.txt"
        output_file.write_text(
            "\n".join(genes) + "\n" if genes else "",
            encoding="utf-8",
        )

        summary.append({
            "cell_line": cell_line,
            "group": group,
            "rows_in_direction": len(subset),
            "unique_gene_symbols": len(genes),
            "output_file": str(output_file),
        })

summary_df = pd.DataFrame(summary)
summary_file = OUTPUT_DIR / "gene_list_summary.csv"
summary_df.to_csv(summary_file, index=False)

print("=== PER-CELL-LINE GENE LIST SUMMARY ===")
print(summary_df.to_string(index=False))
print("\nSaved summary:", summary_file)
