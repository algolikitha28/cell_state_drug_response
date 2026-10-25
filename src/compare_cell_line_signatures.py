from pathlib import Path
import pandas as pd
import numpy as np

input_dir = Path("results/differential_expression")
output_dir = Path("results/evidence")
output_dir.mkdir(parents=True, exist_ok=True)

du145 = pd.read_csv(
    input_dir / "DU145_resistant_vs_sensitive_descriptive.csv"
)
pc3 = pd.read_csv(
    input_dir / "PC3_resistant_vs_sensitive_descriptive.csv"
)

# Use Ensembl gene IDs as the join key.
du145 = du145[
    ["gene", "gene_symbol", "log2_ratio_of_group_means", "direction"]
].rename(columns={
    "gene_symbol": "gene_symbol_DU145",
    "log2_ratio_of_group_means": "effect_DU145",
    "direction": "direction_DU145",
})

pc3 = pc3[
    ["gene", "gene_symbol", "log2_ratio_of_group_means", "direction"]
].rename(columns={
    "gene_symbol": "gene_symbol_PC3",
    "log2_ratio_of_group_means": "effect_PC3",
    "direction": "direction_PC3",
})

merged = du145.merge(
    pc3,
    on="gene",
    how="inner",
    validate="one_to_one"
)

merged["same_direction"] = np.sign(
    merged["effect_DU145"]
) == np.sign(merged["effect_PC3"])

merged["shared_direction"] = np.select(
    [
        (merged["effect_DU145"] > 0) & (merged["effect_PC3"] > 0),
        (merged["effect_DU145"] < 0) & (merged["effect_PC3"] < 0),
    ],
    [
        "Higher in resistant in both",
        "Higher in sensitive in both",
    ],
    default="Mixed or zero effect"
)

merged["minimum_absolute_effect"] = np.minimum(
    merged["effect_DU145"].abs(),
    merged["effect_PC3"].abs()
)

merged["mean_absolute_effect"] = (
    merged["effect_DU145"].abs()
    + merged["effect_PC3"].abs()
) / 2

merged.to_csv(
    output_dir / "cross_cell_line_descriptive_comparison.csv",
    index=False
)

shared = merged[merged["same_direction"]].copy()
shared = shared.sort_values(
    "minimum_absolute_effect",
    ascending=False
)

shared.to_csv(
    output_dir / "shared_direction_candidates.csv",
    index=False
)

print("Genes present in both comparisons:", len(merged))
print("Same direction in both cell lines:", len(shared))
print(
    "Higher in resistant in both:",
    int((shared["shared_direction"] == "Higher in resistant in both").sum())
)
print(
    "Higher in sensitive in both:",
    int((shared["shared_direction"] == "Higher in sensitive in both").sum())
)
print(
    "Different or zero direction:",
    int((~merged["same_direction"]).sum())
)

print("\nTop 15 shared-direction genes:")
print(
    shared[
        [
            "gene",
            "gene_symbol_DU145",
            "effect_DU145",
            "effect_PC3",
            "shared_direction",
            "minimum_absolute_effect",
        ]
    ].head(15).to_string(index=False)
)

print("\nSaved:")
print(output_dir / "cross_cell_line_descriptive_comparison.csv")
print(output_dir / "shared_direction_candidates.csv")
