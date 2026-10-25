from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


INPUT_FILE = Path("data/processed/GSE140440_gene_filtered.h5ad")
OUTPUT_DIR = Path("results/differential_expression")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CELL_LINES = ["DU145", "PC3"]
PSEUDOCOUNT = 0.1


def mean_expression(adata, mask):
    """Return mean normalized expression for the selected samples."""
    X = adata.X[mask]

    if hasattr(X, "toarray"):
        X = X.toarray()

    return np.asarray(X).mean(axis=0).ravel()


def main():
    adata = sc.read_h5ad(INPUT_FILE)

    required = {"cell_line", "docetaxel_response"}
    missing = required - set(adata.obs.columns)

    if missing:
        raise ValueError(f"Missing required metadata: {sorted(missing)}")

    if not adata.var_names.is_unique:
        raise ValueError("Gene identifiers are not unique.")

    if not np.isfinite(adata.X.data).all():
        raise ValueError("Expression matrix contains non-finite values.")

    if "counts" not in adata.layers:
        raise ValueError("Raw counts layer is missing.")

    all_results = []

    for cell_line in CELL_LINES:
        line_mask = adata.obs["cell_line"].eq(cell_line).to_numpy()
        response = adata.obs["docetaxel_response"]

        resistant_mask = (
            line_mask & response.eq("Resistant").to_numpy()
        )
        sensitive_mask = (
            line_mask & response.eq("Sensitive").to_numpy()
        )

        n_resistant = int(resistant_mask.sum())
        n_sensitive = int(sensitive_mask.sum())

        if n_resistant == 0 or n_sensitive == 0:
            raise ValueError(
                f"{cell_line}: both response groups must contain samples."
            )

        resistant_mean = mean_expression(adata, resistant_mask)
        sensitive_mean = mean_expression(adata, sensitive_mask)

        # These are descriptive differences in the stored normalized,
        # log-transformed expression values, not inferential test results.
        log2_fold_change = np.log2(
            (resistant_mean + PSEUDOCOUNT)
            / (sensitive_mean + PSEUDOCOUNT)
        )

        results = pd.DataFrame({
            "gene": adata.var_names.astype(str),
            "cell_line": cell_line,
            "n_resistant": n_resistant,
            "n_sensitive": n_sensitive,
            "mean_resistant_normalized": resistant_mean,
            "mean_sensitive_normalized": sensitive_mean,
            "log2_ratio_of_group_means": log2_fold_change,
            "direction": np.where(
                log2_fold_change > 0,
                "Higher in resistant",
                np.where(
                    log2_fold_change < 0,
                    "Higher in sensitive",
                    "Equal"
                )
            ),
        })

        results["absolute_effect"] = results[
            "log2_ratio_of_group_means"
        ].abs()

        results = results.sort_values(
            "absolute_effect", ascending=False
        )

        output_file = OUTPUT_DIR / (
            f"{cell_line}_resistant_vs_sensitive_descriptive.csv"
        )
        results.to_csv(output_file, index=False)
        all_results.append(results)

        print(f"\n=== {cell_line}: RESISTANT VS SENSITIVE ===")
        print(f"Resistant samples: {n_resistant}")
        print(f"Sensitive samples: {n_sensitive}")
        print("Genes compared:", len(results))
        print(
            "Higher in resistant:",
            int((results["log2_ratio_of_group_means"] > 0).sum())
        )
        print(
            "Higher in sensitive:",
            int((results["log2_ratio_of_group_means"] < 0).sum())
        )
        print("\nTop 10 genes by absolute descriptive effect:")
        print(
            results[
                [
                    "gene",
                    "mean_resistant_normalized",
                    "mean_sensitive_normalized",
                    "log2_ratio_of_group_means",
                    "direction",
                ]
            ].head(10).to_string(index=False)
        )
        print("Saved:", output_file)

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(
        OUTPUT_DIR / "combined_descriptive_results.csv",
        index=False
    )
    print("\nCombined results saved.")
    print("Output directory:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
