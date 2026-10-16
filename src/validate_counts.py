import csv
import gzip
import io
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = ROOT / "data" / "raw" / "GSE140440_RAW.tar"
MANIFEST = ROOT / "data" / "processed" / "validated_sample_manifest.csv"
OUTPUT = ROOT / "data" / "processed" / "count_file_validation.csv"

manifest = pd.read_csv(MANIFEST, dtype=str)
expected_files = set(manifest["archive_filename"].dropna())

results = []
reference_genes = None
reference_file = None

with tarfile.open(ARCHIVE, "r:*") as tar:
    members = {
        Path(member.name).name: member
        for member in tar.getmembers()
        if member.isfile()
        and Path(member.name).name in expected_files
    }

    for filename in sorted(expected_files):
        result = {
            "filename": filename,
            "status": "FAIL",
            "gene_count": 0,
            "duplicate_gene_ids": 0,
            "invalid_rows": 0,
            "negative_counts": 0,
            "non_integer_counts": 0,
            "gene_set_matches_reference": False,
            "gene_order_matches_reference": False,
            "error": "",
        }

        try:
            member = members.get(filename)
            if member is None:
                raise ValueError("File not found in archive")

            raw = tar.extractfile(member)
            if raw is None:
                raise ValueError("Could not read archive member")

            with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
                text = io.TextIOWrapper(gz, encoding="utf-8-sig")
                reader = csv.reader(text, delimiter="\t")

                gene_ids = []
                invalid_rows = 0
                negative_counts = 0
                non_integer_counts = 0

                for row in reader:
                    if not row or all(not field.strip() for field in row):
                        continue

                    if len(row) != 2:
                        invalid_rows += 1
                        continue

                    gene_id = row[0].strip()
                    count_text = row[1].strip()

                    if not gene_id:
                        invalid_rows += 1
                        continue

                    gene_ids.append(gene_id)

                    try:
                        count = int(count_text)
                        if count < 0:
                            negative_counts += 1
                    except ValueError:
                        non_integer_counts += 1

                if not gene_ids:
                    raise ValueError("No valid gene rows found")

                duplicate_gene_ids = len(gene_ids) - len(set(gene_ids))

                gene_set = set(gene_ids)
                gene_order = tuple(gene_ids)

                if reference_genes is None:
                    reference_genes = gene_set
                    reference_file = filename

                result.update({
                    "gene_count": len(gene_ids),
                    "duplicate_gene_ids": duplicate_gene_ids,
                    "invalid_rows": invalid_rows,
                    "negative_counts": negative_counts,
                    "non_integer_counts": non_integer_counts,
                    "gene_set_matches_reference": (
                        gene_set == reference_genes
                    ),
                    "gene_order_matches_reference": (
                        gene_order == tuple(
                            gene for gene in gene_order
                            if gene in reference_genes
                        )
                        and gene_order == reference_order
                    ) if "reference_order" in globals() else False,
                })

                if "reference_order" not in globals():
                    reference_order = gene_order
                    result["gene_order_matches_reference"] = True

                passed = (
                    invalid_rows == 0
                    and negative_counts == 0
                    and non_integer_counts == 0
                    and duplicate_gene_ids == 0
                    and gene_set == reference_genes
                )

                result["status"] = "PASS" if passed else "CHECK"
                if not passed:
                    result["error"] = (
                        "Review row format, counts, duplicates, or gene set"
                    )

        except Exception as exc:
            result["error"] = str(exc)

        results.append(result)

validation = pd.DataFrame(results)
validation.to_csv(OUTPUT, index=False)

print("\n=== GENE-COUNT FILE VALIDATION ===")
print(f"Files expected: {len(expected_files)}")
print(f"Files found in archive: {len(members)}")
print(f"Reference gene file: {reference_file}")
print(f"Reference gene count: {len(reference_genes) if reference_genes else 0}")
print("\nStatus counts:")
print(validation["status"].value_counts(dropna=False).to_string())

print("\nGene-count summary:")
print(validation["gene_count"].describe().to_string())

print("\nFiles needing review:")
issues = validation[validation["status"] != "PASS"]
if issues.empty:
    print("None")
else:
    print(
        issues[
            [
                "filename",
                "status",
                "gene_count",
                "duplicate_gene_ids",
                "invalid_rows",
                "negative_counts",
                "non_integer_counts",
                "gene_set_matches_reference",
                "error",
            ]
        ].head(20).to_string(index=False)
    )

print(f"\nDetailed results saved to: {OUTPUT}")
