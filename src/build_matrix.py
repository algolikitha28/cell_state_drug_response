import csv
import gzip
import io
import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parent.parent
ARCHIVE = ROOT / "data" / "raw" / "GSE140440_RAW.tar"
MANIFEST = ROOT / "data" / "processed" / "validated_sample_manifest.csv"
OUTPUT = ROOT / "data" / "processed" / "GSE140440_counts.h5ad"

manifest = pd.read_csv(MANIFEST, dtype=str).fillna("")

required = {"archive_filename", "cell_line", "docetaxel_response"}
missing = required - set(manifest.columns)
if missing:
    raise ValueError(
        f"Manifest is missing columns: {sorted(missing)}. "
        f"Available columns: {manifest.columns.tolist()}"
    )

manifest = manifest.sort_values("archive_filename").reset_index(drop=True)

data = []
indices = []
indptr = [0]
gene_ids = None
gene_symbols = None

print(f"Building matrix from {len(manifest)} samples...")

with tarfile.open(ARCHIVE, "r:*") as tar:
    members = {
        Path(m.name).name: m
        for m in tar.getmembers()
        if m.isfile()
    }

    for i, row in manifest.iterrows():
        filename = row["archive_filename"]
        member = members.get(filename)

        if member is None:
            raise FileNotFoundError(
                f"Archive member not found: {filename}"
            )

        raw = tar.extractfile(member)
        if raw is None:
            raise OSError(f"Could not read {filename}")

        current_ids = []
        current_symbols = []
        counts = []

        with gzip.GzipFile(fileobj=raw, mode="rb") as gz:
            text = io.TextIOWrapper(gz, encoding="utf-8-sig")
            reader = csv.reader(text, delimiter="\t")

            for fields in reader:
                if not fields:
                    continue

                if len(fields) != 2:
                    raise ValueError(
                        f"Unexpected row format in {filename}: {fields}"
                    )

                feature = fields[0].strip()
                count = int(fields[1])

                if count < 0:
                    raise ValueError(
                        f"Negative count in {filename}: {feature}"
                    )

                parts = feature.split("|", 1)
                current_ids.append(parts[0])
                current_symbols.append(
                    parts[1] if len(parts) > 1 else ""
                )
                counts.append(count)

        current_ids = np.asarray(current_ids, dtype=str)
        counts = np.asarray(counts, dtype=np.int32)

        if gene_ids is None:
            gene_ids = current_ids
            gene_symbols = np.asarray(current_symbols, dtype=str)
        elif not np.array_equal(gene_ids, current_ids):
            raise ValueError(
                f"Gene order differs from reference in {filename}"
            )

        nonzero = np.flatnonzero(counts)
        indices.extend(nonzero.tolist())
        data.extend(counts[nonzero].tolist())
        indptr.append(len(data))

        if (i + 1) % 25 == 0 or i + 1 == len(manifest):
            print(f"Processed {i + 1}/{len(manifest)} samples")

X = csr_matrix(
    (
        np.asarray(data, dtype=np.int32),
        np.asarray(indices, dtype=np.int32),
        np.asarray(indptr, dtype=np.int64),
    ),
    shape=(len(manifest), len(gene_ids)),
    dtype=np.int32,
)

accessions = [
    re.search(r"(GSM\d+)", name).group(1)
    for name in manifest["archive_filename"]
]

obs = manifest.copy()
obs.index = pd.Index(accessions, name="sample_id")

var = pd.DataFrame(
    {
        "gene_id": gene_ids,
        "gene_symbol": gene_symbols,
        "mitochondrial": np.char.startswith(gene_symbols, "MT-"),
    },
    index=pd.Index(gene_ids, name="gene_id_index"),
)

adata = ad.AnnData(X=X, obs=obs, var=var)
adata.uns["data_description"] = (
    "Raw gene counts from GEO GSE140440; "
    "no normalization or filtering applied."
)

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUTPUT, compression="gzip")

print("\n=== MATRIX BUILD COMPLETE ===")
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"Nonzero counts: {adata.X.nnz}")
print(f"Cell-line counts:\n{adata.obs['cell_line'].value_counts().to_string()}")
print(f"Response counts:\n{adata.obs['docetaxel_response'].value_counts().to_string()}")
print(f"Mitochondrial genes: {int(adata.var['mitochondrial'].sum())}")
print(f"Saved to: {OUTPUT}")
