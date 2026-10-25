import scanpy as sc
import numpy as np
from scipy import sparse

original = sc.read_h5ad("data/processed/GSE140440_counts.h5ad")
qc = sc.read_h5ad("data/processed/GSE140440_qc_filtered.h5ad")
genes = sc.read_h5ad("data/processed/GSE140440_gene_filtered.h5ad")

print("=== SAMPLE ALIGNMENT ===")
print("QC sample IDs match gene-filtered IDs:",
      qc.obs_names.equals(genes.obs_names))

print("\n=== GENE ALIGNMENT ===")
missing_genes = genes.var_names.difference(qc.var_names)
print("Genes in gene-filtered file missing from QC file:",
      len(missing_genes))

if len(missing_genes) == 0:
    gene_positions = qc.var_names.get_indexer(genes.var_names)

    qc_counts = qc.X[:, gene_positions]
    saved_counts = genes.layers["counts"]

    if sparse.issparse(qc_counts):
        qc_counts = qc_counts.tocsr()
    else:
        qc_counts = np.asarray(qc_counts)

    if sparse.issparse(saved_counts):
        saved_counts = saved_counts.tocsr()
    else:
        saved_counts = np.asarray(saved_counts)

    difference = qc_counts - saved_counts

    if sparse.issparse(difference):
        difference.eliminate_zeros()
        matches = difference.nnz == 0
    else:
        matches = np.array_equal(difference, np.zeros_like(difference))

    print("\n=== RAW COUNT VALIDATION ===")
    print("Gene-filtered counts exactly match QC values:", matches)
else:
    print("Cannot compare counts until gene alignment is resolved.")

print("\n=== SAMPLE IDENTIFIERS ===")
print("Original sample IDs unique:", original.obs_names.is_unique)
print("QC sample IDs unique:", qc.obs_names.is_unique)
print("Gene-filtered sample IDs unique:", genes.obs_names.is_unique)
