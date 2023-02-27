from datetime import datetime
import time

from scipy.sparse import csr_matrix
import numpy as np
from auth import authenticate_google_credentials
from extraction import extract_from_google_sheets
from anndata import AnnData
import scanpy as sc


def get_anndata(rows):
    data_start_column = 2  # I.e. column C
    column_headers = rows[0][data_start_column:]
    adata_input = []
    for row in rows[1:]:
        adata_input.append(row[data_start_column:])
    counts = csr_matrix(adata_input, dtype=np.float32)
    adata = AnnData(counts)
    adata.obs_names = [f"person_{i}" for i in range(adata.n_obs)]
    adata.var_names = [f"{column_headers[i]}" for i in range(adata.n_vars)]
    return adata


if __name__ == "__main__":
    while True:
        creds = authenticate_google_credentials()
        rows = extract_from_google_sheets(creds, "")

        if rows:
            adata = get_anndata(rows)
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=10, min_disp=0.5)
            # sc.pl.highly_variable_genes(adata, show=False)
            # adata.raw = adata
            # adata = adata[:, adata.var.highly_variable]

            # sc.pp.regress_out(adata, ['total_counts'])
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, svd_solver='arpack')
            sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
            sc.tl.umap(adata)
            sc.tl.leiden(adata, resolution=1)

            sc.pl.umap(adata, color=['leiden'], legend_loc='on data')

        time.sleep(60)
