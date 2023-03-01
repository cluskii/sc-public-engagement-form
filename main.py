from datetime import datetime
import time

from scipy.sparse import csr_matrix
import numpy as np
from auth import authenticate_google_credentials
from extraction import extract_from_google_sheets
from anndata import AnnData
import scanpy as sc
import tkinter
import itertools


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

    cells_on_graph = []

    window = tkinter.Tk()
    window.title("PCA Demonstration")
    canvas_dimension = 500
    window.geometry(f'{canvas_dimension}x{canvas_dimension}')

    canvas = tkinter.Canvas(window)
    canvas.configure(bg="white")
    canvas.pack(fill="both", expand=True)

    last_update_time = 0.0
    creds = None
    rows = None
    previous_anndata = AnnData()

    while True:

        current_time = time.time()
        if current_time > last_update_time + 30:
            print("Fetching data again")
            last_update_time = current_time
            creds = authenticate_google_credentials()
            rows = extract_from_google_sheets(creds, "")

            # two colors for new vs old
            # Draw automatically
            if rows:
                adata = get_anndata(rows)

                if len(previous_anndata.obs.index.to_list()) == 0:
                    previous_anndata = adata

                print(adata)
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
                sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=10, min_disp=0.5)

                sc.pp.scale(adata, max_value=10)
                sc.tl.pca(adata, svd_solver='arpack')
                sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)

                all_coords = list(itertools.chain.from_iterable(adata.obsm['X_pca']))
                cell_coord_min_abs = abs(min(all_coords))
                cell_coord_max = max(all_coords)
                del all_coords

                for cell in cells_on_graph:
                    canvas.delete(cell)

                for i, cell in enumerate(adata.obsm['X_pca']):
                    cell_x = (((cell_coord_min_abs + cell[0]) / cell_coord_max) * canvas_dimension / 2) + 50
                    cell_y = (((cell_coord_min_abs + cell[1]) / cell_coord_max) * canvas_dimension / 2) + 50

                    color = 'black'
                    if adata.obs.index.to_list()[i] not in previous_anndata.obs.index.to_list():
                        color = 'red'

                    cells_on_graph.append(canvas.create_oval(cell_x - 5, cell_y - 5, cell_x + 5, cell_y + 5, fill=color))

                previous_anndata = adata
        window.update()
        time.sleep(0.1)
