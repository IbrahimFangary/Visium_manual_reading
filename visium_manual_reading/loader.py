# manual_visium_loader/loader.py

import scipy.io
from anndata import AnnData
import scanpy as sc
import pandas as pd
import numpy as np
import json
from PIL import Image
import os 
import gzip

def load_manual_spatial_data(
    adata,
    spatial_dir="spatial",
    sample_id="my_sample",
    hires_image_name="tissue_hires_image.png",
    lowres_image_name="tissue_lowres_image.png",
    positions_file="tissue_positions_list.csv",
    scalefactors_file="scalefactors_json.json"
):
    """
    Load spatial coordinates and image data into an AnnData object manually,
    mimicking 10x Genomics Visium spatial structure.

    Parameters:
    - adata: AnnData object to attach spatial data to.
    - spatial_dir: Path to the directory containing spatial files.
    - sample_id: Key name for the sample in adata.uns['spatial'].
    - hires_image_name: High-resolution tissue image filename.
    - lowres_image_name: Low-resolution tissue image filename.
    - positions_file: CSV file with barcode and pixel positions.
    - scalefactors_file: JSON file with image scalefactors.

    Returns:
    - AnnData object with spatial info added.
    """
    
    # Load tissue positions
    positions_path = f"{spatial_dir}/{positions_file}"
    positions = pd.read_csv(positions_path, header=None)
    positions.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_col_in_fullres', 'pxl_row_in_fullres']
    positions = positions[positions['in_tissue'] == 1]
    positions.set_index('barcode', inplace=True)

    # Attach positions to AnnData.obs
    if not positions.index.isin(adata.obs_names).all():
        raise ValueError("tissue_positions_list.csv index must match all barcodes in the Visium data (adata.obs_names).")
    adata.obs = adata.obs.join(positions)
    
    # Store spatial coordinates (correct order: X then Y)
    adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()

    # Load scalefactors
    scalefactors_path = f"{spatial_dir}/{scalefactors_file}"
    with open(scalefactors_path) as f:
        scalefactors = json.load(f)

    # Load tissue images
    hires_img = np.array(Image.open(f"{spatial_dir}/{hires_image_name}")).astype("float32") / 255.
    lowres_img = np.array(Image.open(f"{spatial_dir}/{lowres_image_name}")).astype("float32") / 255.

    # Store in AnnData.uns
    adata.uns["spatial"] = {
        sample_id: {
            "images": {
                "hires": hires_img,
                "lowres": lowres_img
            },
            "scalefactors": scalefactors
        }
    }

    return adata
    
def load_manual_visium_data(path="filtered_count_matrix/", metadata=None):
    """
    Load Visium count matrix data manually from the given path.

    Parameters:
    - path: str
        Path to the folder containing 'matrix.mtx', 'barcodes.tsv', and 'features.tsv'.
    - metadata: pd.DataFrame, optional
        Additional metadata to add to the AnnData object (adata.obs).
        The index of the DataFrame **must match the Visium cell/spot barcodes** 
        (typically from 'barcodes.tsv') exactly.

    Returns:
    - AnnData
        AnnData object containing the spatial gene expression matrix.
    """

    # Load matrix and transpose it (genes x spots -> spots x genes)
    matrix_file = os.path.join(path, "matrix.mtx")
    if os.path.exists(matrix_file + ".gz"):
        with gzip.open(matrix_file + ".gz", "rt") as f:
            X = scipy.io.mmread(f).tocsr().T
    else:
        X = scipy.io.mmread(matrix_file).tocsr().T

    # Load barcodes
    barcodes_file = os.path.join(path, "barcodes.tsv")
    if os.path.exists(barcodes_file + ".gz"):
        obs = pd.read_csv(barcodes_file + ".gz", header=None, sep="\t")
    else:
        obs = pd.read_csv(barcodes_file, header=None, sep="\t")
    obs.columns = ["barcode"]
    obs.set_index('barcode', inplace= True)  # Use barcodes as the index

    # Load features
    features_file = os.path.join(path, "features.tsv")
    if os.path.exists(features_file + ".gz"):
        var = pd.read_csv(features_file + ".gz", header=None, sep="\t")
    else:
        var = pd.read_csv(features_file, header=None, sep="\t")
    var.columns = ['gene_ids', 'gene_names', 'feature_types']

    # Create AnnData object
    adata = AnnData(X, obs=obs, var=var)

    # Merge metadata if provided
    if metadata is not None:
        if not metadata.index.isin(adata.obs_names).all():
            raise ValueError("Metadata index must match all barcodes in the Visium data (adata.obs_names).")
        adata.obs = adata.obs.join(metadata)
    print(adata)
    return adata
