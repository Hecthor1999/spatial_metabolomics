from functions_leiden_clust import load_cut_imzml_data, cluster_leiden, create_UMAP_in_RGB
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os

def process_image(file_path, save_path, image_name):
    """
    input: file path to imzML obtained from the cardinal pipeline
    output: 
    - csv cluster file
    - png file for clustering image
    - png file for UMAP

    """

    adata = load_cut_imzml_data(file_path)

    # normalise and log transfrom
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace = True)
    sc.pp.log1p(adata)

    # cluster
    cluster_leiden(adata, resolution = 0.3, save_path = save_path, file_path = file_path)

    # Save UMAP as PNG (assuming UMAP has been computed in create_UMAP_in_RGB)
    create_UMAP_in_RGB(adata, save_path=save_path, file_path=file_path) 


if __name__ == "__main__":

    file_path = sys.argv[1]
    save_path = sys.argv[2]
    image_name = sys.argv[3]
    image = process_image(file_path, save_path, image_name) 
    print(f"Processed {file_path}")
