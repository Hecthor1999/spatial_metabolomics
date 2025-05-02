import os, glob, ast
import scanpy as sc
import seaborn as sns
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.sparse import csr_matrix
from matplotlib.colors import to_rgb
import anndata as ad
from pyimzml.ImzMLParser import ImzMLParser
from matplotlib.colors import rgb2hex


def cluster_leiden(adata, n_iterations = 2, resolution = 0.5, spot_size = 1, show_plots = False, save_path = None, file_path = None):
    """
    Perform leiden clustering on an AnnData object.
    The argument spot_size describes the resolution in um / how far apart the coordinates are.
    """
    sc.pp.pca(adata, n_comps = 50)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, n_components = 3)
    sc.tl.leiden(adata, key_added = "clusters", flavor = "leidenalg", directed = False, n_iterations = n_iterations, resolution = resolution)
    
    # sc.pl.spatial(adata, color = "clusters", img_key = None, spot_size = spot_size)
    sc.pl.umap(adata, color = ["clusters"])

    # get spatial coordinates
    x = adata.obsm['spatial'][:, 0]
    y = adata.obsm['spatial'][:, 1]

    # get cluster labels and convert them to numbers
    clusters = adata.obs["clusters"].astype(str).astype("category")
    cluster_labels = clusters.cat.codes  # Convert categories to numerical codes

    # generate a colormap using Seaborn or Matplotlib
    num_clusters = len(clusters.cat.categories)
    palette = sns.color_palette("tab10", num_clusters)  
    colors = np.array([palette[i] for i in cluster_labels])  

    # scatter plot with square markers
    plt.figure(figsize=(6, 6))
    plt.scatter(x, y, c=colors, marker="s", s=100)  
    plt.axis("equal") 
    plt.gca().invert_yaxis()

      # Save outputs
    if save_path and file_path:
        os.makedirs(save_path, exist_ok=True)
        image_name = os.path.splitext(os.path.basename(file_path))[0]

        # Save image
        out_file = os.path.join(save_path, f"{image_name}_umap.png")
        try:
            plt.savefig(out_file, dpi=300)
            print(f"Saved plot to: {out_file}")
        except Exception as e:
            print("Failed to save image:", e)

        # Save CSV
        cluster_colors = [rgb2hex(palette[i]) for i in cluster_labels]
        rgb_colors = [palette[i] for i in cluster_labels]

        # Save full cluster info
        cluster_df = pd.DataFrame({
            "x": x.astype(int),
            "y": y.astype(int),
            "cluster": clusters.astype(str),
            "cluster_color": cluster_colors,
            "RGB_color": rgb_colors
        })

        csv_file = os.path.join(save_path, f"{image_name}_leiden_clusters.csv")
        try:
            cluster_df.to_csv(csv_file, index=False)
            print(f"Saved cluster CSV to: {csv_file}")
        except Exception as e:
            print("Failed to save CSV:", e)

    if show_plots:
        plt.show()
    plt.close()

def create_UMAP_in_RGB(adata, show_plots = True, save_path=None, file_path = None):
    """
    Rescale the values of the 3 axes of the UMAP to be shown as RGB colours.
    """

    # rescale using the individual max and min for each colour axis (global max and min gives less varied colours)
    old_min = np.min(adata.obsm['X_umap'], axis = 0)
    old_max = np.max(adata.obsm['X_umap'], axis = 0)
    new_min = 0
    new_max = 1
    umap_RGB = (adata.obsm['X_umap'] - old_min) / (old_max - old_min) * (new_max - new_min)

    # convert the rows to tuples defining a RGB colour
    umap_RGB = [tuple(row) for row in umap_RGB]

    # add the UMAP in RGB scale as observation to adata
    adata.obs['umap_RGB'] = umap_RGB

    # plot
    plt.figure(figsize=(6,6))
    plt.scatter(x = adata.obsm['spatial'][:,0], y = adata.obsm['spatial'][:,1], c = umap_RGB, marker="s", s=50)
    plt.gca().invert_yaxis()

    if save_path and file_path:
        os.makedirs(save_path, exist_ok=True)
        image_name = os.path.splitext(os.path.basename(file_path))[0]
        out_file = os.path.join(save_path, f"{image_name}_umap_rgb.png")
        plt.savefig(out_file, dpi=300)

    if show_plots:
        plt.show()
    plt.close()

def filter_m_over_z(adata, min_m_z = 0, max_m_z = 10000):
    """
    Filter the peaks for a certain m/z cutoff. The defaults are chosen so they are outside of the range of m/z in this data.
    """

    adata_m_z = adata[:, (adata.var['m_over_z'] > min_m_z) & (adata.var['m_over_z'] < max_m_z)].copy()

def subcluster_islets(adata, islets_cluster_id):
    """
    Subcluster the pancreatic islet cluster to find morphology
    """
    islets_cluster_id = str(islets_cluster_id)
    sc.tl.leiden(adata, resolution = 0.4, restrict_to = ('clusters',[islets_cluster_id]), key_added ='islets_subcluster')
    
    #sc.pl.spatial(adata, color = "islets_subcluster", img_key = None, spot_size = 1)
    #sc.pl.umap(adata, color = ["islets_subcluster"])

    x = adata.obsm['spatial'][:, 0]
    y = adata.obsm['spatial'][:, 1]
    # Get cluster labels and convert them to numbers
    clusters = adata.obs["islets_subcluster"].astype(str).astype("category")
    cluster_labels = clusters.cat.codes  # Convert categories to numerical codes
    num_clusters = len(clusters.cat.categories)
    palette = sns.color_palette("tab10", num_clusters)  
    colors = np.array([palette[i] for i in cluster_labels])  
    # Scatter plot with square markers
    plt.figure(figsize=(6, 6))
    scatter = plt.scatter(x, y, c=colors, marker="s", s=50) 
    
    # Add a legend
    plt.gca().invert_yaxis()
    plt.axis("equal")  # Ensures proper aspect ratio
    plt.title(f"Subclustering of Islet Cluster {islets_cluster_id}, 030622_B71407")
    plt.show()

    # UMAP Plot
    sc.pl.umap(adata, color=["islets_subcluster"])
    

def load_cut_imzml_data(imzml_path: str):
    # Parse the imzML file
    parser = ImzMLParser(imzml_path)

    # Create a list to store the spectra and coordinates
    my_spectra = []

    for idx, (x, y, z) in enumerate(parser.coordinates):
        mzs, intensities = parser.getspectrum(idx)
        my_spectra.append([mzs, intensities, (x, y, z)])

    my_spectra = np.array(my_spectra, dtype="object")
    coordinates = np.array([spectrum[2] for spectrum in my_spectra])

    # determine bounds
    x_vals = coordinates[:, 0]
    y_vals = coordinates[:, 1]

    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()

    # Filter out spectra on the outer edges
    valid_indices = [
        idx for idx, (x, y, z) in enumerate(coordinates)
        if x_min < x < x_max and y_min < y < y_max
    ]

    my_spectra = my_spectra[valid_indices]

    # Extract intensity values and pad
    intensities = [spectrum[1] for spectrum in my_spectra]
    max_len = max(len(spec) for spec in intensities)
    padded_intensities = [np.pad(spec, (0, max_len - len(spec))) for spec in intensities]

    intensity_matrix = csr_matrix(padded_intensities)

    # Create AnnData object
    adata = sc.AnnData(X=intensity_matrix)

    # Add spatial coordinates
    adata.obs['x'] = [coord[0] for _, _, coord in my_spectra]
    adata.obs['y'] = [coord[1] for _, _, coord in my_spectra]
    adata.obs['z'] = [coord[2] for _, _, coord in my_spectra]
    adata.obsm['spatial'] = adata.obs[['x', 'y']].to_numpy()

    # Add m/z values as variables
    mz_values = [spectrum[0] for spectrum in my_spectra]
    adata.var['mz'] = np.array(mz_values[0])  # assuming consistent m/z axis

    return adata
