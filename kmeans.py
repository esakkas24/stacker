from sklearn.cluster import KMeans
import numpy as np
from numpy import typing
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

### VARIABLES ###

N_RESIDUES = 127
N_CLUSTERS = 3
dataset_names = ['tGGG_aCCU_+1GCU', 'tGGG_aCCU_+1CGU', 'tUAG_aCUA_+1GCU', 'tUAG_aCUA_+1CGU', 'tGGG_aCCC_+1GCU', 'tGGG_aCCC_+1CGU']  # Add more dataset names as needed
indir = '/home66/esakkas/STACKER/DATA/' # Directory with data.txt output from StACKER (created with -d flag)
outdir = '/home66/esakkas/STACKER/DATA/' # Outdir for clustering results and kmeans plot

##################


def read_and_preprocess_data(dataset_names, data_path, n_residues) -> dict:
    '''Reads and preprocesses data for each dataset

    Reads data from text files for each dataset, decompresses the data, 
    and attaches each Trajectory to its frame-wise SSF results. The values
    are flattened SSF lists, so rather than a 3200 frames x 127 res x 127 res, 
    it's a 3200 frames x 16129 res-res pairs. For example, an SSF of 

    [ [[1, 2],
       [3, 4]],
       
       [[5, 6],
       [7, 8]] ]

    has been flattened to:

    [[1, 2, 3, 4],
     [5, 6, 7, 8]]

    Args:
        dataset_names : list of str
            List of dataset names to read and preprocess.
        data_path : str
            Base path to the directory containing the data files.
        n_residues : int
            Number of residues in each frame of the data.

    Returns:
        data_arrays : dict
            Dictionary where keys are dataset names and values are the processed data arrays.

    Example:
        >>> dataset_names = ['dataset1', 'dataset2'] # 3200 frames, SSFs of 127 x 127 residues
        >>> data_path = '/path/to/data/'
        >>> n_residues = 127
        >>> data_arrays = read_and_preprocess_data(dataset_names, data_path, n_residues)
        >>> print(data_arrays['dataset1'].shape)
        (3200, 16129)
    '''
    data_arrays = {}
    for name in dataset_names:
        print('Reading data:', name)
        data = np.loadtxt(f'{data_path}5JUP_N2_{name}_data.txt')
        data_arrays[name] = data
    return data_arrays

def create_kmeans_input(data_arrays : dict) -> typing.ArrayLike:   
    '''Stacks SSF data into a single 2D Numpy array.

    Creates a 2D numpy array from all frames of all trajectories without
    labels for each frame.

    Args:
        data_arrays : dict
            output of read_and_preprocess_data().Dictionary where keys are dataset 
            names and values are the processed data arrays.

    Returns:
        np.typing.ArrayLike
            A 2D numpy array containing all frames stacked together.

    Example:
        >>> data_arrays = {
        ...     'dataset1': np.random.rand(3200, 16129),
        ...     'dataset2': np.random.rand(3200, 16129)
        ... }
        >>> kmeans_input = create_kmeans_input(data_arrays)
        >>> print(kmeans_input.shape)
        (6400, 16129)
    '''
    blind_data = list(data_arrays.values())
    data = np.vstack(blind_data)
    print(data.shape)
    return data

def run_kmeans(blinded_data : typing.ArrayLike, N_CLUSTERS : int = N_CLUSTERS,
                max_iter : int = 1000, n_init : int =20, random_state : int =1, outdir : str = ''):
    '''Performs KMeans clustering on the provided data and saves the results.

    This function applies the KMeans clustering algorithm to the provided
    blinded data, assigns each frame to a cluster, and counts the number of
    frames in each cluster for each dataset. The results are printed and
    saved to a file.

    Args:
        blinded_data : np.typing.ArrayLike
            A 2D numpy array containing all frames stacked together.
        N_CLUSTERS : int
            The number of clusters to form (default is N_CLUSTERS).
        max_iter : int, default = 1000
            Maximum number of iterations of the k-means algorithm for a single run
        n_init : int, default = 20
            Number of time the k-means algorithm will be run with different centroid seeds 
        random_state : int, default = 1
            Determines random number generation for centroid initialization 

    Returns:
        None

    Example:
        >>> data_arrays = {
        ...     'dataset1': np.random.rand(3200, 16129),
        ...     'dataset2': np.random.rand(3200, 16129)
        ... }
        >>> blinded_data = create_kmeans_input(data_arrays)
        >>> run_kmeans(blinded_data, N_CLUSTERS=4)
        Reading data: dataset1
        Reading data: dataset2
        (6400, 16129)
        {'dataset1': array([800, 800, 800, 800]), 'dataset2': array([800, 800, 800, 800])}
        Dataset: dataset1
            Cluster 1: 800 matrices
            Cluster 2: 800 matrices
            Cluster 3: 800 matrices
            Cluster 4: 800 matrices
        Dataset: dataset2
            Cluster 1: 800 matrices
            Cluster 2: 800 matrices
            Cluster 3: 800 matrices
            Cluster 4: 800 matrices

    '''
    kmeans_func_instance = KMeans(n_clusters=N_CLUSTERS, max_iter=max_iter, n_init=n_init, random_state=random_state)
    kmeans_func_instance.fit(blinded_data)
    blindframes_labelled_by_cluster = kmeans_func_instance.labels_

    counts = {}
    for name, arr in data_arrays.items():
        counts[name] = np.bincount(blindframes_labelled_by_cluster[:len(arr)], minlength=N_CLUSTERS)
        blindframes_labelled_by_cluster = blindframes_labelled_by_cluster[len(arr):]  # Move to the next dataset

    # Print the results
    for name, count in counts.items():
        print(f'Dataset: {name}')
        for cluster in range(N_CLUSTERS):
            print(f'\tCluster {cluster+1}: {count[cluster]} matrices')


    # Save results to file
    if outdir:
        outfile = open(outdir + 'clustering_results_' + str(N_CLUSTERS) + '.txt', 'w')
        outfile.write('cluster trj number\n')
        for name, count in counts.items():
            for cluster in range(N_CLUSTERS):
                outfile.write(f'{cluster+1} {name} {count[cluster]}\n')
        outfile.close()

def plot_cluster_trj_data(input_file : str, n_cluster: int, outfile : str, seeded: bool = False) -> None:
    '''Plots the output of run_kmeans() to a PNG file
    
    Creates a grouped bar plot of the number of frames from each trajectory in each cluster
    following KMeans clustering. Writes the plot output to a PNG file.

    Args:
        input_file : str
             Path to the input file containing clustering results.
        n_cluster : int
             The number of clusters in this KMeans clustering sequence
        outfile : str
             Filepath where the output PNG file will be saved.
        seeded : bool, default = False
             Set to True if a random seed was assigned during clustering using kseed, False otherwise.
             Used to label plot title, no other impact.
        

    Returns:
        None

    Example:
        >>> plot_cluster_trj_data('clustering_results.txt', 4, True, '/path/to/output/')
        # This will read the clustering results from 'clustering_results.txt',
        # create a bar plot, and save it as 'kmeans_plot.cluster_4.png' in the specified output directory.
    '''
    cluster_data = pd.read_table(input_file, sep=' ', header=0, quotechar="\"")
    
    g = sns.FacetGrid(cluster_data.dropna(subset=['trj']), col="cluster", col_wrap=2, height = 6)

    colors = sns.color_palette("husl", len(cluster_data['trj'].unique()))

    g.map(plt.bar, 'trj', 'number', color=colors) 

    for ax in g.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)
            label.set_ha('right')
            
    g.set_titles(col_template="{col_name}")
    plt.tight_layout()
    plt.savefig(outfile)
    print(f"Plot Outputted to {outfile}")

if __name__ == "__main__":
    data_arrays = read_and_preprocess_data(dataset_names, indir, N_RESIDUES)
    blinded_data = create_kmeans_input(data_arrays)
    run_kmeans(blinded_data, N_CLUSTERS, outdir = outdir)
    cluster_file = outdir + 'clustering_results_' + str(N_CLUSTERS) + '.txt'
    outfile = f"{outdir}kmeans_plot.cluster_{N_CLUSTERS}.png"
    plot_cluster_trj_data(cluster_file, n_cluster = N_CLUSTERS, seeded = False, outfile = outfile)