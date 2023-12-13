import numpy as np
from numpy import typing
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

class NoResidues(Exception):
    pass

def create_axis_labels(res_indicies : typing.ArrayLike, tick_distance : int = 10) -> list:
    '''Designates the axis labels to use in the pairwise plot
    
    Returns the x-axis tick positions and labels to use on the ticks based on the 
        residues used in a specific pairwise analysis. Meant to be used when many 
        disjoint sets of residue indices are used. Ticks will be present every tick_distance 
        residues in a collection of adjacent residues, and a tick will exist at both
        ends of any consecutive residue sequence.

    Args:
        res_indicies : list
            The list of residue indices used in the pairwise analysis.
        tick_distance : int, default = 10
            Distance between ticks in blocks of consecutive residues
    Returns:
        tick_locations : array_like
            List of tick positions (0-based) to place labels on in the axes
        tick_labels : array_like
            List of labels to place at the adjacent tick locations

    Examples:
    >>> create_axis_labels([0,1,2,3,4,5,6,7,8,9,10,11,12,98,99,100])
    [0,10,12,13,15], [0,10,12,98,100]

    >>> create_axis_labels([94,95,96,97,98,99,100,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428])
    [0,6,7,17,27], [94,100,408,418,428]
    '''
    n_residues = len(res_indicies)

    if n_residues < 1: raise NoResidues("pairwise analysis must include at least one residue")

    tick_locations = [0]
    tick_labels = [res_indicies[0]]

    res_sequence_length = 1

    for i in range(1, n_residues):
        if res_indicies[i] == res_indicies[i-1] + 1:
            res_sequence_length += 1

        if res_indicies[i] > res_indicies[i-1] + 1:
            tick_locations += [i-1, i]
            tick_labels += [res_indicies[i-1], res_indicies[i]]
            res_sequence_length = 1
        elif res_sequence_length == tick_distance+1:
            tick_locations += [i]
            tick_labels += [res_indicies[i]]
            res_sequence_length = 1
    
    if n_residues-1 not in tick_locations:
        tick_locations += [n_residues-1]
        tick_labels += [res_indicies[n_residues-1]]
    
    return tick_locations, tick_labels
        
def display_arrays_as_video(numpy_arrays : list | typing.ArrayLike, res_indicies : typing.ArrayLike, seconds_per_frame : int = 10, tick_distance : int = 10) -> None:
    '''Displays list/array of 2D NumPy arrays as matrix heatmaps

    Takes list/array of 2D NumPy arrays and treats them as frames 
    in a video, filling in a grid at position i,j by the value 
    at i,j in the array.

    Args:
        numpy_arrays: array_like
            List or array of 2D NumPy arrays
        res_indicies : list
            The list of residue indices used in the pairwise analysis.
        seconds_per_frame : int, default = 10
            Number of seconds to display each matrix for.
        tick_distance : int, default = 10
            Distance between ticks in blocks of consecutive residues
    Returns:
        None
        Displays video of NumPy arrays
    '''
    orange_colormap = mpl.colormaps['Oranges_r'].resampled(100)

    newcolors = np.vstack((orange_colormap(np.linspace(1, 0, 128)),
                       orange_colormap(np.linspace(0, 1, 128))))
    newcmp = mpl.colors.ListedColormap(newcolors, name='OrangeBellcurve')
    
    fig , ax = plt.subplots(figsize=(8,8))
    plt.ion()
    for hist in numpy_arrays:
        ax.clear()
        neg = ax.imshow(hist, cmap = newcmp, vmin=2, vmax=5, interpolation = 'nearest')
        ax.set_title('Distance Between Residues Center of Geometries')
        colorbar = fig.colorbar(neg, ax=ax, location='right', anchor=(0, 0.3), shrink=0.7)
        ticks, labels = create_axis_labels(res_indicies)
        plt.xticks(ticks, labels, rotation = 'vertical')
        plt.yticks(ticks, labels)
        ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
        plt.pause(seconds_per_frame)
        colorbar.remove()

def visualize_two_residue_movement_scatterplot(csv_filepath : str) -> None:
    '''Creates scatterplot of two-residue movement relative to each other.

    Args:
        csv_filepath (str) : filepath to csv file containing data on the movement
            of two residues relative to each other (r, rho, and theta values). Created
            in residue_movement
    Returns:
        None
    
    takes the data created in residue_movement and visualizes it as a polar coordinate
        scatterplot similar to the Figure D link in Proposal Feature 4.
    '''
    bottaro_values = pd.read_csv(csv_filepath, sep=',')

    theta_values = bottaro_values['theta']

    # convert rho values from nm to Angstroms
    rho_values = bottaro_values['rho_dist'] * 10

    theta_values_rad = np.radians(theta_values)

    fig = plt.figure()
    ax = fig.add_subplot(polar=True)
    ax.scatter(theta_values_rad, rho_values, s=1)
    ax.set_ylim(0,15)
    ax.set_xticks(np.pi/180. * np.linspace(0,  360, 3, endpoint=False))
    ax.set_xticklabels([r"$\theta=0^\circ$",r"$\theta=120^\circ$",r"$\theta=240^\circ$"])
    ax.set_yticks(np.linspace(0,  15, 4, endpoint=True))

    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.show()

def visualize_two_residue_movement_heatmap(csv_filepath : str) -> None:
    '''Creates heatmap of two-residue movement relative to each other.

    Args:
        csv_filepath (str) : filepath to csv file containing data on the movement
            of two residues relative to each other (r, rho, and theta values). Created
            in residue_movement
    Returns:
        None
    
    Should run visualize_two_residue_movement_scatterplot() to set the current pyplot as
        a scatterplot of the data and then adjust this scatterplot to be a heatmap instead.
        Should now look identical to Figure D link in the proposal.
    '''
    pass

if __name__ == '__main__':
    visualize_two_residue_movement_scatterplot('tUAG_aCUA_+1GCU_GC_plot.csv')