import matplotlib.pyplot as plt
import residue_movement
import pairwise_distance
import numpy as np

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
    pass

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

def visualize_pairwise_one_frame(pairwise_matrix : np.typing.ArrayLike) -> None:
    '''Creates matrix heatmap of the matrix created in the pairwise comparison.

    Args:
        pairwise_matrix (NumPy array) : matrix where position i,j represents the distance from 
            residue i to residue j. Output of pairwise_distance.get_residue_distance_for_frame()

    Returns:
        None

    Implementation very similar to game of life where a cell is colored by the presence of a 1 value
        in matrix position i,j. Here values can range >1 and we care about values 0-8 Angstroms, as 
        this is where stacking occurs.
    '''
    pass

def visualize_pairwise_trajectory(heatmap_frames : list) -> None:
    '''Creates movie of heatmaps across all frames.

    Displays the pairwise heatmaps for every frame of a trajectory to visualize
        stacking fingerprint overtime.

    Args:
        heatmap_frames (list) : list of NumPy arrays where each array is matrix where position i,j 
            represents the distance from residue i to residue j. Each is the output of 
            pairwise_distance.get_residue_distance_for_frame(). 
    
    Returns:
        None

    heatmap_frames[i] is the heatmap created at frame i+1 (since frames are 1-indexed). 
    Implementation is similar to the updating game of life.
    '''
    pass
