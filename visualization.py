import matplotlib.pyplot as plt
import residue_movement

def visualize_two_residue_movement_scatterplot(csv_filepath : str) -> None:
    '''Creates scatterplot of two-residue movement relative to each other

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
    '''Creates heatmap of two-residue movement relative to each other

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

