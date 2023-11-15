import mdtraj as md
import pandas
import numpy as np

def calculate_residue_distance(res1_num : int, res2_num : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> float:
    '''Calculates the distance between two residues in Angstroms

    Calcualtes the distance between the center of two residues. The center is denoted
        by the average x,y,z position of three passed atoms for each residue.

    Args:
        res1_num (int) : the residue number of the first residue (PDB Column 5)
        res2_num (int) : the residue number of the second residue (PDB Column 5)
        res1_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 [("C2","C4","C6")]
        res2_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 [("C2","C4","C6")]
    '''
    pass

def get_residue_distance_for_frame(trajectory  : md.Trajectory, frame : int) -> np.array:
    '''Calculates pairwise the distance between all residues in a given frame

    Args:
        trajectory (md.Trajectory) : trajectory to analyze (must have topology aspect)
        frame (int) : frame to analyze
    
    Returns:
        distances (2D NumPy array) : matrix where position i,j represents the distance from 
            residue i to residue j
    
    Only need to implement for the top right triangle of the matrix since distance i -> j
        is the same as distance j -> i. 
    '''
    pass

if __name__ == "__main__":
    ########JOB VARIABLES#######
    #load pdb file
    topology_filename = '5JUP_N2_wCCC_+1GCU_nowat_posttleap.pdb'

    trajectory_filename = '5JUP_N2_wCCC_+1GCU_nowat.mdcrd'
    ############################

    ########OPTIONAL VARS#######
    perspective_atom1_name = "C2"
    perspective_atom2_name = "C4"
    perspective_atom3_name = "C6"
    viewed_atom1_name = "C2"
    viewed_atom2_name = "C4"
    viewed_atom3_name = "C6"
    ############################
