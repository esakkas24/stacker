import mdtraj as md
import pandas
import numpy as np
from residue_movement import calc_center_3pts
from vector import *

def calculate_residue_distance(trajectory : md.Trajectory, 
                               res1_num : int, res2_num : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> float:
    '''Calculates the distance between two residues in Angstroms

    Calcualtes the distance between the center of two residues. The center is denoted
        by the average x,y,z position of three passed atoms for each residue (typically
        every other carbon on the 6-C ring of the nucleotide base).

    Args:
        trajectory (md.Trajectory) : single frame trajectory
        res1_num (int) : the residue number of the first residue (PDB Column 5)
        res2_num (int) : the residue number of the second residue (PDB Column 5)
        res1_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 [("C2","C4","C6")]
        res2_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 [("C2","C4","C6")]
    
    Returns:
        distance_res12 (float) : scalar distance between the center of geometry of the two residues
    '''
    # Correct for mdtraj 0-indexing
    res1_num = res1_num - 1 
    res2_num = res2_num - 1

    topology = trajectory.topology
    res1_atom_idx_query = "(name " + res1_atoms[0] + " or name " + res1_atoms[1] + " or name " + res1_atoms[2] + ") and residue " + str(res1_num)
    res1_atom_idices = topology.select(res1_atom_idx_query)
    res2_atom_idx_query = "(name " + res2_atoms[0] + " or name " + res2_atoms[1] + " or name " + res2_atoms[2] + ") and residue " + str(res2_num)
    res2_atom_idices = topology.select(res2_atom_idx_query)

    # convert nanometer units in trajectory.xyz to Angstroms
    res1_atom_xyz = trajectory.xyz[0, res1_atom_idices,:] * 10
    vectorized_res1_atom_xyz = [Vector(x,y,z) for [x,y,z] in res1_atom_xyz]
    res1_center_of_geometry = calc_center_3pts(*vectorized_res1_atom_xyz)

    res2_atom_xyz = trajectory.xyz[0, res2_atom_idices,:] * 10
    vectorized_res2_atom_xyz = [Vector(x,y,z) for [x,y,z] in res2_atom_xyz]
    res2_center_of_geometry = calc_center_3pts(*vectorized_res2_atom_xyz)

    distance_res12 = (res1_center_of_geometry - res2_center_of_geometry).magnitude()
    return distance_res12

def get_residue_distance_for_frame(trajectory : md.Trajectory, frame : int) -> np.array:
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
    # np.array[i,j] = dist(nres1, nres2)
    pass

if __name__ == "__main__":
    trj = md.load('first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', top = '5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop')
    print(calculate_residue_distance(trj[0], 426, 427))
