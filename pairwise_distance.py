import mdtraj as md
import pandas
import numpy as np
from numpy import typing
from residue_movement import calc_center_3pts
from vector import *

import matplotlib.pyplot as plt

def display_arrays_as_video(numpy_arrays : list | typing.ArrayLike) -> None:
    '''Displays list/array of NumPy arrays as video

    Takes list/array of 2D NumPy arrays and treats them as frames 
    in a video, filling in a grid at position i,j by the value 
    at i,j in the array.

    Args:
        numpy_arrays: list or array of NumPy arrays
    Returns:
        None
        Displays video of NumPy arrays
    '''
    _ , ax = plt.subplots()
    plt.ion()
    for hist in numpy_arrays:
        ax.clear()
        ax.imshow(hist, cmap = 'Greys', vmin=0.01, vmax=4.5)
        ax.set_title('The Game of Life')
        plt.pause(10)

class MultiFrameTraj(Exception):
    pass

_NUCLEOTIDE_NAMES = {"A", "G", "C", "T", "U"}

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
    if len(trajectory) > 1: raise MultiFrameTraj("calculate_residue_distance() expects a 1-frame trajectory")

    # Correct for mdtraj 0-indexing
    res1_num = res1_num - 1 
    res2_num = res2_num - 1

    topology = trajectory.topology
    res1_name = topology.residue(res1_num).name # 
    res2_name = topology.residue(res2_num).name
    if (res1_name not in _NUCLEOTIDE_NAMES) or (res2_name not in _NUCLEOTIDE_NAMES):
        return -1

    res1_atom_idx_query = "(name " + res1_atoms[0] + " or name " + res1_atoms[1] + " or name " + res1_atoms[2] + ") and residue " + str(res1_num)
    res2_atom_idx_query = "(name " + res2_atoms[0] + " or name " + res2_atoms[1] + " or name " + res2_atoms[2] + ") and residue " + str(res2_num)
    res1_atom_idices = topology.select(res1_atom_idx_query)
    res2_atom_idices = topology.select(res2_atom_idx_query)

    # convert nanometer units in trajectory.xyz to Angstroms
    res1_atom_xyz = trajectory.xyz[0, res1_atom_idices,:] * 10
    res2_atom_xyz = trajectory.xyz[0, res2_atom_idices,:] * 10
    vectorized_res1_atom_xyz = [Vector(x,y,z) for [x,y,z] in res1_atom_xyz]
    vectorized_res2_atom_xyz = [Vector(x,y,z) for [x,y,z] in res2_atom_xyz]
    res1_center_of_geometry = calc_center_3pts(*vectorized_res1_atom_xyz)
    res2_center_of_geometry = calc_center_3pts(*vectorized_res2_atom_xyz)

    distance_res12 = (res1_center_of_geometry - res2_center_of_geometry).magnitude()
    return distance_res12

def get_residue_distance_for_frame(trajectory : md.Trajectory, frame : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> typing.ArrayLike:
    '''Calculates pairwise the distance between all residues in a given frame

    Args:
        trajectory (md.Trajectory) : trajectory to analyze (must have topology aspect)
        frame (int) : 1-indexed frame to analyze
        res1_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 [("C2","C4","C6")]
        res2_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 [("C2","C4","C6")]
    
    Returns:
        pairwise_distances (2D NumPy array) : matrix where position i,j represents the distance from 
            residue i to residue j
    
    Only need to implement for the top right triangle of the matrix since distance i -> j
        is the same as distance j -> i. 
    '''
    trajectory = trajectory[frame-1]
    n_residues = trajectory.n_residues

    pairwise_distances = np.zeros([n_residues,n_residues])

    for i in range(0, n_residues):
        print(i)
        for j in range(0, n_residues):
            if i == j: 
                pairwise_distances[i,j] = 0
            elif pairwise_distances[j,i] != 0:
                pairwise_distances[i,j] = pairwise_distances[j,i]
            else:
                pairwise_distances[i,j] = calculate_residue_distance(trajectory, i+1, j+1, res1_atoms, res2_atoms)
    
    return(pairwise_distances)

if __name__ == "__main__":
    # Load test trajectory and topology
    trj = md.load('first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', top = '5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop')
    trj_sub = trj.atom_slice(trj.top.select('resi 0 to 100'))
    frames = [get_residue_distance_for_frame(trj_sub, i) for i in range(1,3)]
    print(frames)
    display_arrays_as_video(frames)

    # "Correct" residue distances determined using PyMOL, a standard interface
    # for visualizing 3D molecules (distances limited to 3 decimal places)
    assert (round(calculate_residue_distance(trj[0], 426, 427), 3) == 7.525)
    assert (round(calculate_residue_distance(trj[0], 3, 430), 3) == 22.043)

    # Multi-frame exception
    assert (round(calculate_residue_distance(trj[0:10], 3, 430), 3) == 22.043)

