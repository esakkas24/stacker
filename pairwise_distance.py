import mdtraj as md
import numpy as np
from numpy import typing
from residue_movement import calc_center_3pts
from vector import *
from visualization import NoResidues, create_axis_labels, display_arrays_as_video

class MultiFrameTraj(Exception):
    pass

_NUCLEOTIDE_NAMES = {"A", "G", "C", "T", "U"}

def calculate_residue_distance(trajectory : md.Trajectory, 
                               res1_num : int, res2_num : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> Vector:
    '''Calculates the vector between two residues with x,y,z units in Angstroms

    Calcualtes the distance between the center of two residues. The center is denoted
        by the average x,y,z position of three passed atoms for each residue (typically
        every other carbon on the 6-C ring of the nucleotide base).

    Args:
        trajectory : md.Trajectory 
            single frame trajectory
        res1_num : int
            the residue number of the first residue (PDB Column 5)
        res2_num : int
            the residue number of the second residue (PDB Column 5)
        res1_atoms : tuple, default = ("C2","C4","C6")
            a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 
        res2_atoms : tuple, default = ("C2","C4","C6")
            a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 [("C2","C4","C6")]
    
    Returns:
        distance_res12 : Vector
            Vector from center of geometry of residue 1 to center of geometry of residue 2
    '''
    if len(trajectory) > 1: raise MultiFrameTraj("calculate_residue_distance() expects a 1-frame trajectory")

    # Correct for mdtraj 0-indexing
    res1_num = res1_num - 1 
    res2_num = res2_num - 1

    topology = trajectory.topology
    res1_atom_indices = topology.select("resSeq " + str(res1_num))
    res2_atom_indices = topology.select("resSeq " + str(res2_num))
    res1_name = topology.atom(res1_atom_indices[0]).residue.name
    res2_name = topology.atom(res2_atom_indices[0]).residue.name

    if (res1_name not in _NUCLEOTIDE_NAMES) or (res2_name not in _NUCLEOTIDE_NAMES):
        return Vector(0,0,0)
    
    desired_res1_atom_indices = topology.select("(name " + res1_atoms[0] + " or name " + res1_atoms[1] + " or name " + res1_atoms[2] + ") and residue " + str(res1_num))
    desired_res2_atom_indices = topology.select("(name " + res2_atoms[0] + " or name " + res2_atoms[1] + " or name " + res2_atoms[2] + ") and residue " + str(res2_num))

    # convert nanometer units in trajectory.xyz to Angstroms
    res1_atom_xyz = trajectory.xyz[0, desired_res1_atom_indices, :] * 10
    res2_atom_xyz = trajectory.xyz[0, desired_res2_atom_indices, :] * 10
    vectorized_res1_atom_xyz = [Vector(x,y,z) for [x,y,z] in res1_atom_xyz]
    vectorized_res2_atom_xyz = [Vector(x,y,z) for [x,y,z] in res2_atom_xyz]
    res1_center_of_geometry = calc_center_3pts(*vectorized_res1_atom_xyz)
    res2_center_of_geometry = calc_center_3pts(*vectorized_res2_atom_xyz)

    distance_res12 = res2_center_of_geometry - res1_center_of_geometry
    return distance_res12

def get_residue_distance_for_frame(trajectory : md.Trajectory, frame : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> typing.ArrayLike:
    '''Calculates pairwise the distance between all residues in a given frame

    Args:
        trajectory : md.Trajectory
            trajectory to analyze (must have topology aspect)
        frame : int
            1-indexed frame to analyze
        res1_atoms : tuple, default = ("C2","C4","C6")
            a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 
        res2_atoms : tuple, default = ("C2","C4","C6")
            a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 
    
    Returns:
        pairwise_distances : array_like
            matrix where position i,j represents the distance from 
            residue i to residue j
    '''
    trajectory = trajectory[frame-1]
    n_residues = trajectory.n_residues
    res_indices = [res.resSeq for res in trajectory.topology.residues]
    zero_vector = Vector(0,0,0)

    pairwise_distances = np.full((n_residues, n_residues), zero_vector)

    mat_i = 0
    for i in res_indices:
        print(i)
        mat_j = 0
        res1_atom_indices = trajectory.topology.select("resSeq " + str(i))
        res1_name = trajectory.topology.atom(res1_atom_indices[0]).residue.name
        for j in res_indices:
            if i == j: 
                pairwise_distances[mat_i,mat_j] = zero_vector
            elif pairwise_distances[mat_j,mat_i] != zero_vector:
                pairwise_distances[mat_i,mat_j] = pairwise_distances[mat_j,mat_i]
            elif any(np.logical_and(pairwise_distances[:mat_i, mat_i] != zero_vector,
                                       pairwise_distances[:mat_i, mat_j] != zero_vector)):
                for intermediate_res in range(0, mat_i):
                    if (pairwise_distances[intermediate_res, mat_i] != zero_vector and pairwise_distances[intermediate_res, mat_j] != zero_vector):
                        pairwise_distances[mat_i,mat_j] = pairwise_distances[intermediate_res, mat_i].scale(-1) + pairwise_distances[intermediate_res, mat_j]
                        break
            else:
                if (res1_name not in _NUCLEOTIDE_NAMES):
                    pairwise_distances[mat_i,:] = Vector(0,0,0)
                    break
                else:
                    pairwise_distances[mat_i,mat_j] = calculate_residue_distance(trajectory, i+1, j+1, res1_atoms, res2_atoms)
            mat_j+=1
        mat_i+=1

    get_magnitude = np.vectorize(Vector.magnitude)
    pairwise_res_magnitudes = get_magnitude(pairwise_distances)
    return(pairwise_res_magnitudes)

if __name__ == "__main__":
    # Load test trajectory and topology
    trj = md.load('first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', top = '5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop')

    # "Correct" residue distances determined using PyMOL, a standard interface
    # for visualizing 3D molecules (distances limited to 3 decimal places)

    # calculate_residue_distance() tests
    assert (round(calculate_residue_distance(trj[0], 426, 427).magnitude(), 3) == 7.525)
    assert (round(calculate_residue_distance(trj[0], 3, 430).magnitude(), 3) == 22.043)
    ### Multi-frame exception
    try:
        round(calculate_residue_distance(trj[0:10], 3, 430).magnitude(), 3) == 22.043
    except MultiFrameTraj:
        print("MultiFrameTraj: calculate_residue_distance_vector() fails on multiple-frame trajectory")

    # create_axis_labels() test
    assert(create_axis_labels([0,1,2,3,4,5,6,7,8,9,10,11,12,98,99,100]) == ([0,10,12,13,15], [0,10,12,98,100]))
    assert(create_axis_labels([94,95,96,97,98,99,100,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428]) == ([0,6,7,17,27], [94,100,408,418,428]))
    ### No passed in residues exception
    try:
        assert(create_axis_labels([]) == ([],[]))
    except NoResidues:
        print("NoResidues: create_axis_labels() fails on empty residue list")

    # get_residue_distance_for_frame() test
    trj_three_residues = trj.atom_slice(trj.top.select('resi 407 or resi 425 or resi 426'))
    assert(np.all(np.vectorize(round)(get_residue_distance_for_frame(trj_three_residues, 2), 3) == np.array([[0,      8.231,   11.712], 
                                                                                                              [8.231,  0,       6.885], 
                                                                                                               [11.712, 6.885,   0]])))

    # display_arrays_as_video() tests
    trj_sub = trj.atom_slice(trj.top.select('resi 50 to 100 or resi 150 to 200'))
    resSeqs = [res.resSeq for res in trj_sub.topology.residues]
    frames = [get_residue_distance_for_frame(trj_sub, i) for i in range(1,10)]
    display_arrays_as_video(frames, resSeqs, seconds_per_frame=1)

    resSeqs = [res.resSeq for res in trj.topology.residues]
    frames = [get_residue_distance_for_frame(trj, i) for i in range(1,2)]
    display_arrays_as_video(frames, resSeqs, seconds_per_frame=60, tick_distance=20)