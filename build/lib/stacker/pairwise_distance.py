import mdtraj as md
import numpy as np
from numpy import typing
from .residue_movement import calc_center_3pts
from .vector import *
from .visualization import NoResidues, create_axis_labels, display_arrays_as_video
import sys
import concurrent.futures

class MultiFrameTraj(Exception):
    pass

_NUCLEOTIDE_NAMES = {"A", "A5", "A3", "G", "G5", "G3", "C", "C5", "C3",
                     "T" "T5", "T3", "U", "U5", "U3", "INO"}

def calculate_residue_distance(trajectory: md.Trajectory, 
                               res1_num: int, 
                               res2_num: int, 
                               res1_atoms: tuple = ("C2", "C4", "C6"),
                               res2_atoms: tuple = ("C2", "C4", "C6"),
                               frame: int = 1) -> Vector:
    """
    Calculates the vector between two residues with x, y, z units in Angstroms.

    Calculates the distance between the center of two residues. The center is defined
    by the average x, y, z position of three passed atoms for each residue (typically
    every other carbon on the 6-carbon ring of the nucleotide base).

    Parameters
    ----------
    trajectory : md.Trajectory
        Single frame trajectory.
    res1_num : int
        1-indexed residue number of the first residue (PDB Column 5).
    res2_num : int
        1-indexed residue number of the second residue (PDB Column 5).
    res1_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 1.
    res2_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 2.
    frame : int, default=1
        1-indexed frame number of trajectory to calculate the distance.

    Returns
    -------
    distance_res12 : Vector
        Vector from the center of geometry of residue 1 to residue 2.
    
    See Also
    --------
    get_residue_distance_for_frame : Calculates pairwise distances between all residues in a given frame.
    """
    trajectory = trajectory[frame-1]

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

def get_residue_distance_for_frame(trajectory: md.Trajectory, 
                                   frame: int, 
                                   res1_atoms: tuple = ("C2", "C4", "C6"),
                                   res2_atoms: tuple = ("C2", "C4", "C6"),
                                   write_output: bool = True) -> typing.ArrayLike:
    """
    Calculates pairwise distances between all residues in a given frame.

    Parameters
    ----------
    trajectory : md.Trajectory
        Trajectory to analyze (must have a topology).
    frame : int
        1-indexed frame to analyze.
    res1_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 1.
    res2_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 2.
    write_output : bool, default=True
        If True, displays a loading screen to standard output.

    Returns
    -------
    pairwise_distances : array_like
        Matrix where position (i, j) represents the distance from residue i to residue j.
    
    """
    trajectory = trajectory[frame-1]
    topology = trajectory.topology
    n_residues = trajectory.n_residues
    res_indices = [res.resSeq for res in trajectory.topology.residues]
    zero_vector = Vector(0,0,0)

    pairwise_distances = np.full((n_residues, n_residues), zero_vector)

    mat_i = 0
    for i in res_indices:
        if write_output:
            percent_done = round((mat_i+1) / n_residues * 100, 2)
            sys.stdout.write(f'\rLoading: [{"#" * int(percent_done)}{" " * (100 - int(percent_done))}] Current Residue: {mat_i+1}/{n_residues} ({percent_done}%)')
        mat_j = 0
        res1_name = topology.residue(mat_i).name
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
        sys.stdout.flush()
    print(f"\nFrame {frame} done.")
    get_magnitude = np.vectorize(Vector.magnitude)
    pairwise_res_magnitudes = get_magnitude(pairwise_distances)
    return(pairwise_res_magnitudes)

def get_residue_distance_for_trajectory(trajectory: md.Trajectory, 
                                        frames: typing.ArrayLike,
                                        res1_atoms: tuple = ("C2", "C4", "C6"),
                                        res2_atoms: tuple = ("C2", "C4", "C6"),
                                        threads: int = 1) -> typing.ArrayLike:
    """
    Calculates pairwise distances for all residues across all frames of a trajectory.

    Parameters
    ----------
    trajectory : md.Trajectory
        Trajectory to analyze (must have a topology).
    frames : array_like
        Frame indices to analyze (1-indexed).
    res1_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 1.
    res2_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 2.
    threads : int, default=1
        Number of threads to use for parallel processing.

    Returns
    -------
    ssf_per_frame : array_like
        List where `pairwise_distances[f]` is the output of
        `get_residue_distance_for_frame(trajectory, f, res1_atoms, res2_atoms)`.
    
    """
    write_output = False
    if threads <= 1:
        write_output = True

    with concurrent.futures.ProcessPoolExecutor(max_workers = threads) as executor:            
        SSFs = np.array(list(executor.map(get_residue_distance_for_frame, [trajectory]*len(frames), frames,
                                    [res1_atoms]*len(frames),[res2_atoms]*len(frames), [write_output]*len(frames))))
    return SSFs


def increment_residue(residue_id : str) -> str:
    '''
    Increments residue ID by 1
    
    Useful when converting from mdtraj 0-index residue naming to 1-indexed
    
    Parameters
    ----------
    residue_id : str
        The residue id given by trajectory.topology.residue(i)

    Returns
    -------
        incremented_id : str
            The residue id with the sequence number increased by 1

    Examples
    --------
    >>> increment_residue('G43')
    'G44'

    '''
    letter_part = ''.join(filter(str.isalpha, residue_id))
    number_part = ''.join(filter(str.isdigit, residue_id))
    incremented_number = str(int(number_part) + 1)
    return letter_part + incremented_number

def get_top_stacking(trajectory : md.Trajectory, matrix : typing.ArrayLike, output_csv : str = '',
                     n_events : int = 5, include_adjacent : bool = False) -> None:
    '''
    Returns top stacking events for a given stacking fingerprint

    Given a trajectory and a stacking fingerprint made from get_residue_distance_for_frame(),
    prints the residue pairings with the strongest stacking events (ie. the residue pairings
    with center of geometry distance closest to 3.5Å)

    Parameters
    ----------    
    trajectory : md.Trajectory
        trajectory used to get the stacking fingerprint
    matrix : typing.ArrayLike
        stacking fingerprint matrix created by get_residue_distance_for_frame()
    output_csv : str, default = '',
        output filename of the tab-separated txt file to write data to. If empty, data printed to standard output
    n_events : int, default = 5
        maximum number of stacking events to display, if -1 display all residue pairings
    include_adjacent : bool, default = False
        True if adjacent residues should be included in the printed output

    '''
    top_stacking_indices = np.argsort(np.abs(matrix - 3.5), axis = None)
    rows, cols = np.unravel_index(top_stacking_indices, matrix.shape)
    closest_values = matrix[rows, cols]

    if include_adjacent:
        # non_adjacent_indices includes adjacent indices in this case
        non_adjacent_indices = [(row, col, value) for row, col, value in zip(rows, cols, closest_values) if abs(row - col) > 0]
    else:
        non_adjacent_indices = [(row, col, value) for row, col, value in zip(rows, cols, closest_values) if abs(row - col) > 1]

    no_mirrored_indices = [] # keep only one side of x=y line, since mat[i,j] = mat[j,i]
    for row, col, value in non_adjacent_indices:
        if (col, row, value) not in no_mirrored_indices:
            no_mirrored_indices += [(row, col, value)]
    if n_events == -1: n_events = len(no_mirrored_indices) 
    no_mirrored_indices = no_mirrored_indices[:n_events]

    if output_csv:
        with open(output_csv, 'w') as csv_file:
            csv_file.write('Row\tColumn\tValue\n')
            for row, col, value in no_mirrored_indices:
                res1 = increment_residue(str(trajectory.topology.residue(row).resSeq))
                res2 = increment_residue(str(trajectory.topology.residue(col).resSeq))
                csv_file.write(f"{res1}\t{res2}\t{value:.2f}\n")
    else:
        print('\nRow\tColumn\tValue')
        for row, col, value in no_mirrored_indices:
            res1 = increment_residue(str(trajectory.topology.residue(row).resSeq))
            res2 = increment_residue(str(trajectory.topology.residue(col).resSeq))
            print(f"{res1}\t{res2}\t{value:.2f}")
    
def get_frame_average(frames : typing.ArrayLike) -> typing.ArrayLike:
    '''
    Calculates an average pairwise matrix across multiple frames of a trajectory

    Parameters
    ----------
    frames : numpy.typing.ArrayLike
        List or array of 2D NumPy arrays representing a pairwise distance matrix
        of an MD structure. All 2D NumPy arrays must be of the same dimenstions.
        
    Returns
    -------
    avg_frame : numpy.typing.ArrayLike
        A single 2D NumPy array representing a pairwise distance matrix where each
        position i,j is the average distance from residue i to j across all matrices
        in frames.
    '''
    avg_frame = np.mean(frames, axis = 0)
    return avg_frame

if __name__ == "__main__":
    trajectory_file = 'testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    topology_file = 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    # Load test trajectory and topology
    trj = md.load(trajectory_file, top = topology_file)

    # "Correct" residue distances determined using PyMOL, a standard interface
    # for visualizing 3D molecules (distances limited to 3 decimal places)

    # calculate_residue_distance() tests
    tolerance = 1e-6
    assert round(calculate_residue_distance(trj[0], 426, 427).magnitude(), 3) - 7.525 < tolerance
    assert (round(calculate_residue_distance(trj[0], 3, 430).magnitude(), 3) - 22.043 < tolerance)
    ### Multi-frame exception
    try:
        round(calculate_residue_distance(trj[0:10], 3, 430).magnitude(), 3) - 22.043 < tolerance
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
    residue_selection_query = 'resi 90 to 215'
    frames_to_include = [1,2,3,4,5]

    trj_sub = trj.atom_slice(trj.top.select(residue_selection_query))
    resSeqs = [res.resSeq for res in trj_sub.topology.residues]
    frames = get_residue_distance_for_trajectory(trj_sub, frames_to_include, threads = 5)
    get_top_stacking(trj_sub, frames[0])
    display_arrays_as_video([get_frame_average(frames)], resSeqs, seconds_per_frame=10)

    display_arrays_as_video(frames, resSeqs, seconds_per_frame=10)

    # All Residues one large matrix
    resSeqs = [res.resSeq for res in trj.topology.residues]
    print('\n')
    frames = [get_residue_distance_for_frame(trj, i) for i in range(1,2)]
    display_arrays_as_video(frames, resSeqs, seconds_per_frame=10, tick_distance=20)