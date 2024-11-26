"""
Create and Analyze System Stacking Fingerprints (SSFs)

This module contains the functions to create and analyze System
Stacking Fingerprints. It allows the user to generate pairwise
distance data across a trajectory and compare these matrices across
different trajectories.
"""

import mdtraj as md
import numpy as np
from numpy import typing
from .residue_movement import calc_center_3pts
from .vector import *
from .file_manipulation import SmartIndexingAction
from .visualization import NoResidues, create_axis_labels, display_arrays_as_video
import sys
import concurrent.futures
import functools
import math

_NUCLEOTIDE_NAMES = {"A", "A5", "A3", "G", "G5", "G3", "C", "C5", "C3",
                     "T" "T5", "T3", "U", "U5", "U3", "INO"}

def calculate_residue_distance(trj: md.Trajectory, 
                               res1: int, 
                               res2: int, 
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
    trj : md.Trajectory
        Single frame trajectory.
    res1 : int
        1-indexed residue number of the first residue (PDB Column 5).
    res2 : int
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

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                              'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                              residues = {426,427}, 
    ...                              atoms = {'C2','C4','C6'})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    >>> distance_vec = st.calculate_residue_distance(
    ...     trajectory=filtered_traj, 
    ...     res1_num=426, 
    ...     res2_num=427, 
    ...     res1_atoms=("C2", "C4", "C6"), 
    ...     res2_atoms=("C2", "C4", "C6"), 
    ...     frame=1
    ... )
    >>> distance_vec.magnitude()
    7.5253396
    """
    trj = trj[frame-1]

    # Correct for mdtraj 0-indexing
    res1 = res1 - 1 
    res2 = res2 - 1

    topology = trj.topology
    res1_atom_indices = topology.select("resSeq " + str(res1))
    res2_atom_indices = topology.select("resSeq " + str(res2))
    res1_name = topology.atom(res1_atom_indices[0]).residue.name
    res2_name = topology.atom(res2_atom_indices[0]).residue.name

    if (res1_name not in _NUCLEOTIDE_NAMES) or (res2_name not in _NUCLEOTIDE_NAMES):
        return Vector(0,0,0)
    
    desired_res1_atom_indices = topology.select("(name " + res1_atoms[0] + " or name " + res1_atoms[1] + " or name " + res1_atoms[2] + ") and residue " + str(res1))
    desired_res2_atom_indices = topology.select("(name " + res2_atoms[0] + " or name " + res2_atoms[1] + " or name " + res2_atoms[2] + ") and residue " + str(res2))

    # convert nanometer units in trajectory.xyz to Angstroms
    res1_atom_xyz = trj.xyz[0, desired_res1_atom_indices, :] * 10
    res2_atom_xyz = trj.xyz[0, desired_res2_atom_indices, :] * 10
    vectorized_res1_atom_xyz = [Vector(x,y,z) for [x,y,z] in res1_atom_xyz]
    vectorized_res2_atom_xyz = [Vector(x,y,z) for [x,y,z] in res2_atom_xyz]
    res1_center_of_geometry = calc_center_3pts(*vectorized_res1_atom_xyz)
    res2_center_of_geometry = calc_center_3pts(*vectorized_res2_atom_xyz)

    distance_res12 = res2_center_of_geometry - res1_center_of_geometry
    return distance_res12

def get_residue_distance_for_frame(trj: md.Trajectory, 
                                   frame: int, 
                                   res1_atoms: tuple = ("C2", "C4", "C6"),
                                   res2_atoms: tuple = ("C2", "C4", "C6"),
                                   write_output: bool = True) -> typing.ArrayLike:
    """
    Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.

    Calculates the distances between all pairs of residues in a given frame
    of a trajectory. Outputs as a square matrix with all residues on each
    side. This is the data behind a System Stacking Fingerprint (SSF)

    Parameters
    ----------
    trj : md.Trajectory
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
    pairwise_distances : np.typing.ArrayLike
        Matrix where position (i, j) represents the distance from residue i to residue j.
    
    See Also
    --------
    get_residue_distance_for_trajectory : Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory
    filter_traj : Filters an input trajectory to only the specified atoms and residues
    mdtraj.load : Load a trajectory+topology file

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             residues = '2-5,13-16,23-31,46-51,65-76,88-104,122-141,164-175,184-198,288-289,401-415,420-430', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssf = st.get_residue_distance_for_frame(filtered_traj, frame = 2, write_output = False)
    >>> ssf.shape
    (127, 127)

    """
    trj = trj[frame-1]
    topology = trj.topology
    n_residues = trj.n_residues
    res_indices = [res.resSeq for res in trj.topology.residues]
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
                    pairwise_distances[mat_i,mat_j] = calculate_residue_distance(trj, i+1, j+1, res1_atoms, res2_atoms)
            mat_j+=1
        mat_i+=1
        sys.stdout.flush()
    print(f"\nFrame {frame} done.")
    get_magnitude = np.vectorize(Vector.magnitude)
    pairwise_res_magnitudes = get_magnitude(pairwise_distances)
    return(pairwise_res_magnitudes)

def get_residue_distance_for_trajectory(trj: md.Trajectory, 
                                        frames : typing.ArrayLike | str | set = {},
                                        res1_atoms: tuple = ("C2", "C4", "C6"),
                                        res2_atoms: tuple = ("C2", "C4", "C6"),
                                        threads: int = 1,
                                        write_output: bool = True) -> typing.ArrayLike:
    """
    Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory.

    Calculates the distances between all pairs of residues for all frames
    of a trajectory. Outputs as a square matrix with all residues on each
    side. This is the data behind a System Stacking Fingerprint (SSF)

    Parameters
    ----------
    trj : md.Trajectory
        Trajectory to analyze (must have a topology).
        Output of st.filter_traj() or mdtraj.load()
    frames : array_like or str
        list of frame indices to analyze (1-indexed).
        Accepts smart-indexed str representing a list of frames (e.g '1-5,6,39-48')
        If empty, uses all frames.
    res1_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 1.
    res2_atoms : tuple, default=("C2", "C4", "C6")
        Atom names whose positions are averaged to find the center of residue 2.
    threads : int, default=1
        Number of threads to use for parallel processing.
    write_output : bool, default=True
        If True, displays a loading screen to standard output. Do not use if
        threads > 1.

    Returns
    -------
    ssf_per_frame : array_like
        List where `pairwise_distances[f]` is the output of
        `get_residue_distance_for_frame(trajectory, f, res1_atoms, res2_atoms)`.
    
    See Also
    --------
    system_stacking_fingerprints : Alias for this function
    get_residue_distance_for_frame : Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.
    filter_traj : Filters an input trajectory to only the specified atoms and residues
    mdtraj.load : Load a trajectory+topology file
    display_arrays_as_video : Displays this data as an SSF.

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             residues = '2-5,13-16,23-31,46-51,65-76,88-104,122-141,164-175,184-198,288-289,401-415,420-430', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssfs = st.get_residue_distance_for_trajectory(filtered_traj, frames = '1-3', write_output = False)
    >>> ssfs.shape
    (3, 127, 127)
    """
    frames = SmartIndexingAction.parse_smart_index(frames)
    
    if (frames == {}) or (frames == []):
        frames = [i for i in range(1, trj.n_frames + 1)]

    with concurrent.futures.ProcessPoolExecutor(max_workers = threads) as executor:            
        ssf_per_frame = np.array(list(executor.map(get_residue_distance_for_frame, [trj]*len(frames), frames,
                                    [res1_atoms]*len(frames),[res2_atoms]*len(frames), [write_output]*len(frames))))
    return ssf_per_frame

@functools.wraps(get_residue_distance_for_trajectory)
def system_stacking_fingerprints(*args, **kwargs):
    return get_residue_distance_for_trajectory(*args, **kwargs)

system_stacking_fingerprints.__doc__ = f"""
Alias for `get_residue_distance_for_trajectory()`.

{get_residue_distance_for_trajectory.__doc__}
"""

def get_frame_average(ssfs : typing.ArrayLike) -> typing.ArrayLike:
    '''
    Calculates an average System Stacking Fingerprint (SSF) across multiple SSFs

    Used to calculate an average SSF across multiple frames of a trajectory. Can
    average the result of `get_residue_distance_for_trajectory`

    Parameters
    ----------
    ssfs : numpy.typing.ArrayLike
        List or array of 2D NumPy arrays representing a pairwise distance matrix
        of an MD structure. All 2D NumPy arrays must be of the same dimenstions.
        Output of ``get_residue_distance_for_trajectory()``
        
    Returns
    -------
    avg_frame : numpy.typing.ArrayLike
        A single 2D NumPy array representing a pairwise distance matrix where each
        position i,j is the average distance from residue i to j across all matrices
        in frames.

    See Also
    --------
    get_residue_distance_for_trajectory : Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory
    system_stacking_fingerprints : Alias for get_residue_distance_for_trajectory

    Examples
    --------
    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             residues = '2-5,13-16,23-31,46-51,65-76,88-104,122-141,164-175,184-198,288-289,401-415,420-430', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssfs = st.get_residue_distance_for_trajectory(filtered_traj, frames = '1-3', write_output = False)
    >>> avg_ssf = st.get_frame_average(ssfs)
    >>> avg_ssf.shape
    [FILL]
    '''
    avg_frame = np.mean(ssfs, axis = 0)
    return avg_frame 

def get_top_stacking(trj : md.Trajectory, matrix : typing.ArrayLike, csv : str = '',
                     n_events : int = 5, include_adjacent : bool = False) -> None:
    '''
    Returns top stacking residue pairs for a given System Stacking Fingerprint (SSF)

    Given a trajectory and a SSF made from `get_residue_distance_for_frame()` or `get_frame_average()`
    prints the residue pairings with the strongest stacking events (ie. the residue pairings
    with center of geometry distance closest to 3.5Å). 

    Parameters
    ----------    
    trj : md.Trajectory
        trajectory used to get the stacking fingerprint
    matrix : typing.ArrayLike
        Single-frame SSF created by ``get_residue_distance_for_frame()`` or ``get_frame_average()``
    csv : str, default = '',
        output filename of the tab-separated txt file to write data to. If empty, data printed to standard output
    n_events : int, default = 5
        maximum number of stacking events to display, if -1 display all residue pairings
    include_adjacent : bool, default = False
        True if adjacent residues should be included in the printed output

    See Also
    --------
    get_residue_distance_for_frame : Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.
    get_residue_distance_for_trajectory : Calculates System Stacking Fingerprints (SSFs) for all residues across all frames of a trajectory
    system_stacking_fingerprints : Alias for get_residue_distance_for_trajectory

    Examples
    --------
    We can calculate the stacking events for a single frame:

    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssf = st.get_residue_distance_for_frame(filtered_traj, frame = 2, write_output = False)
    >>> ssf.shape
    (252,252)
    >>> st.get_top_stacking(filtered_traj, ssf)
    Row     Column  Value
    197     195     3.50
    420     413     3.51
    94      127     3.51
    93      130     3.53
    117     108     3.38

    Or we can get most residue pairs that had the most stacking across many frames of a trajectory:

    >>> import stacker as st
    >>> filtered_traj = st.filter_traj('stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', 
    ...                             'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', 
    ...                             atoms = {'C2','C4','C6'})
    >>> ssfs = st.get_residue_distance_for_trajectory(filtered_traj, frames = '1-3', write_output = False)
    >>> avg_ssf = st.get_frame_average(ssfs)
    >>> avg_ssf.shape
    (252, 252)
    >>> st.get_top_stacking(filtered_traj, ssf)
    Row     Column  Value
    130     93      3.56
    108     117     3.44
    195     197     3.61
    127     94      3.65
    47      167     3.66

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

    if csv:
        with open(csv, 'w') as csv_file:
            csv_file.write('Res1\tRes2\tAvg_Dist\n')
            for row, col, value in no_mirrored_indices:
                res1 = increment_residue(str(trj.topology.residue(row).resSeq))
                res2 = increment_residue(str(trj.topology.residue(col).resSeq))
                csv_file.write(f"{res1}\t{res2}\t{value:.2f}\n")
    else:
        print('Res1\tRes2\tAvg_Dist')
        for row, col, value in no_mirrored_indices:
            res1 = increment_residue(str(trj.topology.residue(row).resSeq))
            res2 = increment_residue(str(trj.topology.residue(col).resSeq))
            print(f"{res1}\t{res2}\t{value:.2f}")
    
def increment_residue(residue : str) -> str:
    '''
    Increments residue ID by 1
    
    Useful when converting from mdtraj 0-index residue naming to 1-indexed
    
    Parameters
    ----------
    residue : str
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
    letter_part = ''.join(filter(str.isalpha, residue))
    number_part = ''.join(filter(str.isdigit, residue))
    incremented_number = str(int(number_part) + 1)
    return letter_part + incremented_number

def load_ssfs(file : str) -> typing.ArrayLike:
    """
    Loads a list of SSFs created by ``stacker -s ssf -d OUTFILE``

    Loads a list of SSFs where each element is an SSF from a given
    frame. This is ``OUTFILE`` when running ``stacker -s ssf -d OUTFILE``
    and quickly provides saved SSF data rather than recalculating SSFs
    for a given trajectory.

    Parameters
    ----------
    file : str 
        outfile from ``stacker -s ssf -d OUTFILE``. 
        Must be ``.txt`` or ``.txt.gz``.

    Returns
    -------
    ssfs : numpy.typing.ArrayLike
        List or array of 2D NumPy arrays representing a pairwise distance matrix
        of an MD structure. All 2D NumPy arrays must be of the same dimenstions.
        Output of ``get_residue_distance_for_trajectory()``

    See Also
    --------
    system_stacking_fingerprints : calculates ``ssfs`` rather than loading it like this function.
    get_residue_distance_for_trajectory : Alias for :func:`system_stacking_fingerprints()`
    get_frame_average : calculates average SSF for a trajectory. ``matrix`` parameter is the output of this
    """
    flattened_ssf = np.loadtxt(file)
    ssfs = flattened_ssf.reshape(flattened_ssf.shape[0], math.isqrt(flattened_ssf.shape[1]), math.isqrt(flattened_ssf.shape[1]))
    return ssfs

class MultiFrameTraj(Exception):
    """
    A multi-frame trajectory is passed to a one-frame function

    Raised if a multi-frame trajectory is passed to a function that
    only works on one trajectory (eg. calculate_residue_distance_vector())
    """
    pass

if __name__ == "__main__":
    trajectory_file = '../testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    topology_file = '../testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
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