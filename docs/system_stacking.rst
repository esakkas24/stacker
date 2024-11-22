.. _system_stacking:

********************************************
Making System Stacking Fingerprints (SSFs)
********************************************

This pipeline creates System Stacking Fingerprints (SSFs) which can be used 
to analyze system-wide pi-stacking in a molecular dynamics trajectory.

Load Trajectory
---------------

To analyze pi-stacking between two nucleotides, we compare their center of
geometry (COG) distance, the distance between the centers of their 6-membered
(pyrimidine) rings. For this, the only atoms we need are C2, C4, and C6 to triangulate
the COG of each ring. For reference, the Carbon Numbering of each nucleotide:

.. image:: images/nucleotide_numbering.png

In this example pipeline, we use a trajectory of the 
`Ribosome CAR-mRNA Interaction Surface <https://www.mdpi.com/1422-0067/23/3/1417>`_
found in the `StACKER GitHub Repository <https://github.com/esakkas24/stacker/tree/main/testing>`_.
We filter this trajectory to the pi-stacking residue pair.

MD Files are provided for testing convenience in the testing folder:

- ``first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd``: A 10-frame trajectory file with all atoms/residues.
- ``5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop``: The associated Topology File with the above trajectory.

.. currentmodule:: stacker.file_manipulation

We filter to atoms C2, C4, and C6 using :func:`filter_traj`::

    >>> import stacker as st
    >>> filtered_traj = st.filter_traj("first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd",
    ...                 topology_filename = "5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop",
    ...                 atomnames_desired = {"C2", "C4", "C6"})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    >>> filtered_traj
    <mdtraj.Trajectory with 10 frames, 756 atoms, 252 residues, without unitcells at 0x1164eab10>

Now the Python variable ``filtered_traj`` contains 10 frames of C2, C4, C6 information for all nucleotides.

Calculate Distance Between Residues
------------------------------------

.. currentmodule:: stacker.pairwise_distance

Immediately, we can check if a residue pair is pi-stacking in a given frame. We calculate 
the distance between the COG of two residues using :func:`calculate_residue_distance`.
COG distance close to 3.5 Å indicates pi-stacking. For instance, the A-site mRNA codons 
(residues 422, 423, and 424) are likely pi-stacking::

    >>> distance_vec = st.calculate_residue_distance(
    ...     filtered_traj,
    ...     res1_num = 422,
    ...     res2_num = 423,
    ...     frame = 3
    ... )
    >>> distance_vec.magnitude()
    3.509042

Create a System Stacking Fingerprint (SSF)
------------------------------------------

An SSF calculates the COG distance for all nucleotide pairs in a trajectory frame.
The result is a square matrix where position (i, j) represents the distance 
from residue i to residue j. This is done with :func:`get_residue_distance_for_frame`

The ``filtered_traj`` object has 252 residues, so we create a 252 x 252 SSF::

    >>> ssf = st.get_residue_distance_for_frame(filtered_traj, frame = 2, write_output = True)
    Loading: [####################################################################################################] Current Residue: 252/252 (100.0%)
    Frame 2 done.
    >>> ssf.shape
    (252, 252)

We can calculate the SSF for multiple frames of a trajectory using :func:`system_stacking_fingerprints`, which
accepts smart indexing of frames (eg. 1-3,15-17 = 1,2,3,15,16,17)::

    >>> ssfs = st.system_stacking_fingerprints(filtered_traj, frames = '1-3')
    >>> ssfs.shape
    (3, 252, 252)

``ssfs`` is now a list, where ``ssfs[i]`` is the SSF for frame ``i`` (0-indexed frame). If ``frames`` is empty,
the SSF will be calculated for all frames. For multi-frame
trajectories, it is recommended to use the ``threads`` option to parallelize, calculating the SSF for multiple
frames at once. When parralelizing, turn off standard output with ``write_output``::

    >>> ssfs = st.system_stacking_fingerprints(
    ...     filtered_traj,
    ...     frames = '1-10',
    ...     threads = 10,
    ...     write_output = False
    ... )
    >>> ssfs.shape
    (10, 252, 252)

Get the Average SSF for a Trajectory
------------------------------------------

Single-frame SSFs are rarely as illuminating as the average SSF for all frames
of a trajectory. Users can create this using :func:`get_frame_average`, using the 
output from :func:`system_stacking_fingerprints` in the step above::

    >>> ssfs = st.system_stacking_fingerprints(
    ...     filtered_traj,
    ...     frames = '1-10',
    ...     threads = 10,
    ...     write_output = False
    ... )
    >>> avg_ssf = st.get_frame_average(ssfs)
    >>> avg_ssf.shape
    (252, 252)

``avg_ssf`` contains averaged stacking information throughout the trajectory.
The next step shows how to analyze these results.

How to Use an SSF
-----------------

:func:`get_top_stacking` will give the stacking pairs with the most
pi-stacking (ie. closest to 3.5Å).




.. currentmodule:: stacker.pairwise_distance

It is recommended that you parallelize for large trajectories, it calculates the SSF
per frame in parallel.

Comparing Pi-Stacking Pairs
----------------------------