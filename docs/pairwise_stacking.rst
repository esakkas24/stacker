.. _pairwise_stacking:

********************************************
Making Pairwise Stacking Fingerprints (PSFs)
********************************************

.. currentmodule:: stacker.pairwise_distance

Identifying a Pi-Stacking Pair
------------------------------
Users can identify a residue pair engaging in pi-stacking
by using the :ref:`System Stacking Fingerprint Pipeline <system_stacking>`.

Filter Molecular Dynamics Trajectory
------------------------------------
In this example Pipeline, we use a trajectory of the 
`Ribosome CAR-mRNA Interaction Surface <https://www.mdpi.com/1422-0067/23/3/1417>`_
found in the `StACKER GitHub Repository <https://github.com/esakkas24/stacker/tree/main/testing>`_.
We filter this trajectory to the pi-stacking residue pair.

MD Files are provided for testing convenience in the testing folder:

- ``first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd``: A 10-frame trajectory file with all atoms/residues.
- ``5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop``: The associated Topology File with the above trajectory.

.. currentmodule:: stacker.file_manipulation

Using the :ref:`System Stacking Fingerprint Pipeline <system_stacking>` we found that residues 426 and 427
are pi-stacking. We filter the data from all residues to 426 and 427 and output it to a trajectory ``.pdb`` file
using :func:`filter_traj_to_pdb`::

    >>> import stacker as st
    >>> st.filter_traj_to_pdb("first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd", 
    ...                       topology_filename = "5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop", 
    ...                       output_pdb_filename = "first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.pdb",
    ...                       residues_desired = {426, 427}, 
    ...                       atomnames_desired = {"C2", "C4", "C6"})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    WARNING: Output file atom, residue, and chain indices are zero-indexed
    Filtered trajectory written to:  first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.pdb

So the file ``first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.pdb`` contains the trajectory information for only
residues 426 and 427 for the first 10 frames. 

Write PSF Data to CSV
---------------------

Instead of the 10-frame file we've provided above, that same 
`GitHub testing <https://github.com/esakkas24/stacker/tree/main/testing>`_
repository includes:

- 5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb : A trajectory PDB of residues 426+427 with 3200 frames.

.. currentmodule:: stacker.residue_movement

The data for a Pairwise Stacking Fingerprint (PSF) is outlined in `Bottaro et al <>`_.

We use :func:`write_bottaro_to_csv` to extract this data from the PDB to a CSV.::
    >>> write_bottaro_to_csv


Comparing Pi-Stacking Pairs
----------------------------