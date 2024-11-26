.. _command_line_options:

***********************
Command Line Options
***********************

After installing StACKER, the command ``stacker`` is added
to terminal. Different routines can be run in the command Line
using ``stacker -s ROUTINE``. Running ``stacker --help`` 
will give the subroutine options::

    [user]$ stacker --help
    usage: stacker -s ROUTINE [-h]

    Wrapper to run stacker subroutines using the -s flag.
    More info on each routine given by `stacker -s ROUTINE -h`

    options:
    -s, --script ROUTINE  Name of command to use. Options for ROUTINE:
                            
                            filter_traj:
                                    filters trajectory and topology files to desired residue numbers and atom names
                            bottaro OR pairwise OR psf:
                                    Create Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)
                            res_distance:
                                    Get the distance between two residues in a given frame
                            system OR ssf:
                                    Create a System Stacking Fingerprint of distances by residue
                            stack_events:
                                    Get list of residues with most stacking events (distance closest to 3.5Å)
                            compare:
                                    Get the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events
    -h, --help            show this help message and exit

Filtering Trajectory
--------------------
Running ``stacker -s filter_traj`` will filter a trajectory by
residue/atom indices and output the trajectory to a ``.pdb``.
Running ``stacker -s filter_traj --help`` will give more flag options::

    [user]$ stacker -s filter_traj --help
    usage: stacker -s ROUTINE -trj TRAJECTORY_FILENAME -top TOPOLOGY_FILENAME -o OUTPUT_FILE [-r RESIDUES] [-a ATOM_NAMES] [-h]

    Filters trajectory and topology files to desired residue numbers and atom names and outputs to a PDB

    Examples:
    [user]$ stacker -s filter_traj -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 426,427 -a C2,C4,C6

    options:
    -s, --script ROUTINE  Name of command to use. Options for ROUTINE:
                            
                            filter_traj:
                                    filters trajectory and topology files to desired residue numbers and atom names
                            bottaro OR pairwise OR psf:
                                    Create Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)
                            res_distance:
                                    Get the distance between two residues in a given frame
                            system OR ssf:
                                    Create a System Stacking Fingerprint of distances by residue
                            stack_events:
                                    Get list of residues with most stacking events (distance closest to 3.5Å)
                            compare:
                                    Get the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events
    -r, --residues RESIDUES
                            Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)
    -a, --atom_names ATOM_NAMES
                            Comma-separated list of atom names to filter
    -h, --help            show this help message and exit

    Required Arguments:
    -trj, --trajectory TRAJECTORY_FILENAME
                            Filepath to trajectory file for the MD simulation
    -top, --topology TOPOLOGY_FILENAME
                            Filepath to Topology file for the MD simulation
    -o, --output OUTPUT_FILE
                            Filepath of PDB to output to

.. currentmodule:: stacker.file_manipulation

This is the command line equivalent of :func:`filter_traj_to_pdb`. So, the following are equivalent::
        
        >>> st.filter_traj_to_pdb(
        ...     trajectory_filename = "../testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd",
        ...     topology_filename = "../testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop",
        ...     output_pdb_filename = "../testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb",
        ...     residues_desired = "426,427",
        ...     atomnames_desired = {"C2","C4","C6"}
        ... )

Get Distance Between Residues
-----------------------------
Running ``stacker -s res_distance`` will give the distance in Angstroms
between the center of geometry of two residues for a given frame.
Running ``stacker -s res_distance --help`` will give more flag options::

    [user]$ stacker -s res_distance --help
    usage: stacker -s ROUTINE -trj TRAJECTORY_FILENAME -top TOPOLOGY_FILENAME [-f FRAME_NUM] -r RESIDUES [-b N_FRAMES] [-a ATOM_NAMES] [-h]

    Get the distance between two residues in a given frame

    Examples:
    [user]$ stacker -s res_distance -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -f 2 --residues 426,427 --atom_names C2,C4,C6

    options:
    -s, --script ROUTINE  Name of command to use. Options for ROUTINE:
                            
                            filter_traj:
                                    filters trajectory and topology files to desired residue numbers and atom names
                            bottaro OR pairwise OR psf:
                                    Create Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)
                            res_distance:
                                    Get the distance between two residues in a given frame
                            system OR ssf:
                                    Create a System Stacking Fingerprint of distances by residue
                            stack_events:
                                    Get list of residues with most stacking events (distance closest to 3.5Å)
                            compare:
                                    Get the most changed stacking events between two fingerprints using the outputs of stacker -s stack_events
    -f, --frame FRAME_NUM
                            1-indexed Frame Number within trajectory to analyze
    -b, --bootstrap N_FRAMES
                            Run bootstrap analysis on this residue pairing, sampling N_FRAMES with replacement
    -a, --atom_names ATOM_NAMES
                            Comma-separated list of atom names. Three required to get center of geometry for a residue. default = C2,C4,C6
    -h, --help            show this help message and exit

    Required Arguments:
    -trj, --trajectory TRAJECTORY_FILENAME
                            Filepath to trajectory file for the MD simulation
    -top, --topology TOPOLOGY_FILENAME
                            Filepath to Topology file for the MD simulation
    -r, --residues RESIDUES
                            Smart-indexed list of 1-indexed residues, must provide only 2 residues, accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)

