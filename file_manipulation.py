import mdtraj as md

def filter_traj(trajectory_filename : str, topology_filename : str, 
                        residues_desired : set = {}, atomnames_desired : set = {}) -> md.Trajectory:
    '''Filters an input trajectory to only the specified atoms and residues

    Filteres an input trajectory that contains all of the atoms in a trajectory to only
        the desired atoms at the desired residues (eg. the atoms necessary to find the 
        center of geometry of a residue). If residues_desired or atomnames_desired are
        empty, all residues or atoms are included respectively.

    Args:
        trajectory_filename : str
            path to file of the concatenated trajectory. Should be resampled to the
            1 in 50 frames sampled trajectories for each replicate.
        topology_filename : str
            path to file of the topology of the molecule
        residues_desired : set
            1-indexed residue numbers of residues to keep in the trajectory
        atomnames_desired : set 
            atomnames to keep in the trajectory

    Returns:
        filtered_trajectory : md.Trajectory
            a trajectory object representing the filtered structure across all frames
    '''
    print("WARNING: Residue Indices are expected to be 1-indexed")
    
    print("Reading trajectory...")
    trajectory = md.load(trajectory_filename, top = topology_filename)
    
    print("Reading topology...")
    topology = trajectory.topology
    
    print("Filtering trajectory...")
    # make resSeq 0-indexed for mdtraj query
    residues_desired = {resnum-1 for resnum in residues_desired} 

    atomnames_query = " or ".join([f"name == '{atom}'" for atom in atomnames_desired])
    residues_query = " or ".join([f"residue == {resnum}" for resnum in residues_desired])

    if len(atomnames_query) == 0:
        if len(residues_query) == 0:
            filtered_trajectory = trajectory
        else:
            atom_indices_selection = topology.select(residues_query)
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
    else:
        if len(residues_query) == 0:
            atom_indices_selection = topology.select(atomnames_query)
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
        else:
            atom_indices_selection = topology.select('(' + atomnames_query + ') and (' + residues_query + ')')
            filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
    print("WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed")

    return filtered_trajectory

def filter_traj_to_pdb(trajectory_filename : str, topology_filename : str, 
                       output_pdb_filename : str,
                        residues_desired : set = {},
                        atomnames_desired : set = {}) -> None:
    '''Filters an input trajectory to only the specified atoms and residues and outputs to pdb

    Filteres an input trajectory that contains all of the atoms in a trajectory to only
        the desired atoms at the desired residues (eg. the atoms necessary to find the 
        center of geometry of a residue) and writes the output to a specified pdb file.
        If residues_desired or atomnames_desired are empty, all residues or atoms are included respectively.

    Args:
        trajectory_filename : str
            path to file of the concatenated trajectory. Should be resampled to the
            1 in 50 frames sampled trajectories for each replicate.
        topology_filename : str
            path to file of the topology of the molecule
        output_pdb_filename : str
            path to the output pdb file
        residues_desired : set
            1-indexed residue numbers of residues to keep in the trajectory
        atomnames_desired : set 
            atomnames to keep in the trajectory

    Returns:
        None
    '''
    print("WARNING: Residue Indices are expected to be 1-indexed")
    filtered_trajectory = filter_traj(trajectory_filename, topology_filename, residues_desired, atomnames_desired)
    filtered_trajectory.save_pdb(output_pdb_filename)
    print("WARNING: Output file atom, residue, and chain indices are zero-indexed")
    print("Filtered trajectory written to: ", output_pdb_filename)


def file_convert(trajectory_filename : str, topology_filename : str, output_file : str) -> None:
    '''Converts mdcrd trajectory input file to new output type
    
    Output filetype determined by output_file extension. Uses mdtraj.save() commands to convert 
    trajectory files to various filetypes mdtraj.save_mdcrd(), mdtraj.save_pdb(), mdtraj.save_xyz(), etc

    Args:
        trajectory_filename : str
            path to file of the concatenated trajectory (.mdcrd file). Should be resampled to the
            1 in 50 frames sampled trajectories for each replicate.
        topology_filename : str
            path to file of the topology of the molecule (.prmtop file)
        output_file : str
            output filename (include .mdcrd, .pdb, etc.)
    Returns:
        None
    '''
    print("WARNING: Output file atom, residue, and chain indices are zero-indexed")
    trajectory = md.load(trajectory_filename, top = topology_filename)
    trajectory.save(output_file)
    print("Trajectory written to: ", output_file)

if __name__ == "__main__":
    # filter_traj tests
    print('Known Res: 426 = G and 427 = C')
    filtered_traj = filter_traj('first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', '5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', {426,427}, {'C2','C4','C6'})
    table, bonds = filtered_traj.topology.to_dataframe()
    print(table)

    ### No Filtering
    print("No Filtering, known trj has 12089 atoms")
    filtered_traj = filter_traj('first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd', '5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop', residues_desired={}, atomnames_desired={})
    table, bonds = filtered_traj.topology.to_dataframe()
    print(table)