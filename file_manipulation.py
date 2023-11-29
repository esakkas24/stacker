import mdtraj as md

def filter_traj_to_pdb(trajectory_filename : str, topology_filename : str, 
                        atomnames_desired : set = {},
                        residues_desired : set = {}) -> md.Trajectory:
    '''Filters an input trajectory to only the specified atoms and residues

    Filteres an input trajectory that contains all of the atoms in a trajectory to only
        the desired atoms at the desired residues (eg. the atoms necessary to find the 
        center of geometry of a residue) 

    Args:
        trajectory_filename (str) : path to file of the concatenated trajectory
            - Should be concatenation of the 1 in 50 frames sampled trajectories for each replicate
        topology_filename (str) : path to file of the topology of the molecule
        atomnames_desired (set) : atomnames to keep in the trajectory
        residues_desired (set) : residue numbers of residues to keep in the trajectory

    Returns:
        filtered_trajectory (md.Trajectory) : a trajectory object representing the filtered structure across all frames
    '''
    print("Reading trajectory...")
    trajectory = md.load(trajectory_filename, top = topology_filename)
    
    print("Getting topology...")
    topology = trajectory.topology
    
    print("Filtering trajectory...")
    atomnames_query = " or ".join([f"name == '{atom}'" for atom in atomnames_desired])
    residues_query = " or ".join([f"residue == {resnum}" for resnum in residues_desired])
    atom_indices_selection = topology.select('(' + atomnames_query + ') and (' + residues_query + ')')
    filtered_trajectory = trajectory.atom_slice(atom_indices_selection)
    return filtered_trajectory

def file_convert(input_file : str, ouput_filetype : str) -> None:
    '''Converts trajectory input file to new output type
    
    Args:
        input file (str) : path name of the input file to be converted
        output_file (str) : desired filetype to convert to (mdcrd, pdb, etc.)
    Returns:
        None

    Uses mdtraj.save() commands to save to convert trajectories to various filetypes
        mdtraj.save_mdcrd(), mdtraj.save_pdb(), mdtraj.save_xyz(), etc'''
    pass

if __name__ == "__main__":
    filtered_traj = filter_traj_to_pdb('first10_5JUP_N2_wCCC_+1GCU_nowat.mdcrd', '5JUP_N2_wCCC_+1GCU_nowat.prmtop', {'C2','C4','C6'}, {426,427})
    print([residue for residue in filtered_traj.topology.residues])
    print([atom for atom in filtered_traj.topology.atoms])
    table, bonds = filtered_traj.topology.to_dataframe()
    print(table)