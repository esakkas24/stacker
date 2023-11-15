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
        filtered_pdb_filename (str) : path to file of the filtered PDB to write to
        atomnames_desired (set) : atomnames to keep if the line represents this type of atom
        residues_desired (set) : residue numbers of residues to keep if the line represents this type of residue

    Returns:
        filtered_trajectory (md.Trajectory) : a trajectory object representing the filtered structure across all frames
    '''
    print("Reading trajectory...")
    trajectory = md.load(trajectory_filename, top = topology_filename)
    first_10_frames = trajectory[0:10]
    first_10_frames.save_mdcrd('first10_'+trajectory_filename)
    print("Getting topology...")
    topology = trajectory.topology
    
    print("Filtering trajectory...")
    atomnames_query = " or ".join([f"name == '{atom}'" for atom in atomnames_desired])
    residues_query = " or ".join([f"residue == {resnum}" for resnum in residues_desired])
    selection_indices = topology.select('(' + atomnames_query + ') and (' + residues_query + ')')
    filtered_trajectory = trajectory.atom_slice(selection_indices)
    filtered_trajectory.save_pdb('test.pdb')
    return filtered_trajectory
