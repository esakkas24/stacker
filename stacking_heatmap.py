import mdtraj as md
import pandas

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
    
filtered_traj = filter_traj_to_pdb('5JUP_N2_wCCC_+1GCU_nowat.mdcrd', '5JUP_N2_wCCC_+1GCU_nowat_posttleap.pdb', {'C2','C4','C6'}, {426,427})
print([residue for residue in filtered_traj.topology.residues])
table, bonds = filtered_traj.topology.to_dataframe()
print(table)



def calculate_residue_distance(res1_num : int, res2_num : int, 
                                res1_atoms : tuple = ("C2","C4","C6"),
                                res2_atoms : tuple = ("C2","C4","C6")) -> float:
    '''Calculates the distance between two residues in Angstroms

    Calcualtes the distance between the center of two residues. The center is denoted
        by the average x,y,z position of three passed atoms for each residue.

    Args:
        res1_num (int) : the residue number of the first residue (PDB Column 5)
        res2_num (int) : the residue number of the second residue (PDB Column 5)
        res1_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 1 [("C2","C4","C6")]
        res2_atoms (tuple) : a tuple of the atom names of the three atoms whose position
            to average to find the center of residue 2 [("C2","C4","C6")]
    '''
    pass

def get_residue_distance_per_frame(pdb_filename : str) -> list:
    pass

if __name__=="__main__":
    ########JOB VARIABLES#######
    #load pdb file
    topology_filename = '5JUP_N2_wCCC_+1GCU_nowat_posttleap.pdb'

    trajectory_filename = '5JUP_N2_wCCC_+1GCU_nowat.mdcrd'
    ############################

    ########OPTIONAL VARS#######
    perspective_atom1_name = "C2"
    perspective_atom2_name = "C4"
    perspective_atom3_name = "C6"
    viewed_atom1_name = "C2"
    viewed_atom2_name = "C4"
    viewed_atom3_name = "C6"
    ############################
