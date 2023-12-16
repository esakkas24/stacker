import argparse
import sys, os
from file_manipulation import filter_traj_to_pdb, filter_traj
from residue_movement import write_bottaro_to_csv
from pairwise_distance import calculate_residue_distance, get_residue_distance_for_frame
from visualization import display_arrays_as_video

class InvalidRoutine(Exception):
    pass

class ResEmpty(Exception):
    pass

class AtomEmpty(Exception):
    pass

class FrameEmpty(Exception):
    pass

def block_printing():
    '''Disable printing to standard output
    
    Adapted from https://stackoverflow.com/a/8391735'''
    sys.stdout = open(os.devnull, 'w')

def enable_printing():
    '''Enable printing to standard output
    
    Adapted from https://stackoverflow.com/a/8391735'''
    sys.stdout = sys.__stdout__

def run_python_command() -> None:
    '''Reads the user's passed in command line and runs the command

    Reads the command line input, runs the associated command with the
        added flags.
    '''
    parser = argparse.ArgumentParser()

    # Organizes all possible arguments
    ## Determines which script/subroutine is to be run
    parser.add_argument("-s", "--script", metavar="ROUTINE", help="Name of command to use (eg. filter_traj, bottaro, res_distance, pairwise)", required=True, default='')

    ## filter_traj requirements (--script filter_traj)
    parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=False)
    parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=False)
    parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath to output to, or prefix of output file if multiple outputs expected.", required=False)
    parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Comma-separated list of residues, also accepts colon (:) list creation (eg. 1:6 = 1,2,3,4,5)", required=False)
    parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names", required=False, default="C2,C4,C6")

    ## bottaro requirements (--script bottaro)
    parser.add_argument("-pdb", "--pdb_input", metavar="PDB_INPUT", help="Filepath to input PDB file containing two residues, the perspective and viewed nucleotide", required=False)
    parser.add_argument("-p", "--pers_res", metavar="PERSPECTIVE_RES", help="Comma-separated list of atomnames to use from residue 1 to find center of geometry for perspective nucleotide", required=False)
    parser.add_argument("-v", "--view_res", metavar="VIEWED_RES", help="Comma-separated list of atomnames to use from residue 2 to find center of geometry for viewed nucleotide", required=False)

    ## res_distance requirements (--script res_distance)
    parser.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze", required=False)

    ## pairwise requirements (--script pairwise)
    parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", help="Comma-separated list of 1-indexed Frame Numbers within trajectory to analyze", required=False)

    global args;
    args = parser.parse_args()

    convert_to_python_command()

def convert_to_python_command() -> None:
    '''Converts a parsed command to use to the correct subroutine and runs the routine

    Converts the specified script to a python command and runs it with the associated inputs
        based on the flags.
    '''
    command = args.script

    if command == 'filter_traj':
        filter_traj_routine()
    elif command == 'bottaro':
        bottaro_routine()
    elif command == 'res_distance':
        res_distance_routine()
    elif command == 'pairwise':
        pairwise_routine()
    else:
        raise InvalidRoutine(args.script + " is not a valid routine")
    
def filter_traj_routine() -> None:
    '''Runs the filter input trajectory to PDB routine
    
    Uses the passed in flags to run the filter_traj_to_pdb() script
        with the determined inputs from passed in flags
        
    Example Usage:
        [user]$ python3 stacker.py -s filter_traj -trj first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top 5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 425,426 -a C2,C4,C6
    '''
    if args.residues is not None:
        residues_desired = {res.strip() for res in args.residues.split(",")}
    else:
        ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names is not None:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        AtomEmpty("Must include a list of atom names to keep in the trajectory")

    create_parent_directories(args.output)
    filter_traj_to_pdb(trajectory_filename=args.trajectory, topology_filename=args.topology, output_pdb_filename=args.output, residues_desired=residues_desired, atomnames_desired=atomnames_desired)

def bottaro_routine() -> None:
    '''Runs the residue movement routine to create a CSV of the Bottaro values for a trajectory
    
    Runs the residue movement subroutine. Parses the inputs and creates a CSV containing the 
        r, rho, and theta values for each frame of the PDB trajectory between the two residues.
        
    Example Usage:
        [user]$ python3 stacker.py -s bottaro -pdb 5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -o command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv -p C2,C4,C6 -v C2,C4,C6
        [user]$ python3 stacker.py -s bottaro -pdb 5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p C2,C4,C6 -v C2,C4,C6
    '''
    if args.pers_res is not None:
        perspective_atom_names = {res.strip() for res in args.pers_res.split(",")}
    else:
        AtomEmpty("Must include a list of atom names to define Pespective Residue center of geometry")
    
    if args.view_res is not None:
        viewed_atom_names = {res.strip() for res in args.view_res.split(",")}
    else:
        AtomEmpty("Must include a list of atom names to define Viewed Residue center of geometry")

    create_parent_directories(args.output)
    write_bottaro_to_csv(pdb_filename=args.pdb_input, output_csv_name=args.output, res1_atom_names=perspective_atom_names, res2_atom_names=viewed_atom_names)

def res_distance_routine() -> None:
    '''Runs the Residue distance routine to determine the distance between the center of masses of two given residues
    
    Example Usage:
        [user]$ python3 stacker.py -s res_distance -trj first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top 5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -f 1 --residues 425,426 --atom_names C2,C4,C6
        '''
    if args.residues is not None:
        residues_desired = {res.strip() for res in args.residues.split(",")}
    else:
        ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names is not None:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        AtomEmpty("Must include a list of atom names to keep in the trajectory")

    block_printing()
    filtered_trj = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired, atomnames_desired=atomnames_desired)
    trj_frame = filtered_trj[args.frame-1]

    residues_desired = list(residues_desired)
    # Correct that calculate_residue_distance res_nums are 1-indexed
    distance_vector = calculate_residue_distance(trajectory=trj_frame, res1_num=int(residues_desired[0])+1, res2_num=int(residues_desired[1])+1, res1_atoms=tuple(atomnames_desired), res2_atoms=tuple(atomnames_desired))
    enable_printing()
    print(distance_vector.magnitude())

def pairwise_routine() -> None:
    '''Runs the Pairwise distance routine to create a single pairwise distance matrix for a chosen frame

    Runs the routine to create a pairwise distance matrix for a passed in trajectory and frame.

    Example Usage:
        [user]$ python3 stacker.py -s pairwise -trj first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top 5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 1,2,3,4,5,6,7,8,9,10 -fl 1,2 -o command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_nowat_pairwise_
    '''
    if args.residues is not None:
        residues_desired = {int(res.strip()) for res in args.residues.split(",")}
    else:
        residues_desired = {}

    if args.frame_list is not None:
        frame_list = [int(frame.strip()) for frame in args.frame_list.split(",")]
    else:
        FrameEmpty("Must include a list of frames to analyze in the trajectory")    
    
    trj_sub = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired)
    residues_desired = list(residues_desired)
    frames = [get_residue_distance_for_frame(trj_sub, i) for i in frame_list]
    create_parent_directories(args.output)
    display_arrays_as_video(frames, residues_desired, seconds_per_frame=1, outfile_prefix=args.output)

def create_parent_directories(outfile_prefix : str) -> None:
    '''Creates necessary parent directories to write an outfile given a prefix'''
    dir_name = os.path.dirname(outfile_prefix)
    os.makedirs(dir_name, exist_ok=True)

if __name__ == '__main__':
    run_python_command()

