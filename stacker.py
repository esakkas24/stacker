import argparse
import sys, os
from file_manipulation import *
from residue_movement import write_bottaro_to_csv
from pairwise_distance import calculate_residue_distance, get_residue_distance_for_frame, get_frame_average, get_top_stacking
from visualization import display_arrays_as_video, visualize_two_residue_movement_heatmap, visualize_two_residue_movement_scatterplot

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
    parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.RawTextHelpFormatter, 
                                     description="Wrapper to run stacker subroutines using the -s flag.\n" + \
                                        "More info on each routine given by `python stacker.py -s ROUTINE -h`")
    global args;
    args, remaining_args = parser.parse_known_args()

    # help when no script specified
    if '-s' not in remaining_args and '--script' not in remaining_args and ('--help' in remaining_args or '-h' in remaining_args): 
        parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. OPTIONS:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro:\n\tCreate polar plots like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  pairwise:\n\tCreate a stacking fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n   ", required=True, default='',
                            choices=['filter_traj', 'bottaro', 'res_distance', 'pairwise', 'stack_events'])
        parser.add_argument("-h", "--help", help="show this help message and exit", action='help')
        args = parser.parse_args()

    parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. OPTIONS:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro:\n\tCreate polar plots like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  pairwise:\n\tCreate a stacking fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n   ", required=True, default='',
                            choices=['filter_traj', 'bottaro', 'res_distance', 'pairwise', 'stack_events'])  
      
    args, remaining_args = parser.parse_known_args()

    # Organizes all possible arguments
    ## Determines which script/subroutine is to be run
    if args.script == 'filter_traj':
        parser.description = 'Filters trajectory and topology files to desired residue numbers and atom names and outputs to a PDB\n\nExamples:\n' +\
                            '[user]$ python3 stacker.py -s filter_traj -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 426,427 -a C2,C4,C6'
        
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        required_group.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath of PDB to output to", required=True)
        
        # optional arguments
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names to filter", required=False, default="C2,C4,C6")

    if args.script == 'bottaro':
        parser.description = 'Create polar plots of the movement of a "viewed residue" from the perspective of a "perspective residue"\nlike those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972). Creates CSV of these values' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ python3 stacker.py -s bottaro -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt scatter\n' +\
                            '\n[user]$ python3 stacker.py -s bottaro -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt heat'
        
        required_group = parser.add_argument_group('Required Arguments')
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-pdb", "--pdb_input", metavar="PDB_INPUT", help="If trajectory provided: filepath to intermediary PDB file containing two residues, the perspective and viewed nucleotide.\nIf no trajectory given, PDB is expected to already be 2-residue (use -s filter_traj if needed).\nIf empty, will use the same prefix as the trajectory file", required=False, default = '')
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath to output Bottaro values (frame, r, rho, theta) to.", required=False)
        required_group.add_argument("-p", "--pers_res", metavar="PERSPECTIVE_RES", help="residue index of the perspective residue whose plane to project onto. 0-/1-indexed based on -i flag (default: 1-indexed)", required=False)
        required_group.add_argument("-v", "--view_res", metavar="VIEWED_RES", help="residue index of the viewed residue whose midpoint will be projected onto perspective residue plane. 0-/1-indexed based on -i flag (default: 1-indexed)", required=False)
        required_group.add_argument("-pa", "--pers_atoms", metavar="PERSPECTIVE_ATOMS", help="Comma-separated list of atomnames to use from residue 1 to find center of geometry for perspective nucleotide", required=False)
        required_group.add_argument("-va", "--view_atoms", metavar="VIEWED_ATOMS", help="Comma-separated list of atomnames to use from residue 2 to find center of geometry for viewed nucleotide", required=False)
        parser.add_argument("-i", "--index", metavar="INDEX", type=int, help="index (0-index or 1-index) for perspective/viewed residue numbers (default: 1-indexed)", required=False, default = 1)
        parser.add_argument("-pt", "--plot_type", metavar="PLOT_TYPE", choices = ['scatter', 'heat', ''], help="plot type (scatter or heat) to visualize Bottaro values. If empty string, then just write to csv with no visualization", required=False, default = '')
        parser.add_argument("-po", "--plot_outfile", metavar="PLOT_OUTFILE", help="filename to output plot png to. If empty string, outputs to standard Python vis", required=False, default = '')

    if args.script == 'res_distance':
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, must provide only 2 residues, accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=True, action = SmartIndexingAction)
        parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names. Three required to get center of geometry for a residue", required=False, default="C2,C4,C6")

    if args.script == 'pairwise':
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze, if empty all frames are used", required=False, action=SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Prefix of output file if multiple outputs expected. If empty, will output displays to Python visual", default = '', required=False)

    if args.script == 'stack_events':
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues to subset trajectory, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Output CSV to write top stacking events to. If empty, will output displays to standard output", default = '', required=False)
        parser.add_argument("-n", "--n_events", type = int, metavar="N_EVENTS", help="Number of stacking events to display. If -1 display all events", default = '', required=False)
        frame_group = parser.add_mutually_exclusive_group()
        frame_group.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze, cannot be used with -fl", required=False)
        frame_group.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -fl", required=False, action=SmartIndexingAction)

    # help for specific scripts
    if '--help' in remaining_args or '-h' in remaining_args:
        parser.add_argument("-h", "--help", help="show this help message and exit", action='help')
        args = parser.parse_args()

    args = parser.parse_args()
    convert_to_python_command()

class SmartIndexingAction(argparse.Action):
    '''
    Custom argparse action to handle smart indexing of frame numbers.

    Parses a comma-separated list of frame numbers with optional ranges (e.g., '1-20, 34, 25, 50-100')
    and generates a list of individual frame numbers. Modifies the namespace by setting the attribute specified by the 'dest' parameter to the
    list of individual frame numbers.

    Args:
        parser: argparse.ArgumentParser
            The argparse parser object.
        namespace: argparse.Namespace
            The argparse namespace.
        values: str
            The input string containing frame numbers and ranges.
        option_string: str, default=None
            The option string.

    Returns:
        None

    Examples:
        >>> parser = argparse.ArgumentParser()
        >>> parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze", required=False, action=SmartIndexingAction)
        >>> args = parser.parse_args(["-fl", "1-20,34,25,50-100"])
        >>> print(args.frame_list)
        [1, 2, ..., 20, 34, 25, 50, 51, ..., 100]
    '''
    def __call__(self, parser, namespace, values, option_string=None):
        frame_list = []
        for item in values.split(','):
            if '-' in item:
                start, end = map(int, item.split('-'))
                frame_list.extend(range(start, end + 1))
            else:
                frame_list.append(int(item))
        frame_list.sort()
        setattr(namespace, self.dest, frame_list)

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
    elif command == 'stack_events':
        stack_events_routine()
    else:
        raise InvalidRoutine(args.script + " is not a valid routine")
    
def filter_traj_routine() -> None:
    '''Runs the filter input trajectory to PDB routine
    
    Uses the passed in flags to run the filter_traj_to_pdb() script
        with the determined inputs from passed in flags
        
    Example Usage:
        [user]$ python3 stacker.py -s filter_traj -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o testing/command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 426,427 -a C2,C4,C6
    '''
    if args.residues:
        residues_desired = set(args.residues)
    else:
        raise ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to keep in the trajectory")

    create_parent_directories(args.output)
    filter_traj_to_pdb(trajectory_filename=args.trajectory, topology_filename=args.topology, output_pdb_filename=args.output, residues_desired=residues_desired, atomnames_desired=atomnames_desired)

def bottaro_routine() -> None:
    '''Runs the residue movement routine to create a CSV of the Bottaro values for a trajectory
    
    Runs the residue movement subroutine. Parses the inputs and creates a CSV containing the 
        r, rho, and theta values for each frame of the PDB trajectory between the two residues.
        
    Example Usage:
        [user]$ python3 stacker.py -s bottaro -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt scatter
        [user]$ python3 stacker.py -s bottaro -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt heat
        '''
    if args.pdb_input == '':
        trj_prefix = args.trajectory.rsplit('.', 1)[0] 
        pdb_filename = trj_prefix + '.pdb'
    else:
        pdb_filename = args.pdb_input

    if args.pers_atoms:
        perspective_atom_names = {res.strip() for res in args.pers_atoms.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to define Pespective Residue center of geometry")
    
    if args.view_atoms:
        viewed_atom_names = {res.strip() for res in args.view_atoms.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to define Viewed Residue center of geometry")

    if args.pers_res is not None:
        pers_res_num = int(args.pers_res)
    else:
        raise AtomEmpty("Must include a 1-indexed residue index for the perspective residue")
    
    if args.view_res is not None:
        view_res_num = int(args.view_res)
    else:
        raise AtomEmpty("Must include a 1-indexed residue index for the perspective residue")

    create_parent_directories(pdb_filename)
    if args.trajectory and args.topology:
        filter_traj_to_pdb(trajectory_filename=args.trajectory, topology_filename=args.topology, output_pdb_filename=pdb_filename,
                           residues_desired={pers_res_num,view_res_num}, atomnames_desired=perspective_atom_names.union(viewed_atom_names))
    
    create_parent_directories(args.output)

    write_bottaro_to_csv(pdb_filename=pdb_filename, 
                         output_csv_name=args.output, perspective_residue_num=pers_res_num, viewed_residue_num=view_res_num,
                         res1_atom_names=tuple(perspective_atom_names), 
                         res2_atom_names=tuple(viewed_atom_names), index = args.index)
    
    if args.plot_type == 'heat':
        visualize_two_residue_movement_heatmap(args.output, plot_outfile=args.plot_outfile)
    elif args.plot_type == 'scatter':
        visualize_two_residue_movement_scatterplot(args.output, plot_outfile=args.plot_outfile)

def res_distance_routine() -> None:
    '''Runs the Residue distance routine to determine the distance between the center of masses of two given residues
    
    Example Usage:
        [user]$ python3 stacker.py -s res_distance -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -f 2 --residues 426,427 --atom_names C2,C4,C6
        '''
    if len(args.residues) != 2:
        raise ResEmpty("Must include only 2 residues")
    elif args.residues:
        residues_desired = set(args.residues)
    else:
        raise ResEmpty("Must include a list of residues to keep in the trajectory")

    if args.atom_names is not None:
        atomnames_desired = {atom.strip() for atom in args.atom_names.split(",")}
    else:
        raise AtomEmpty("Must include a list of atom names to keep in the trajectory")

    block_printing()
    filtered_trj = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired, atomnames_desired=atomnames_desired)
    trj_frame = filtered_trj[args.frame-1]

    residues_desired = list(residues_desired)
    # Correct that calculate_residue_distance res_nums are 1-indexed
    distance_vector = calculate_residue_distance(trajectory=trj_frame, res1_num=int(residues_desired[0]), res2_num=int(residues_desired[1]), res1_atoms=tuple(atomnames_desired), res2_atoms=tuple(atomnames_desired))
    enable_printing()
    print(distance_vector.magnitude())

def pairwise_routine() -> None:
    '''Runs the Pairwise distance routine to create a single pairwise distance matrix for a chosen frame

    Runs the routine to create a pairwise distance matrix for a passed in trajectory and frame.

    Example Usage:
        [user]$ python3 stacker.py -s pairwise -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2 
        '''
    if args.residues:
        residues_desired = set(args.residues)
    else:
        residues_desired = {}

    if args.frame_list:
        frame_list = args.frame_list
    else:
        frame_list = []

    trj_sub = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired)
    if frame_list:
        frames = [get_residue_distance_for_frame(trj_sub, i) for i in frame_list]
    else:
        frames = [get_residue_distance_for_frame(trj_sub, i) for i in range(0,trj_sub.n_frames)]
    create_parent_directories(args.output)
    display_arrays_as_video(frames, list(residues_desired), seconds_per_frame=1, outfile_prefix=args.output)

def stack_events_routine() -> None:
    '''Runs the routine to get the residue pairings with the most pi stacking (center of geometry distance closest to 3.5Å)

    Example Usage:
    [user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -f 1 -n 5
    [user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-100 -fl 1-10 -n 5
    '''
    if args.residues:
        residues_desired = set(args.residues)
    else:
        residues_desired = {}

    trj_sub = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired)

    if args.frame:
        frame = get_residue_distance_for_frame(trj_sub, frame = args.frame)
    elif args.frame_list:
        frames = [get_residue_distance_for_frame(trj_sub, frame_i) for frame_i in args.frame_list]
        frame = get_frame_average(frames)

    get_top_stacking(trj_sub, frame, output_csv = args.output, n_events = args.n_events)

def create_parent_directories(outfile_prefix : str) -> None:
    '''Creates necessary parent directories to write an outfile given a prefix'''
    dir_name = os.path.dirname(outfile_prefix)
    if dir_name == '': dir_name = '.'
    os.makedirs(dir_name, exist_ok=True)

if __name__ == '__main__':
    run_python_command()

