import argparse
import sys, os
import numpy as np
import pandas as pd
from file_manipulation import *
from residue_movement import write_bottaro_to_csv
from pairwise_distance import calculate_residue_distance, get_residue_distance_for_frame, get_residue_distance_for_trajectory, get_frame_average, get_top_stacking
from visualization import display_arrays_as_video, visualize_two_residue_movement_heatmap, visualize_two_residue_movement_scatterplot
import random

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
    
    # If no flags specified at all
    if not any(vars(args)) and not remaining_args:
        print('usage: stacker.py -s ROUTINE [-h]\n\n' + \
            'Wrapper to run stacker subroutines using the -s flag.\n' + \
            'More info on each routine given by `python stacker.py -s ROUTINE -h` or `python stacker.py -s ROUTINE --help`\n\n' + \
            'options:\n' +\
            '-s ROUTINE, --script ROUTINE\n' +\
            '            Name of command to use. Options for ROUTINE:\n\n' + \
            '              filter_traj:\n' +\
            '                    filters trajectory and topology files to desired residue numbers and atom names\n' + \
            '              bottaro OR pairwise OR psf:\n' +\
            '                    Create a Pairwise Stacking Fingerprint (PSF): polar plots like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n' + \
            '              res_distance:\n' + \
            '                    Get the distance between two residues in a given frame\n' + \
            '              system OR ssf:\n' + \
            '                    Create a System Stacking Fingerprint (SSF) with pairwise distances for each residue\n' + \
            '              stack_events:\n' + \
            '                    Get list of residues with most stacking events (distance closest to 3.5Å)\n' + \
            '              compare:\n' +\
            '                    Get the most changed stacking events between two fingerprints using the outputs of python stacker.py -s stack_events\n\n' +\
            '-h, --help            show this help message and exit\n')
        return

    # help when no script specified
    if ('-s' not in remaining_args and '--script' not in remaining_args and ('--help' in remaining_args or '-h' in remaining_args)): 
        parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. Options for ROUTINE:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro OR pairwise OR psf:\n\tCreate Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  system OR ssf:\n\tCreate a System Stacking Fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n" +\
                            "  compare:\n\tGet the most changed stacking events between two fingerprints using the outputs of python stacker.py -s stack_events\n",
                             required=True, default='', choices=['filter_traj', 'bottaro', 'pairwise', 'psf', 'res_distance', 'system', 'ssf', 'stack_events', 'compare'])
        parser.add_argument("-h", "--help", help="show this help message and exit", action='help')
        args = parser.parse_args()

    parser.add_argument("-s", "--script", metavar="ROUTINE", help='Name of command to use. Options for ROUTINE:\n\n' + \
                            "  filter_traj:\n\tfilters trajectory and topology files to desired residue numbers and atom names\n" + \
                            "  bottaro OR pairwise OR psf:\n\tCreate Polar Stacking Fingerprint like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)\n" + \
                            "  res_distance:\n\tGet the distance between two residues in a given frame\n" +\
                            "  system OR ssf:\n\tCreate a System Stacking Fingerprint of distances by residue\n" + \
                            "  stack_events:\n\tGet list of residues with most stacking events (distance closest to 3.5Å)\n" +\
                            "  compare:\n\tGet the most changed stacking events between two fingerprints using the outputs of python stacker.py -s stack_events\n",
                             required=True, default='', choices=['filter_traj', 'bottaro', 'pairwise', 'psf', 'res_distance', 'system', 'ssf', 'stack_events', 'compare'])
      
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

    if args.script == 'bottaro' or args.script == 'pairwise' or args.script == 'psf':
        parser.description = 'Create Polar Stacking Fingerprint (PSF) of the movement of a "viewed residue" from the perspective of a "perspective residue"\nlike those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972). Creates CSV of these values' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ python3 stacker.py -s bottaro -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt scatter\n' +\
                            '\n[user]$ python3 stacker.py -s bottaro -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt heat'
        
        required_group = parser.add_argument_group('Required Arguments')
        parser.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation, if empty then 2-residue PDB expected", required=False, default = '')
        parser.add_argument("-pdb", "--pdb_input", metavar="PDB_INPUT", help="If trajectory provided: filepath to intermediary PDB file containing two residues, the perspective and viewed nucleotide.\nIf no trajectory given, PDB is expected to already be 2-residue (use -s filter_traj if needed).\nIf empty, will use the same prefix as the trajectory file", required=False, default = '')
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filepath to output Bottaro values (frame, r, rho, theta) to. If empty, will use the same prefix as the trajectory file.", required=False, default = '')
        required_group.add_argument("-p", "--pers_res", metavar="PERSPECTIVE_RES", help="residue index of the perspective residue whose plane to project onto. 0-/1-indexed based on -i flag (default: 1-indexed)", required=True)
        required_group.add_argument("-v", "--view_res", metavar="VIEWED_RES", help="residue index of the viewed residue whose midpoint will be projected onto perspective residue plane. 0-/1-indexed based on -i flag (default: 1-indexed)", required=True)
        parser.add_argument("-pa", "--pers_atoms", metavar="PERSPECTIVE_ATOMS", help="Comma-separated list of atomnames to use from residue 1 to find center of geometry for perspective nucleotide", required=False, default="C2,C4,C6")
        parser.add_argument("-va", "--view_atoms", metavar="VIEWED_ATOMS", help="Comma-separated list of atomnames to use from residue 2 to find center of geometry for viewed nucleotide", required=False, default="C2,C4,C6")
        parser.add_argument("-i", "--index", metavar="INDEX", type=int, help="index (0-index or 1-index) for perspective/viewed residue numbers (default: 1-indexed)", required=False, default = 1)
        parser.add_argument("-pt", "--plot_type", metavar="PLOT_TYPE", choices = ['scatter', 'heat', ''], help="plot type (scatter or heat) to visualize Bottaro values. If empty string, then just write to csv with no visualization", required=False, default = '')
        parser.add_argument("-po", "--plot_outfile", metavar="PLOT_OUTFILE", help="filename to output plot png to. If empty string, outputs to standard Python vis", required=False, default = '')
        parser.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)
        parser.add_argument("-n", "--no_inter", help="Delete intermediate files after command runs", action = 'store_true', default=False)

    if args.script == 'res_distance':
        parser.description = 'Get the distance between two residues in a given frame\n\n' + \
                                'Examples:\n' +\
                                '[user]$ python3 stacker.py -s res_distance -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -f 2 --residues 426,427 --atom_names C2,C4,C6'
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze", required=False)
        required_group.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, must provide only 2 residues, accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=True, action = SmartIndexingAction)
        parser.add_argument("-b", "--bootstrap", metavar="N_FRAMES", help="Run bootstrap analysis on this residue pairing, sampling N_FRAMES with replacement", required=False, type = int)
        parser.add_argument("-a", "--atom_names", metavar="ATOM_NAMES", help="Comma-separated list of atom names. Three required to get center of geometry for a residue. default = C2,C4,C6", required=False, default="C2,C4,C6")

    if args.script == 'system' or args.script == 'ssf':
        parser.description = 'Creates a System Stacking Fingerprint of the average structure across the chosen frames of a trajectory.' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ python stacker.py -s system -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2\n' +\
                                '\n[user]$ python stacker.py -s ssf -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2 -g 10 -o testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_nowat_pairwise_avg_1to2.png -d testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt\n'
        
        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-i", "--input", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information, in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `python stacker.py -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        frame_group = parser.add_mutually_exclusive_group()
        frame_group.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze, cannot be used with -fl", required=False)
        frame_group.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Filename of output PNG to write plot to. If empty, will output displays to Python visual", default = '', required=False)
        parser.add_argument("-g", "--get_stacking", metavar="N_EVENTS", help="Get list of N_EVENTS residues with most stacking events (distance closest to 3.5Å) in the average structure across all frames.\nPrint to standard output. Equivalent to -s stack_events -n N_EVENTS", type = int, required=False, default = -1)
        parser.add_argument("-d", "--data_output", metavar="OUTPUT_FILE", help="Output the calculated per-frame numpy arrays that create the stacking fingerprint matrix to a file", default = '', required=False)
        parser.add_argument("-B", "--input_B", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information for a second fingerprint, creates fingerprint where top left is initial input, bottom right is second fingerprint.\n Used in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `python stacker.py -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        parser.add_argument("-l", "--limits", metavar="LIMITS", help="limits of the color scale, default = (0,7)", default = '0,7')
        parser.add_argument("-y", "--scale_style", metavar="SCALE_STYLE", help=" style of color scale. {bellcurve, gradient}", default = 'bellcurve')
        parser.add_argument("-t", "--threads", metavar="N_THREADS", help="Use multithreading with INT worker threads", type = int, required = False, default = 1)


    if args.script == 'stack_events':
        parser.description = 'Get list of residues with most stacking events (distance closest to 3.5Å) in the stacking fingerprint of the average structure across all frames of a trajectory' + \
                                '\n\nExamples:\n' +\
                                '\n[user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -f 1 -n 5\n' +\
                                '[user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-100 -fl 1-10 -n 5\n' +\
                                '[user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -n 5 -i testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt\n' 

        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-trj", "--trajectory", metavar="TRAJECTORY_FILENAME", help="Filepath to trajectory file for the MD simulation", required=True)
        required_group.add_argument("-top", "--topology", metavar="TOPOLOGY_FILENAME", help="Filepath to Topology file for the MD simulation", required=True)
        parser.add_argument("-r", "--residues", metavar="RESIDUES", help="Smart-indexed list of 1-indexed residues to subset trajectory, also accepts dash (-) list creation (eg. 1-5,10 = 1,2,3,4,5,10)", required=False, action = SmartIndexingAction)
        parser.add_argument("-o", "--output", metavar="OUTPUT_FILE", help="Output tab-separated txt file to write top stacking events to. If empty, will output displays to standard output", default = '', required=False)
        parser.add_argument("-n", "--n_events", type = int, metavar="N_EVENTS", help="Number of stacking events to display. If -1 display all events", default = -1, required=False)
        parser.add_argument("-i", "--input", metavar="INPUT_FILE", help="Input .txt file containing per-frame stacking information, in lieu of running stacking fingerprint analysis again.\nTXT file can be created by running `python stacker.py -s system -d OUTPUT_FILE`\n-r flag must match the residues used to create the TXT file")
        parser.add_argument("-j", "--include_adjacent", help="Boolean whether to include adjacent residues in the printed output", action = 'store_true', default=False)
        frame_group = parser.add_mutually_exclusive_group()
        frame_group.add_argument("-f", "--frame", type=int, metavar="FRAME_NUM", help="1-indexed Frame Number within trajectory to analyze, cannot be used with -fl", required=False)
        frame_group.add_argument("-fl", "--frame_list", metavar="FRAME_LIST", default='', help="Smart-indexed list of 1-indexed Frame Numbers within trajectory to analyze,\ngets average distance between residues across these frames\nif empty all frames are used, cannot be used with -f", required=False, action=SmartIndexingAction)

    if args.script == 'compare':
        parser.description = 'Print the most changed stacking events between two fingerprints using the outputs of python stacker.py -s stack_events' +\
                                '\n\nExamples:\n' +\
                                '[user]$ python stacker.py -s compare -A /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418986 -B /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418997 -SA _tUAG_aCUA_+1GCU -SB _tUAG_aCUA_+1CGU\n'

        required_group = parser.add_argument_group('Required Arguments')
        required_group.add_argument("-A", "--file_A", metavar="FILENAME_A", help = "Filepath to the output log of python stacker.py -s stack_events for the first stacking fingerprint", required = True)
        required_group.add_argument('-B', '--file_B', metavar="FILENAME_B", help = 'Filepath to the output log of python stacker.py -s stack_events for the second stacking fingerprint', required = True)
        required_group.add_argument('-SA', '--source_A', metavar="SOURCE_A", help = 'String describing source of file A, e.g. `_tUAG_aCUA_+1GCU`', required = True)
        required_group.add_argument('-SB', '--source_B', metavar="SOURCE_B", help = 'String describing source of file B, e.g. `_tUAG_aCUA_+1CGU`', required = True)

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
    elif command == 'bottaro' or command == 'pairwise' or command == 'psf':
        bottaro_routine()
    elif command == 'res_distance':
        res_distance_routine()
    elif command == 'system' or command == 'ssf':
        system_routine()
    elif command == 'stack_events':
        stack_events_routine()
    elif command == 'compare':
        compare_routine()
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
        [user]$ python3 stacker.py -s pairwise -pdb testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb -o testing/command_line_tests/bottaro/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv -p 426 -v 427 -pa C2,C4,C6 -va C2,C4,C6 -pt heat
        '''
    trj_prefix = args.trajectory.rsplit('.', 1)[0] 

    if args.pdb_input == '':
        pdb_filename = trj_prefix + '.pdb'
    else:
        pdb_filename = args.pdb_input

    if args.output == '':
        prefix = pdb_filename.rsplit('.', 1)[0] 
        output_name = prefix + '.csv'
    else:
        output_name = args.output

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

    if args.frame_list:
        frame_list = set(args.frame_list)
    else:
        frame_list = {}

    create_parent_directories(pdb_filename)
    if args.trajectory and args.topology:
        filter_traj_to_pdb(trajectory_filename=args.trajectory, topology_filename=args.topology, output_pdb_filename=pdb_filename,
                           residues_desired={pers_res_num,view_res_num}, atomnames_desired=perspective_atom_names.union(viewed_atom_names))
    
    create_parent_directories(output_name)

    write_bottaro_to_csv(pdb_filename=pdb_filename, 
                         output_csv_name=output_name, perspective_residue_num=pers_res_num, viewed_residue_num=view_res_num,
                         res1_atom_names=tuple(perspective_atom_names), 
                         res2_atom_names=tuple(viewed_atom_names), index = args.index)
    
    if args.plot_type == 'heat':
        create_parent_directories(args.plot_outfile)
        visualize_two_residue_movement_heatmap(output_name, plot_outfile=args.plot_outfile, frame_list = frame_list)
    elif args.plot_type == 'scatter':
        create_parent_directories(args.plot_outfile)
        visualize_two_residue_movement_scatterplot(output_name, plot_outfile=args.plot_outfile, frame_list = frame_list)

    if args.no_inter:
        if args.trajectory and args.topology: #intermediate pdb created
            os.remove(pdb_filename)
        if os.path.exists(output_name):
            os.remove(output_name)

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
    if args.bootstrap:
        n_frames = len(filtered_trj)
        frames = [random.randint(0, n_frames-1) for _ in range(args.bootstrap)]
    else:
        frames = [args.frame-1]

    i = 0
    residues_desired = list(residues_desired)
    res_distances = []
    for frame in frames:
        trj_frame = filtered_trj[frame]
        # Correct that calculate_residue_distance res_nums are 1-indexed
        print(i)
        i+=1
        distance_vector = calculate_residue_distance(trajectory=trj_frame, res1_num=int(residues_desired[0]), res2_num=int(residues_desired[1]), res1_atoms=tuple(atomnames_desired), res2_atoms=tuple(atomnames_desired))
        enable_printing()
        res_distances.append(distance_vector.magnitude())

    if args.bootstrap:
        print("Bootstrap Mean Distance:", np.mean(res_distances))
        print("Bootstrap Standard Deviation:", np.std(res_distances))
        print("Bootstrap Percentile:", np.percentile(res_distances, [2.5,5.0, 95.0, 97.5]))
    else:
        print(res_distances[0])

def system_routine() -> None:
    '''Runs the System Stacking Fingerprint routine to create a single SSF for a chosen frame

    Runs the routine to create a SSF for a passed in trajectory and frame.

    Example Usage:
        [user]$ python3 stacker.py -s system -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2 
        [user]$ python stacker.py -s ssf -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -fl 1-2 -g 10 -o testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_nowat_pairwise_avg_1to2.png -d testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt
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

    if args.input:
        print("Loaded fingerprint data from:", args.input)
        loaded_arr = np.loadtxt(args.input)
        frames = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // trj_sub.n_residues, trj_sub.n_residues)
    elif args.frame_list:
        frames = get_residue_distance_for_trajectory(trj_sub, frame_list, threads = args.threads)
    elif args.frame:
        frames = np.array([get_residue_distance_for_frame(trj_sub, args.frame)])
    else:
        frames = get_residue_distance_for_trajectory(trj_sub, [i for i in range(1,trj_sub.n_frames+1)], threads = args.threads)

    if args.data_output:
        print(f"Number of SSFs made:{frames.shape}")
        frames_to_save = frames.reshape(frames.shape[0], -1)
        np.savetxt(args.data_output, frames_to_save)

    if args.limits:
        scale_limits = tuple(float(i) for i in args.limits.replace('(','').replace(')','').replace('...', '').split(','))

    avg_frames = [get_frame_average(frames)]

    if args.input_B:
        print("Loaded second fingerprint data from:", args.input_B)
        loaded_arr = np.loadtxt(args.input_B)
        frames_B = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // trj_sub.n_residues, trj_sub.n_residues)
        avg_frames_B = [get_frame_average(frames_B)]
        avg_frames = [combine_frames(avg_frames[0], avg_frames_B[0])]
        print(avg_frames)

    if args.get_stacking:
        get_top_stacking(trj_sub, avg_frames[0], output_csv = '', n_events = args.get_stacking)

    sorted_res = list(residues_desired)
    sorted_res.sort()

    create_parent_directories(args.output)
    display_arrays_as_video(avg_frames, sorted_res, seconds_per_frame=1, outfile=args.output, scale_limits=scale_limits, scale_style=args.scale_style)

def combine_frames(frames_A, frames_B):
    """Combines two 2D numpy arrays (frames_A and frames_B) into a single array.

    This function takes two 2D numpy arrays of the same shape and combines them
    into a new array. The upper triangular part (including the diagonal) of the 
    resulting array is filled with the corresponding elements from frames_A, 
    while the lower triangular part (including the diagonal) is filled with the 
    corresponding elements from frames_B.

    Args:
        frames_A : numpy.ndarray
            A 2D numpy array.
        frames_B : numpy.ndarray 
            A 2D numpy array of the same shape as frames_A.

    Returns:
        numpy.ndarray 
            A new 2D numpy array with combined elements from frames_A and frames_B.

    Example:
    >>> frames_A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> frames_B = np.array([[9, 8, 7], [6, 5, 4], [3, 2, 1]])
    >>> combine_frames(frames_A, frames_B)
    array([[1., 2., 3.],
           [6., 5., 6.],
           [3., 2., 9.]])
    """
    Am, An = frames_A.shape
    array_to_fill = np.zeros((Am,An))
    for i in range(Am):
        array_to_fill[i,i:] = frames_A[i,i:] 
    for j in range(An):
        array_to_fill[j:,j] = frames_B[j:,j]
    return array_to_fill

def stack_events_routine() -> None:
    '''Runs the routine to get the residue pairings with the most pi stacking (center of geometry distance closest to 3.5Å)

    Example Usage:
    [user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -f 1 -n 5
    [user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-100 -fl 1-10 -n 5
    [user]$ python stacker.py -s stack_events -trj testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -r 90-215 -n 5 -i testing/command_line_tests/pairwise/5JUP_N2_tUAG_aCUA_+1GCU_data_1to2.txt
    '''
    if args.residues:
        residues_desired = set(args.residues)
    else:
        residues_desired = {}

    trj_sub = filter_traj(trajectory_filename=args.trajectory, topology_filename=args.topology, residues_desired=residues_desired)

    if args.input:
        loaded_arr = np.loadtxt(args.input)
        frames = loaded_arr.reshape(loaded_arr.shape[0], loaded_arr.shape[1] // trj_sub.n_residues, trj_sub.n_residues)
        frame = get_frame_average(frames)
    elif args.frame:
        frame = get_residue_distance_for_frame(trj_sub, frame = args.frame)
    elif args.frame_list:
        frames = get_residue_distance_for_trajectory(trj_sub, args.frame_list, threads = args.threads)
        frame = get_frame_average(frames)

    get_top_stacking(trj_sub, frame, output_csv = args.output, n_events = args.n_events, include_adjacent = args.include_adjacent)

def compare_routine() -> None:
    '''Runs the routine to return the most changed stacking events—the residue pairs with the largest change in average distance between
        two trajectories.

    Example Usage:
    [user]$ python stacker.py -s compare -A /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418986 -B /home66/esakkas/STACKER/SCRIPTS/slurmLogs_fingerprint/out_fingerprint_2418997 -SA _tUAG_aCUA_+1GCU -SB _tUAG_aCUA_+1CGU
    '''
    def find_row_with_header(filename, header):
        '''Get the row number in file that matches a given header

        Given a filename that represents a spreadsheet, find the row number
            in the file that contains the passed header string.

        Args:
            filename : str
                file path to spreadsheet file
            header : str
                string representng a header
        Returns:
            row_number : int
                row number where the header appears in the file
        '''
        with open(filename, 'r') as file:
            for idx, line in enumerate(file):
                if line.strip() == header:
                    return idx
    
    def preprocess_df(df):
        '''Preprocess the DataFrame by sorting Row and Column alphabetically
        
        Args:
            df : Pandas Dataframe
                inputted df with at least columns Row, Column
        Returns:
            df : Pandas Dataframe
                inputted df with rows organized alphabetically at the Row, Column variable
        '''
        cols_to_sort = ['Row','Column']
        df = pd.concat([pd.DataFrame(np.sort(df[cols_to_sort].values), columns=cols_to_sort, index=df.index), df[df.columns[~df.columns.isin(cols_to_sort)]]], axis=1)
        return df
    
    def process_df(filename : str) -> pd.DataFrame:
        '''Read filename and discard the header to identify the entries.
        '''
        row_number = find_row_with_header(filename, header)
        data = pd.read_csv(filename, sep='\t', skiprows=row_number)
        data = preprocess_df(data)
        return data

    file1 = args.file_A
    file2 = args.file_B
    file1_source = args.source_A
    file2_source = args.source_B

    header = "Row\tColumn\tValue"

    data1 = process_df(file1)
    data2 = process_df(file2)

    merged_data = pd.merge(data1, data2, on=['Row', 'Column'], suffixes=[file1_source, file2_source], how='inner')
    merged_data['Discrepancy'] = abs(merged_data['Value' + file1_source] - merged_data['Value' + file2_source])

    subset_data = merged_data[(merged_data['Value' + file1_source] < 4) | (merged_data['Value' + file2_source] < 4)]
    subset_data = subset_data.sort_values(by='Discrepancy', ascending=False)
    subset_data = subset_data.rename(columns={'Value' + file1_source : 'AvgDist' + file1_source, 'Value' + file2_source : 'AvgDist' + file2_source})
    print(subset_data)

def create_parent_directories(outfile_prefix : str) -> None:
    '''Creates necessary parent directories to write an outfile given a prefix'''
    dir_name = os.path.dirname(outfile_prefix)
    if dir_name == '': dir_name = '.'
    os.makedirs(dir_name, exist_ok=True)

if __name__ == '__main__':
    run_python_command()