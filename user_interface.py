def open_ui() -> None:
    '''Opens an environment where user can input commands to the command line

    Opens in the command line when user inputs the command 'stacker'
    '''
    pass

def parse_command(command : str) -> dict:
    '''Parses a user's input to the user interface

    Separates the user's command into the function, optional flags, and input

    Args:
        command (str) : command inputted to the UI
    Returns:
        parsed_input (dict) : dictionary of input and flags
    
    input: saveAmberParm -r 100 -s chr2 --outfile test.pdb input.pdb
    parsed_input: {command : save_AmberParm, r : 100, s : chr2, outfile : test.pdb, input : input.pdb}
    '''
    pass

