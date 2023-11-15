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
    
    Example:
        input: saveAmberParm -r 100 -s chr2 --outfile test.pdb input.pdb
        parsed_input: {command : save_AmberParm, r : 100, s : chr2, outfile : test.pdb, input : input.pdb}
    '''
    pass

def convert_to_python_command(parsed_input : dict) -> tuple:
    '''Converts a parsed input to a python function from the module

    Takes parsed command line input and converts it to the recognized function
        from the module with the appropriate inputs.

    Args:
        parsed_input (dict) : dictionary of input and flags. Output of parse_command()
    Returns:
        function_and_agrs (tuple) : tuple of (function, args)
            function (Python Function) : function to be run
            args (dict) : dictionary where key is the function argument name and value is the value to be inputted

    Example:
        input: {command : save_AmberParm, r : 100, s : chr2, outfile : test.pdb, input : input.pdb}
        ouput: (save_amber_parm, {region : 100, section : chr2, outfile : test.pdb, input : input.pdb})
    '''

