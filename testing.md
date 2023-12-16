# Testing

## Proposed Feature 1: User Interface

Many of the Features are capable of being run directly from the command line, with various flags indicating the specific functions and inputs to run. The basic method to run a command is to run the following command line prompt:

```
python3 stacker.py -s ROUTINE [FLAGS]
```

### Filter Trajectory

Users can filter a trajectory file and associated topology with multiple residues and atomnames to a smaller PDB file (combined trajectory and topology)

The command line prompt below takes the trajectory `first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd` and topology `5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop` which have around 500 Residues with 20 atoms each and filters to an outputed PDB `command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb` with just Residues 425 and 426 and just atoms C2, C4, and C6

```
[user]$ python3 stacker.py -s filter_traj -trj first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd -top 5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop -o command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb -r 425,426 -a C2,C4,C6
```

Looking into this file we see only residues 425/426 (Column 6) and only atoms C2, C4, and C6 (Column 3):

```
$ head command_line_tests/filter/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb
MODEL        0
ATOM      1  C6    G A 425      76.043  74.000  47.143  1.00  0.00           C  
ATOM      2  C2    G A 425      76.759  75.533  48.812  1.00  0.00           C  
ATOM      3  C4    G A 425      75.138  74.094  49.385  1.00  0.00           C  
ATOM      4  C6    C A 426      69.477  73.317  48.913  1.00  0.00           C  
ATOM      5  C4    C A 426      69.380  71.739  47.088  1.00  0.00           C  
ATOM      6  C2    C A 426      67.419  72.324  48.189  1.00  0.00           C  
TER       7        C A 426
ENDMDL
MODEL        1
```



