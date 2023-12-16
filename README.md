# StACKER
**St**acking **A**nalysis and **C**onformational **K**inetics for **E**xamining **R**esidues

Analyzes pi-stacking interactions between residues in [Molecular Dynamics (MD)](https://github.com/esakkas24/stacker/docs/background.md) trajectories.

Developed by Eric Sakkas ([esakkas@wesleyan.edu](mailto:esakkas@wesleyan.edu)) in the [Weir Lab](https://weirlab.research.wesleyan.edu/) at Wesleyan University.

Runs on Python 3.5.x and has the following dependencies: [mdtraj](https://www.mdtraj.org/1.9.8.dev0/index.html), [pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/doc/stable/index.html), [matplotlib](https://matplotlib.org/stable/), [seaborn](https://seaborn.pydata.org/index.html)

## Overview

Manipulates the outputs of an MD simulation and analyzes the pi-stacking interactions. Creates a "Pi-Stacking Fingerprint" for a structure at each frame. Presents Pi-stacking interactions between two residues through analysis of their relevant movement.

## Installation Instructions

### Clone StACKER repository to local computer

In the command line, run:
```
git clone https://github.com/esakkas24/stacker.git
```

### Download mdtraj
Installing mdtraj presents issues on the newest version of pip. To fix this, run the script get-pip.py to download an older version of pip:
```
python3 get-pip.py
```
The output should come with a warning showing the location of the new pip version:
```
WARNING: The scripts pip, pip3, and pip3.8 are installed in '/Users/ericsakkas/Library/Python/3.8/bin' which is not on PATH
```

Use this path to install `mdtraj`:
```
/Users/ericsakkas/Library/Python/3.8/bin/pip3 install mdtraj
```

If successful, this will also install the NumPy dependency:

```
Successfully installed astunparse-1.6.3 mdtraj-1.9.9 numpy-1.24.4 ...
```

### Install Pandas
```
pip3 install pandas
```
### Install matplotlib
```
pip3 install matplotlib
```
### Install seaborn
```
pip3 install seaborn
```
## Testing Features

All features can be tested by running the unit tests at the end of each Python script, or by running StACKER in the command line. All tests are explained in the testing.md file.

MD Files are provided for testing convenience:
- first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd : A 10-frame trajectory file
- 5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop : The associated Topology File with the above trajectory.
- 5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd_3200frames.pdb : A larger Trajectory combined with a Topology file with 3200 frames.

## Proposed Features

1) A frontend UI where users can call "stacker" in the command line and pull up an interface to input stacker-specific functions (similar to how the "python" command works in the Terminal)
2) A Vector Class to compute distances within the 3D space of the MD simulation.
3) Users can convert the .trj output of an MD simulation (which contains atom position, velocities, and forces per frame with no info on atom identity) to a .pdb file (which has atom idenity, position, the residue they make up, and more).
4) Users can input an MD simulation and two residues and get a map of how those two residues move relative to each other ([Figure D](https://www.mdpi.com/ijms/ijms-23-01417/article_deploy/html/images/ijms-23-01417-g005.png) shows a heatmap of how one residue moves in the perspective of another).
5) Users can get a "Stacking Fingerprint" for a given frame of a structure. A Stacking Fingerprint involves checking every pair of residues for stacking interactions between the two residues, and returning a pairwise comparison as a Matrix.
6) A visualization interface that allows for the display of residue-residue movement from Feature 4 (as in Figure D) and the display of a heatmap for stacking interactions as in Feature 5
    - Heatmap: matrix where x-axis and y-axis are residue idenities within the strucutre (eg. residue 1, residue 48, etc.) and matrix(i,j) is colored by distance between the residues, where stacking interactions occur at around 3-4 Angstroms apart.

## Stakeholders and Intended Users

The package is intended to be used by anyone in academic or computational biology centers running molecular dynamics. No explicit prerequisites, but a conceptual understanding of MD output files is beneficial. The stakeholders include the users, other researchers reliant on the computational data, and any beneficiaries of the overall research conducted using StACKER.

