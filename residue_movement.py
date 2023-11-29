from __future__ import print_function
import math
import csv
import mdtraj as md
from vector import *

'''
Eric Sakkas - Weir Lab, Wesleyan University
Email : esakkas@wesleyan.edu

This code takes a pdb file with multiple frames and two residues. 
It then calculates the r, rho, and theta values at each frame
as defined in the Bottaro paper (), and exports these values to a .csv file.
'''

def collect_atom_locations_by_frame(pdb : md.Trajectory, residue_num : int, atom_id : str) -> list:
    '''Creates a list of all atom locations for a particular residue per frame

    Curates a list coords_by_frame where coords_by_frame[i] is the (x,y,z) positions at 
        of a provided atom_id of a residue residue_num at the ith frame.

    Args:
        topology (topology) : topology object of the pdb created by pdb.topology
        residue_num (int) : the residue number of the residue where atom_id is found (PDB Column 5)
        atom_id (str) : the name of the atom to get coordinates for (PDB Column 2)
    Returns:
        coords_by_frame (list) : list of (x,y,z) coordinates of atom_id in residue_num for each frame
    '''
    topology = pdb.topology
    number_of_frames = pdb.n_frames
    atomic_index = topology.select("name " + atom_id + " and residue " + str(residue_num))[0]
    coords_by_frame = [tuple(pdb.xyz[frame_idx, atomic_index,:]) for frame_idx in range(0, number_of_frames)]
    return coords_by_frame

def calc_center_3pts(a : Vector, b : Vector, c : Vector) -> Vector:
    '''Finds average x,y,z position of three x,y,z tuples
    
    Takes in three tuples generated using the pdb.xyz method and finds their center. Works 
    with three points that make up a triangle (like those within a 6-member ring). 
        
    Args: 
        (a1,a2,a3) (Vector) : x,y,z coordinates of first point
        (b1,b2,b3) (Vector) : x,y,z coordinates of second point
        (c1,c2,c3) (Vector) : x,y,z coordinates of third point
    Returns:
        midpoint (tuple) : one tuple with (x,y,z) coordinates at the center of the three input points.
    '''
    midpoint = Vector((a.x+b.x+c.x)/3,(a.y+b.y+c.y)/3,(a.z+b.z+c.z)/3)
    return midpoint

class Base:
    '''Represents a Base data type with coordinates for C2, C4, C6, and their midpoint at a single frame.

    This class defines a data type 'Base' that consists of coordinates for the atoms
    C2, C4, and C6, as well as the midpoint of these atoms for a single residue at a single frame.

    Attributes:
        c2_coords (Vector) : (x, y, z) coordinates for the atom C2.
        c4_coords (Vector) : (x, y, z) coordinates for the atom C4.
        c6_coords (Vector) : (x, y, z) coordinates for the atom C6.
        midpoint_coords (list) : (x, y, z) coordinates representing the midpoint of C2, C4, and C6.
    '''
    def __init__(self, c2_coords : tuple, c4_coords : tuple, c6_coords : tuple, midpoint_coords : tuple) -> None:
        '''Initialize a Base instance.

        Args:
            c2 (tuple): (x, y, z) coordinates for C2.
            c4 (tuple): (x, y, z) coordinates for C4.
            c6 (tuple): (x, y, z) coordinates for C6.
        '''
        self.c2_coords = Vector(*c2_coords)
        self.c4_coords = Vector(*c4_coords)
        self.c6_coords = Vector(*c6_coords)
        self.midpoint_coords = Vector(*midpoint_coords)

def create_base_from_coords_list(frame : int, C2_coords : list, C4_coords : list, C6_coords : list, midpoint_coords : list) -> list:
    '''Combines C2, C4, C6 positiions with midpoint positions for a given frame

    Takes a frame (indexed at 0) and outputs a list of the x,y,z locations of C2,C4,C6, and midpoint 
    of the same residue at that frame. 
    
    Args: 
        frame (int) : frame number index 0
        C2_coords (list) : list of (x,y,z) coordinates of C2 in some residue for each frame
        C4_coords (list) : list of (x,y,z) coordinates of C4 in some residue for each frame
        C6_coords (list) : list of (x,y,z) coordinates of C6 in some residue for each frame
        midpoint_coords (list) : list of (x,y,z) coordinates of the midpoint of residue for each frame
    Returns:
        coords_at_frame (list): list of 4 tuples; structure:[C2 location,C4 location,C6 location,midpoint location]
    '''
    midpoint_tuple = midpoint_coords[frame].x, midpoint_coords[frame].y, midpoint_coords[frame].z
    return Base(C2_coords[frame], C4_coords[frame], C6_coords[frame], midpoint_tuple)

def correct_theta_sign(rho : Vector, y_axis : Vector, theta : float) -> float:
    '''Corrects sign of Bottaro calculated theta to account for 3D Space

    When calculating the angle theta between two vectors in 3D space, once the vectors move
        >180 degrees apart, the angle becomes the shortest path. To correct this, we calculate
        theta within the plane by checking if it faces the same direction as a vector y, and correcting
        otherwise.

    Args:
        rho (Vector) : the vector compared to the x-axis to form theta
        y_axis (Vector) : directional vector to define a plane with x-axis
        theta (float) : the calculated angle in 3D space to be corrected
    Returns:
        theta (float) : theta as calculated on the plane
    '''
    proj_rho_on_y = rho.calculate_projection(y_axis)
    opposite_direction = (proj_rho_on_y.x/y_axis.x < 0)
    if opposite_direction:
        theta = 360-theta
    return theta
    
def calculate_bottaro_values_for_frame(perspective_base_coords : Base, viewed_midpoint : tuple) -> list:
    '''Calculates the r, rho, and theta values as expressed in the Bottaro paper

    Calculates the r, rho, and theta values between two nucleotides in a single frame
    as presented in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)
    
    Args: 
        perspective_base_coords (Base) : list of the x,y,z coords of C2, C4, C6 and their midpoint for 
            the perspective residue in a single frame
        viewed_midpoint (tuple) : x,y,z position of the midpoint of the viewed nucleotide at the same frame.
    Returns: 
        values (list) : a list containing 3 floats from Bottaro; structure: [r_dist, rho_dist, theta]
    '''
    r_vector = viewed_midpoint-perspective_base_coords.midpoint_coords
    r_magnitude = r_vector.magnitude()
    
    # 2 Vectors define the plane of the perspective nucleotide
    x_axis = perspective_base_coords.c2_coords-perspective_base_coords.midpoint_coords
    vector_midpoint_to_C4 = perspective_base_coords.c4_coords-perspective_base_coords.midpoint_coords

    normal_vector_to_plane = x_axis.calculate_cross_product(vector_midpoint_to_C4)
    
    # y-axis inside plane to get correct theta value in 3D
    y_axis = x_axis.calculate_cross_product(normal_vector_to_plane)

    proj_r_on_normal = r_vector.calculate_projection(normal_vector_to_plane)

    rho = r_vector-proj_r_on_normal
    rho_dist = rho.magnitude()
    
    x_axis_dot_rho = x_axis.x*rho.x+x_axis.y*rho.y+x_axis.z*rho.z
    denominator = x_axis.magnitude()*rho_dist
    cos_theta = x_axis_dot_rho/denominator
    theta = math.degrees(math.acos(cos_theta))
    corrected_theta = correct_theta_sign(rho, y_axis, theta)
    
    values = [r_magnitude,rho_dist,corrected_theta]
    return values

def write_bottaro_to_csv(pdb_filename : str, output_csv_name : str, 
                         res1_atom_names : tuple = ("C2", "C4","C6"), 
                         res2_atom_names : tuple = ("C2", "C4","C6")) -> None:
    '''Write the Bottaro r, rho, and theta values from a trajectory pdb to a CSV

    Calculates the r, rho, and theta values as described in Bottaro et al. from a
    perspective nucleotide to a viewed nucleotide per frame. Writes the results to a CSV file.

    Args:
        pdb_filename (str) : name of pdb containing information for ONLY two residues (perspective and viewed
            nucleotide) at each frame.
        output_csv_name (str) : filename of CSV file to write to
        res1_atom_names (tuple) : tuple of the atom names (eg. "C2", "C4", "C6") to use from
            residue 1 to find center of geometry for perspective nucleotide
        res2_atom_names (tuple) : tuple of the atom names (eg. "C2", "C4", "C6") to use from
            residue 2 to find center of geometry for viewed nucleotide
    '''
    res1_atom1,res1_atom2,res1_atom3 = res1_atom_names
    res2_atom1,res2_atom2,res2_atom3 = res2_atom_names
    
    pdb = md.load(pdb_filename)
    topology = pdb.topology
    number_of_frames = pdb.n_frames
    perspective_residue_num = [residue for residue in topology.residues][0].resSeq
    viewed_residue_num = [residue for residue in topology.residues][1].resSeq

    residue1_C2_list = collect_atom_locations_by_frame(pdb, perspective_residue_num, res1_atom1)
    residue1_C4_list = collect_atom_locations_by_frame(pdb, perspective_residue_num, res1_atom2)
    residue1_C6_list = collect_atom_locations_by_frame(pdb, perspective_residue_num, res1_atom3)
    residue2_C2_list = collect_atom_locations_by_frame(pdb, viewed_residue_num, res2_atom1)
    residue2_C4_list = collect_atom_locations_by_frame(pdb, viewed_residue_num, res2_atom2)
    residue2_C6_list = collect_atom_locations_by_frame(pdb, viewed_residue_num, res2_atom3)
    
    residue1_midpoint_list = [calc_center_3pts(Vector(*residue1_C2_list[i]),
                                               Vector(*residue1_C4_list[i]),
                                               Vector(*residue1_C6_list[i])) for i in range(0, number_of_frames)]
    residue2_midpoint_list = [calc_center_3pts(Vector(*residue2_C2_list[i]),
                                               Vector(*residue2_C4_list[i]),
                                               Vector(*residue2_C6_list[i])) for i in range(0, number_of_frames)]
    
    fields = ['frame','r_dist','rho_dist', 'theta']
    rows=[]
    for i in range(0,number_of_frames):
        residue1_base = create_base_from_coords_list(i, residue1_C2_list, residue1_C4_list, residue1_C6_list, residue1_midpoint_list)
        frame_values = calculate_bottaro_values_for_frame(residue1_base,residue2_midpoint_list[i])
        row = [i]+frame_values
        rows.append(row)
    
    filename = output_csv_name
    with open(filename, 'w') as csvfile:
        csvwriter = csv.writer(csvfile) 
        csvwriter.writerow(fields) 
        csvwriter.writerows(rows)

if __name__ == "__main__":
    ########JOB VARIABLES#######
    #load pdb file
    pdb_filename = 'new_short_wCCC_+1GCU_GC.pdb'

    #name of output file
    output_csv_name = "wCCC_+1GCU_GC_plot.csv"
    ############################

    ########OPTIONAL VARS#######
    perspective_atom1_name = "C2"
    perspective_atom2_name = "C4"
    perspective_atom3_name = "C6"
    viewed_atom1_name = "C2"
    viewed_atom2_name = "C4"
    viewed_atom3_name = "C6"
    ############################

    write_bottaro_to_csv(pdb_filename, 
                         output_csv_name, 
                         (perspective_atom1_name, perspective_atom2_name, perspective_atom3_name), 
                         (viewed_atom1_name,viewed_atom2_name,viewed_atom3_name))