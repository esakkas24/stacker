"""
Visualize the SSFs, PSFs, and other analyses

This module includes the functions used to visualize plots,
including SSFs and PSFs. The data inputs to these plot functions
are provided by the other modules.
"""

import os
import functools
import numpy as np
from numpy import typing
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from seaborn import kdeplot
from .file_manipulation import SmartIndexingAction

def create_parent_directories(outfile_prefix : str) -> None:
    '''
    Creates necessary parent directories to write an outfile given a prefix
    
    Parameters
    ----------
    outfile_prefix : str
        The filepath of the output file, including the path where the file will be saved.

    Examples
    --------
    >>> create_parent_directories('/path/to/output/file.txt')
    # This will create the directories '/path/to/output' if they do not exist.
    >>> create_parent_directories('output/file.txt')
    # This will create the directory 'output' if it does not exist.
    '''
    dir_name = os.path.dirname(outfile_prefix)
    if dir_name == '': dir_name = '.'
    os.makedirs(dir_name, exist_ok=True)

def create_axis_labels(res_indices: typing.ArrayLike, tick_distance: int = 10) -> list:
    """
    Designates the axis labels to use in the SSF plot.

    Helper function for visualizing SSFs.
    Returns the x-axis tick positions and labels to use on the ticks based on the 
    residues used in a specific SSF analysis. Meant to be used when many 
    disjoint sets of residue indices are used. Ticks will be present every `tick_distance` 
    residues in a collection of adjacent residues, and a tick will exist at both
    ends of any consecutive residue sequence.

    Parameters
    ----------
    res_indices : list
        The list of residue indices used in the pairwise analysis.
        Parameter `residue_desired` passed to `filter_traj()`
    tick_distance : int, default = 10
        Distance between ticks in blocks of consecutive residues.

    Returns
    -------
    tick_locations : list
        List of tick positions (0-based) to place labels on in the axes.
    tick_labels : list
        List of labels to place at the adjacent tick locations.

    See Also
    --------
    filter_traj : Filters an input trajectory to desired residues

    Examples
    --------
    Residues 0-12,98-100 were used. The SSF will label 0,10,12,98,100,
    provided in the second returned list. The first returned list gives 
    the positions on the axes to place each label.

    >>> create_axis_labels([0,1,2,3,4,5,6,7,8,9,10,11,12,98,99,100])
    [0, 10, 12, 13, 15], [0, 10, 12, 98, 100]

    >>> create_axis_labels([94,95,96,97,98,99,100,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428])
    [0, 6, 7, 17, 27], [94, 100, 408, 418, 428]
    """
    n_residues = len(res_indices)

    if n_residues < 1: raise NoResidues("pairwise analysis must include at least one residue")

    tick_locations = [0]
    tick_labels = [res_indices[0]]

    res_sequence_length = 1

    for i in range(1, n_residues):
        if res_indices[i] == res_indices[i-1] + 1:
            res_sequence_length += 1

        if res_indices[i] > res_indices[i-1] + 1:
            tick_locations += [i-1, i]
            tick_labels += [res_indices[i-1], res_indices[i]]
            res_sequence_length = 1
        elif res_sequence_length == tick_distance+1:
            tick_locations += [i]
            tick_labels += [res_indices[i]]
            res_sequence_length = 1
    
    if n_residues-1 not in tick_locations:
        tick_locations += [n_residues-1]
        tick_labels += [res_indices[n_residues-1]]
    
    return tick_locations, tick_labels
        
def display_arrays_as_video(numpy_arrays : list | typing.ArrayLike, res_indices: typing.ArrayLike | str, 
                            seconds_per_frame: int = 10, tick_distance: int = 10,
                            outfile_prefix: str = '', scale_limits: tuple = (0, 7), outfile: str = '',
                            scale_style: str = 'bellcurve', xy_line: bool = True) -> None:
    """
    display_arrays_as_video(
        ssfs,
        res_indices,
        seconds_per_frame = 10,
        tick_distance = 10,
        outfile_prefix = '',
        scale_limits = (0,7),
        outfile = '',
        scale_style = 'bellcurve',
        xy_line = True
    )

    Displays SSF data to output or writes SSF as a PNG

    Visualizes the data for an SSF for a trajectory or a single frame.
    Takes an SSF array outputted from `get_residue_distance_for_frame`,
    `get_residue_distance_for_trajectory`,
    or `system_stacking_fingerprints` and treats them as frames 
    in a video, filling in a grid at position i, j by the value 
    at i, j in the array.

    Parameters
    ----------
    ssfs : array_like
        List or array of 2D NumPy arrays representing SSFs, 
        output of ``system_stacking_fingerprints``
    res_indices : list or str
        The list of residue indices used in the pairwise analysis.
        Parameter `residue_desired` passed to `filter_traj()`
        Accepts smart-indexed str representing a list of residues (e.g '1-5,6,39-48')
    seconds_per_frame : int, default = 10
        Number of seconds to display each matrix for.
    tick_distance : int, default = 10
        Distance between ticks in blocks of consecutive residues.
    outfile : str
        Image output filepath to write a single SSF to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    outfile_prefix : str
        Prefix for image filepath to write multiple frames to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    scale_limits : tuple, default = (0, 7)
        Limits of the color scale.
    scale_style : {'bellcurve', 'gradient'}, default = 'bellcurve'
        Style of color scale. 
    xy_line : bool, default = True
        Draw x = y line to separate matrix halves.

    Returns
    -------
    None
        Displays video of NumPy arrays.

    See Also
    --------
    create_axis_labels : Designates the axis labels to use in the SSF plot.
    get_residue_distance_for_frame : Calculates System Stacking Fingerprint (SSF) between all residues in a given frame.
    get_residue_distance_for_trajectory : get SSF data for all frames of a trajectory
    system_stacking_fingerprints : Alias for this `get_residue_distance_for_trajectory`
    display_ssfs : Alias for this function.
    
    Examples
    --------
    >>> import stacker as st
    >>> import mdtraj as md
    >>> trajectory_file = 'stacker/testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    >>> topology_file = 'stacker/testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    >>> trj = md.load(trajectory_file, top = topology_file)
    >>> residue_selection_query = 'resi 90 to 215'
    >>> frames_to_include = [1,2,3,4,5]
    >>> trj_sub = trj.atom_slice(trj.top.select(residue_selection_query))
    >>> resSeqs = [res.resSeq for res in trj_sub.topology.residues]
    >>> frames = st.get_residue_distance_for_trajectory(trj_sub, frames_to_include, threads = 5)
    Frame 2 done.
    Frame 3 done.
    Frame 1 done.
    Frame 5 done.
    Frame 4 done.
    >>> st.display_arrays_as_video([st.get_frame_average(frames)], resSeqs, seconds_per_frame=10)
    # Displays SSF for each frame of this trajectory to standard output

    """
    orange_colormap = mpl.colormaps['Oranges_r'].resampled(100)

    if scale_style == 'gradient':
        newcolors = np.vstack((orange_colormap(np.linspace(1, 1, 1)), orange_colormap(np.linspace(0, 0, 128)),
                        orange_colormap(np.linspace(0, 1, 128))))
        newcmp = mpl.colors.ListedColormap(newcolors, name='OrangeBellcurve')
    elif scale_style == 'bellcurve':
        newcolors = np.vstack((orange_colormap(np.linspace(1, 0, 128)), orange_colormap(np.linspace(0, 1, 128))))
        newcmp = mpl.colors.ListedColormap(newcolors, name='OrangeBellcurve')
    
    fig, ax = plt.subplots(figsize=(8,8))
    plt.ion()
    frame_num = 1
    for hist in numpy_arrays:
        ax.clear()
        vmin, vmax = scale_limits
        neg = ax.imshow(hist, cmap = newcmp, vmin=vmin, vmax=vmax, interpolation = 'nearest')
        ax.set_title('Distance Between Residues Center of Geometries')
        ax.set_title('Distance Between Residues Center of Geometries')
        ax.set_xlabel('Residue Index')  
        ax.xaxis.set_label_position('top')  
        ax.set_ylabel('Residue Index')  
        colorbar = fig.colorbar(neg, ax=ax, location='right', anchor=(0, 0.3), shrink=0.7)
        colorbar.ax.set_title('Center of\nGeometry\nDist. (Å)')

        res_indices = SmartIndexingAction.parse_smart_index(res_indices)
        res_indices = list(res_indices)
        ticks, labels = create_axis_labels(res_indices, tick_distance)
        plt.xticks(ticks, labels, rotation = 'vertical')
        plt.yticks(ticks, labels)
        ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

        # seperate ticks by region
        last_res = 1
        last_label = 0
        long_tick_region = False
        for res_i, label_i, tick_object, y_object in zip(ticks, labels, ax.xaxis.get_major_ticks(), ax.yaxis.get_major_ticks()): 
            if res_i == last_res+1 and not label_i == last_label+1:
                long_tick_region = not long_tick_region

            if long_tick_region:
                tick_object.set_pad(22.)
                tick_object.tick2line.set_markersize(22.)
                y_object.set_pad(22.)
                y_object.tick1line.set_markersize(22.)
            else:
                tick_object.set_pad(2.)
                tick_object.tick2line.set_markersize(3.)
                y_object.set_pad(2.)
                y_object.tick1line.set_markersize(3.)
            
            last_res = res_i
            last_label = label_i

        if xy_line:
            lims = [np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
                    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
                    ]
            ax.plot(lims, lims, 'k-', zorder = 0, alpha = 0.75, linewidth = 0.5) 

        plt.pause(seconds_per_frame)
        if outfile_prefix:
            plt.savefig(outfile_prefix + "frame" + str(frame_num) + ".png")
        elif outfile:
            plt.savefig(outfile)
        colorbar.remove()
        frame_num+=1

@functools.wraps(display_arrays_as_video)
def display_ssfs(*args, **kwargs):
    return display_arrays_as_video(*args, **kwargs)

display_ssfs.__doc__ = f"""
Alias for `display_arrays_as_video()`.

{display_arrays_as_video.__doc__}
"""

def set_polar_grid() -> mpl.projections.polar.PolarAxes:
    """
    Set up axes for PSF.

    Creates polar plot background for two-residue movement comparison
    with theta 0 to 360, a radial maximum of 15 Angstroms, and a visualization 
    of the perspective residue at the center.

    Parameters
    ----------
    None

    Returns
    -------
    ax : matplotlib.projections.polar.PolarAxes
        Axis object for the created polar plot.

    """
    fig = plt.figure()
    ax = fig.add_subplot(polar=True)

    ax.set_xticks(np.pi/180. * np.linspace(0,  360, 3, endpoint=False))
    ax.set_xticklabels([r"$\theta=0^\circ$",r"$\theta=120^\circ$",r"$\theta=240^\circ$"])
    ax.tick_params(pad = -10)

    ax.set_rlim(0,10)
    ax.set_rticks(np.linspace(0,  10, 3, endpoint=True))
    ax.set_rlabel_position(180)
    plt.text(x=np.radians(178), y=12, s=r"$\rho\text{ }(\AA)$", ha="center",va='center',fontsize=11)
    plt.text(x=0, y=2, s="C2", ha="center",va='center',fontsize=7.5)
    plt.text(x=np.radians(240), y=2, s="C4", ha="center",va='center',fontsize=7.5)
    plt.text(x=np.radians(120), y=1.8, s="C6", ha="center",va='center',fontsize=7.5)

    ax.grid(color='gray', linestyle='--', linewidth=0.5)
    return ax

def visualize_two_residue_movement_scatterplot(csv: str, plot_outfile: str = '', frame_list: set = {}) -> None:
    """
    Creates scatterplot of two-residue movement relative to each other.

    Takes the data created in residue_movement and visualizes it as a polar coordinate
    scatterplot similar to the Figure D link in Proposal Feature 4.

    Parameters
    ----------
    csv : str
        Filepath to csv file containing data on the movement
        of two residues relative to each other (r, rho, and theta values). Created
        in residue_movement.
    plot_outfile : str
        Filepath of the image file to write to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    frame_list : set, default = {}
        Set of frames to use in csv, if empty use all frames.

    Returns
    -------
    None

    See Also
    --------
    write_bottaro_to_csv : Creates CSV file that is inputted here

    Examples
    --------
    >>> import stacker as st
    >>> trajectory_file = 'testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    >>> topology_file = 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    >>> pdb_filename = 'testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb'
    >>> output_csv_name = "testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv"
    >>> perspective_residue = 426 # 1-indexed
    >>> viewed_residue = 427 # 1-indexed
    >>> st.filter_traj_to_pdb(trj_file=trajectory_file, top_file=topology_file, output_pdb_filename=pdb_filename,
    ...                        residues_desired={perspective_residue,viewed_residue}, atomnames_desired={"C2", "C4", "C6"})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    WARNING: Output file atom, residue, and chain indices are zero-indexed
    Filtered trajectory written to:  testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb
    >>> st.write_bottaro_to_csv(pdb_filename, output_csv_name, perspective_residue_num=perspective_residue, viewed_residue_num=viewed_residue)
    Output values written to testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv
    >>> st.visualize_two_residue_movement_scatterplot('testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv', 
    ...                                                 plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_10frames_scatter.png')

    """
    bottaro_values = pd.read_csv(csv, sep=',')

    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values = bottaro_values['theta']
    theta_values_rad = np.radians(theta_values)

    rho_values = bottaro_values['rho_dist']

    ax = set_polar_grid()
    ax.scatter(theta_values_rad, rho_values, color = 'purple', s=1, alpha = 0.5)

    # Draw nucleotide ring in polar plot
    r_for_ring = np.ones(7)*1.3
    theta_for_ring = np.linspace(0, 2 * np.pi, 7)   
    ax.fill(theta_for_ring,r_for_ring, color = 'black', fill=False)

    r_for_purine = np.array([1.3, 1.3, 2.58575693, 3.075, 2.58575693, 3.49, 2.58575693, 1.3])
    theta_for_purine = np.array([4*np.pi/3, np.pi, -3.04534743, -2.61795,-2.19064032, -1.88, -2.19064032, 4*np.pi/3])  
    ax.fill(theta_for_purine, r_for_purine, color = 'black', fill=False, alpha = 0.5)

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()


def visualize_two_residue_movement_heatmap(csv: str, plot_outfile: str = '', frame_list: set = {}) -> None:
    """
    Creates heatmap of two-residue movement relative to each other.

    2D shaded contour plot of the density of points in the 
    visualize_two_residue_movement_scatterplot() scatterplot.

    Parameters
    ----------
    csv : str
        Filepath to csv file containing data on the movement
        of two residues relative to each other (r, rho, and theta values). Created
        in residue_movement.
    plot_outfile : str
        Filepath of the image file to write to. Format inferred from file extension.
        png, pdf, ps, eps, and svg supported.
    frame_list : set, default = {}
        Set of frames to use in csv, if empty use all frames.

    Returns
    -------
    None

    See Also
    --------
    write_bottaro_to_csv : Creates CSV file that is inputted here
    
    Examples
    --------
    >>> import stacker as st
    >>> trajectory_file = 'testing/first10_5JUP_N2_tUAG_aCUA_+1GCU_nowat.mdcrd'
    >>> topology_file = 'testing/5JUP_N2_tUAG_aCUA_+1GCU_nowat.prmtop'
    >>> pdb_filename = 'testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb'
    >>> output_csv_name = "testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv"
    >>> perspective_residue = 426 # 1-indexed
    >>> viewed_residue = 427 # 1-indexed
    >>> st.filter_traj_to_pdb(trj_file=trajectory_file, top_file=topology_file, output_pdb_filename=pdb_filename,
    ...                        residues_desired={perspective_residue,viewed_residue}, atomnames_desired={"C2", "C4", "C6"})
    WARNING: Residue Indices are expected to be 1-indexed
    Reading trajectory...
    Reading topology...
    Filtering trajectory...
    WARNING: Output filtered traj atom, residue, and chain indices are zero-indexed
    WARNING: Output file atom, residue, and chain indices are zero-indexed
    Filtered trajectory written to:  testing/script_tests/residue_movement/5JUP_N2_tUAG_aCUA_+1GCU_nowat_mdcrd.pdb
    >>> st.write_bottaro_to_csv(pdb_filename, output_csv_name, perspective_residue_num=perspective_residue, viewed_residue_num=viewed_residue)
    Output values written to testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv
    >>> st.visualize_two_residue_movement_heatmap('testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv', 
    ...                                                 plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_10frames_heat.png')

    """
    bottaro_values = pd.read_csv(csv, sep=',')

    if frame_list:
        bottaro_values = bottaro_values[bottaro_values['frame'].isin(frame_list)]

    theta_values = bottaro_values['theta']
    theta_values_rad = np.radians(theta_values)

    # convert rho values from nm to Angstroms
    rho_values = bottaro_values['rho_dist']

    ax = set_polar_grid()
    ax = kdeplot(x=theta_values_rad, y=rho_values, fill=True, levels = [0.1*i for i in range(1,11)], cbar = True, cmap = 'gist_earth_r')
    plt.xlabel('')
    plt.ylabel('')

    cbar = ax.collections[-1].colorbar
    levels = [0.1*i for i in range(1,11)]
    formatted_labels = ['{:.1f}'.format(level) for level in levels]
    cbar.set_ticklabels(formatted_labels)

    # Draw nucleotide ring in polar plot
    r_for_ring = np.ones(7)*1.3
    theta_for_ring = np.linspace(0, 2 * np.pi, 7)   
    ax.fill(theta_for_ring,r_for_ring, color = 'black', fill=False)

    r_for_purine = np.array([1.3, 1.3, 2.58575693, 3.075, 2.58575693, 3.49, 2.58575693, 1.3])
    theta_for_purine = np.array([4*np.pi/3, np.pi, -3.04534743, -2.61795,-2.19064032, -1.88, -2.19064032, 4*np.pi/3])  
    ax.fill(theta_for_purine, r_for_purine, color = 'black', fill=False, alpha = 0.5)

    if plot_outfile:
        plt.savefig(plot_outfile)
    else:
        plt.show()

@functools.wraps(visualize_two_residue_movement_scatterplot)
def display_psf_scatter(*args, **kwargs):
    return visualize_two_residue_movement_scatterplot(*args, **kwargs)

display_psf_scatter.__doc__ = f"""
Alias for `visualize_two_residue_movement_scatterplot()`.

{visualize_two_residue_movement_scatterplot.__doc__}
"""

@functools.wraps(visualize_two_residue_movement_heatmap)
def display_psf_heatmap(*args, **kwargs):
    return visualize_two_residue_movement_heatmap(*args, **kwargs)

display_psf_heatmap.__doc__ = f"""
Alias for `visualize_two_residue_movement_heatmap()`.

{visualize_two_residue_movement_heatmap.__doc__}
"""

class NoResidues(Exception):
    """Raised if user tries to make SSF with a trajectory of <1 residue"""
    pass

if __name__ == '__main__':
    # 10 frame test
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot.csv')

    # Multi-frame test
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv')

    # Multi-frame heatmap test
    visualize_two_residue_movement_heatmap('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv')

    # Write to outfile tests
    output = "stacker/testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_scatter.png"
    create_parent_directories(output)
    visualize_two_residue_movement_scatterplot('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv', plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_scatter.png')
    visualize_two_residue_movement_heatmap('stacker/testing/script_tests/residue_movement/tUAG_aCUA_+1GCU_GC_plot_3200frames.csv', plot_outfile='testing/script_tests/visualization/tUAG_aCUA_+1GCU_GC_plot_3200frames_heatmap.png')

