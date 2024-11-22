
***************
Getting Started
***************

.. currentmodule:: stacker

Welcome to getting started with StACKER!

StACKER (**St**\ acking **A**\ nalysis and **C**\ onformational **K**\ inetics for **E**\ xamining 
**R**\ esidues) is an open source Python library for nalyzing pi-stacking interactions 
between residues in Molecular Dynamics (MD) trajectories.

Installing StACKER
------------------

StACKER can be installed through GitHub or through PyPi:

In the command line, run::

    pip install pistacker

This will install StACKER, activate the command line option ``stacker``, 
and install all the necessary dependencies.

Alternatively, you can clone the GitHub repository to your local machine.

In the command line, run::

    git clone https://github.com/esakkas24/stacker.git

Descend into the directory and download the neccessary dependencies::

    cd stacker
    pip install -r requirements.txt
    pip install setuptools
    python setup.py install

This will install StACKER, activate the command line option `stacker`, 
and install all the necessary dependencies.

Read more about StACKER's :ref:`Command Line Options <command_line_options>`.

StACKER in the Command Line
---------------------------

StACKER's functions are available through a Python import and through
the Command Line. After installing StACKER, running ``stacker`` in the 
command line will show the Command Line Options for StACKER::

    [user]$ stacker
    usage: stacker -s ROUTINE [-h]

    Wrapper to run stacker subroutines using the -s flag.
    More info on each routine given by `stacker -s ROUTINE -h` or `python stacker.py -s ROUTINE --help`

    options:
    -s ROUTINE, --script ROUTINE
                Name of command to use. Options for ROUTINE:

                filter_traj:
                        filters trajectory and topology files to desired residue numbers and atom names
                bottaro OR pairwise OR psf:
                        Create a Pairwise Stacking Fingerprint (PSF): polar plots like those in Figure 1 of Bottaro et. al (https://doi.org/10.1093/nar/gku972)
                res_distance:
                        Get the distance between two residues in a given frame
                system OR ssf:
                        Create a System Stacking Fingerprint (SSF) with pairwise distances for each residue
                stack_events:
                        Get list of residues with most stacking events (distance closest to 3.5Å)
                compare:
                        Get the most changed stacking events between two fingerprints using the outputs of python stacker.py -s stack_events

    -h, --help            show this help message and exit


Importing StACKER to Python
---------------------------

After installing StACKER, you can import it into Python
code like::

    import stacker as st

Although StACKER's functions are further subsetted into :ref:`modules <index>`,
all functions are accessible through ``st.________``::

    st.<Tab>
    st.AtomEmpty(                                   st.create_axis_labels(                          st.os                                           
    st.Base(                                        st.create_base_from_coords_list(                st.pairwise_distance                            
    st.FrameEmpty(                                  st.create_parent_directories(                   st.pd                                           
    st.InvalidRoutine(                              st.csv                                          st.plt                                          
    st.MultiFrameTraj(                              st.display_arrays_as_video(                     st.print_function                               
    st.NoResidues(                                  st.enable_printing()                            st.random                                       
    st.ResEmpty(                                    st.file_convert(                                st.res_distance_routine()                       
    st.SmartIndexingAction(                         st.file_manipulation                            st.residue_movement                             
    st.Vector(                                      st.filter_traj(                                 st.run_python_command()                         
    st.argparse                                     st.filter_traj_routine()                        st.set_polar_grid()                             
    st.block_printing()                             st.filter_traj_to_pdb(                          st.stack_events_routine()                       
    st.bottaro_routine()                            st.get_frame_average(                           st.sys                                          
    st.calc_center_3pts(                            st.get_residue_distance_for_frame(              st.system_routine()                             
    st.calculate_bottaro_values_for_frame(          st.get_residue_distance_for_trajectory(         st.typing                                       
    st.calculate_residue_distance(                  st.get_top_stacking(                            st.vector                                       
    st.collect_atom_locations_by_frame(             st.increment_residue(                           st.visualization                                
    st.combine_frames(                              st.kdeplot(                                     st.visualize_two_residue_movement_heatmap(      
    st.compare_routine()                            st.math                                         st.visualize_two_residue_movement_scatterplot(  
    st.concurrent                                   st.md                                           st.write_bottaro_to_csv(                        
    st.convert_to_python_command()                  st.mpl                                                                                          
    3 more... 

Reading the example code
------------------------

Throughout the NumPy documentation, you will find blocks that look like::

  >>> a = np.array([[1, 2, 3],
  ...               [4, 5, 6]])
  >>> a.shape
  (2, 3)

Text preceded by ``>>>`` or ``...`` is **input**, the code that you would
enter in a script or at a Python prompt. Everything else is **output**, the
results of running your code. Note that ``>>>`` and ``...`` are not part of the
code and may cause an error if entered at a Python prompt.

Why use NumPy?
--------------

Python lists are excellent, general-purpose containers. They can be
"heterogeneous", meaning that they can contain elements of a variety of types,
and they are quite fast when used to perform individual operations on a handful
of elements.

Depending on the characteristics of the data and the types of operations that
need to be performed, other containers may be more appropriate; by exploiting
these characteristics, we can improve speed, reduce memory consumption, and
offer a high-level syntax for performing a variety of common processing tasks.
NumPy shines when there are large quantities of "homogeneous" (same-type) data
to be processed on the CPU.
