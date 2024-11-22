.. _stacker_docs_mainpage:

##########################
StACKER documentation
##########################


StACKER manipulates the outputs of an MD simulation and analyzes 
the pi-stacking interactions. Creates a "Pi-Stacking Fingerprint" 
for a structure at each frame. Presents Pi-stacking interactions between 
two residues through analysis of their relevant movement.

.. grid:: 1 1 2 2
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: _static/index_images/getting_started.svg
        :text-align: center

        Getting Started
        ^^^

        Learn how to install and use
        the basic StACKER functions

        +++

        .. button-ref:: getting_started
            :expand:
            :color: secondary
            :click-parent:

            To the installation instructions

    .. grid-item-card::
        :img-top: _static/index_images/python_logo.webp
        :text-align: center

        Python Modules
        ^^^

        See detailed instructions on StACKER's 
        functions and modules when using Python

        +++

        .. button-ref:: modules
            :expand:
            :color: secondary
            :click-parent:

            To the Python user guide

    .. grid-item-card::
        :img-top: _static/index_images/system_stacking_fingerprint.png
        :text-align: center

        Create System Stacking Fingerprints (SSFs)
        ^^^

        Pipeline to create System Stacking Fingerprints (SSFs)
        using the Python and Command Line Options

        +++

        .. button-ref:: stacker.pairwise_distance
            :expand:
            :color: secondary
            :click-parent:

            To the SSF creation guide

    .. grid-item-card::
        :img-top: _static/index_images/system_stacking_fingerprint.png
        :text-align: center

        Create Pairwise Stacking Fingerprints (PSFs)
        ^^^

        Pipeline to create Pairwise Stacking Fingerprints (PSFs)
        using the Python and Command Line Options

        +++

        .. button-ref:: stacker/create_psf
            :expand:
            :color: secondary
            :click-parent:

            To the PSF creation guide

    .. grid-item-card::
        :img-top: _static/index_images/visualize_data.svg
        :text-align: center

        Command Line Options
        ^^^

        See detailed instructions on StACKER's 
        functions and modules when using 
        the Command Line
        +++

        .. button-ref:: stacker/visualize_data
            :expand:
            :color: secondary
            :click-parent:

            To the Command Line user guide

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Getting Started <getting_started>
   Command Line Options <command_line_options>
   Python Options <python_options>
   Load and Filter Trajectories <stacker.file_manipulation>
   Create System Stacking Fingerprints (SSFs) <stacker.pairwise_distance>
   Create Pairwise Stacking Fingerprints (PSFs) <stacker.residue_movement>
   Visualize Data <stacker.visualization>

