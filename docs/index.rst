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
        :img-top: _static/images/command_line_options.png
        :text-align: center

        Command Line Options
        ^^^

        Learn how to use the command line options for StACKER.

        +++

        .. button-ref:: stacker/command_line_options
            :expand:
            :color: secondary
            :click-parent:

            To the command line options

    .. grid-item-card::
        :img-top: _static/images/load_filter_trajectories.png
        :text-align: center

        Load and Filter Trajectories
        ^^^

        Instructions on loading and filtering trajectories.

        +++

        .. button-ref:: stacker/load_filter_trajectories
            :expand:
            :color: secondary
            :click-parent:

            To the loading and filtering trajectories

    .. grid-item-card::
        :img-top: _static/images/create_ssf.png
        :text-align: center

        Create System Stacking Fingerprints (SSFs)
        ^^^

        Steps to create System Stacking Fingerprints.

        +++

        .. button-ref:: stacker/create_ssf
            :expand:
            :color: secondary
            :click-parent:

            To the SSF creation guide

    .. grid-item-card::
        :img-top: _static/images/create_psf.svg
        :text-align: center

        Create Pairwise Stacking Fingerprints (PSFs)
        ^^^

        Steps to create Pairwise Stacking Fingerprints.

        +++

        .. button-ref:: stacker/create_psf
            :expand:
            :color: secondary
            :click-parent:

            To the PSF creation guide

    .. grid-item-card::
        :img-top: _static/images/visualize_data.svg
        :text-align: center

        Visualize Data
        ^^^

        Learn how to visualize your data with StACKER.

        +++

        .. button-ref:: stacker/visualize_data
            :expand:
            :color: secondary
            :click-parent:

            To the data visualization guide

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   Command Line Options <stacker>
   Load and Filter Trajectories <stacker.file_manipulation>
   Create System Stacking Fingerprints (SSFs) <stacker.pairwise_distance>
   Create Pairwise Stacking Fingerprints (PSFs) <stacker.residue_movement>
   Visualize Data <stacker.visualization>


.. autosummary::
      :toctree: _autosummary
      :template: custom-module-template.rst
      :recursive:

      stacker
      stacker.file_manipulation
      stacker.kmeans
      stacker.pairwise_distance
      stacker.residue_movement
      stacker.vector
      stacker.visualization

