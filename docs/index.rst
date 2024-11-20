.. StACKER documentation master file, created by
   sphinx-quickstart on Wed Nov 20 10:14:40 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

StACKER documentation
=====================

Add your content using ``reStructuredText`` syntax. See the
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
documentation for details.

.. autosummary::
      :toctree: _autosummary
      :template: custom-module-template.rst
      :recursive:

      stacker
      stacker.file_manipulation
      stacker.kmeans
      stacker.pairwise_distance
      stacker.residue_movment
      stacker.vector
      stacker.visualization

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   File Manipulation <stacker.file_manipulation>
   KMeans <stacker.kmeans>
   Pairwise Distance <stacker.pairwise_distance>
   Residue Movement <stacker.residue_movement>
   Vector <stacker.vector>
   Visualization <stacker.visualization>