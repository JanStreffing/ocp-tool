********
OCP-Tool
********

OCP-Tool, the OpenIFS Coupling Preparation tool, makes the OpenIFS side of a coupled setup agree with an ocean grid. You give it an OpenIFS truncation and a FESOM2 mesh (or a NEMO grid, or no ocean at all for AMIP), and it writes the files that have to share one land-sea mask: the modified OpenIFS initial files, the OASIS3-MCT ``grids.nc``, ``masks.nc`` and ``areas.nc``, the runoff mapper input and the LPJ-GUESS soil type file.

You need it whenever the atmosphere grid or the ocean mesh changes. A new FESOM2 mesh, a new truncation, a mesh with or without ice shelf cavities and a paleo time slice are all the same job from the tool's point of view. A model whose atmosphere mask and coupler mask disagree runs, and exchanges fluxes across points that do not exist on the other side, so do not reuse files from a different grid pair to save the hour.

What it does not do is generate the OASIS ``rmp_`` remapping weights for a production grid pair. It tries (step 4b), and in the normal case reports that it could not. Weights come from the offline weight generation described in :doc:`awiesm3:how_to/coupling_oasis`, driven by the ``tools/submit_weightgen_*.sh`` scripts in this repository.

The tool is part of the AWI-CM3 and AWI-ESM3 workflow. The model side of things (where the files go in a runscript, how to build a paleo setup, what to do with the weights) is documented at https://awi-esm3.readthedocs.io, and this site links there rather than repeating it.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   api

Indices
=======

* :ref:`genindex`
* :ref:`modindex`
