************
Installation
************

The tool is run from a clone of the repository rather than installed as a package, because it reads its default inputs from ``input/`` and writes to ``output/`` relative to the checkout.

.. code-block:: bash

  git clone https://github.com/JanStreffing/ocp-tool.git
  cd ocp-tool

Besides the Python packages it calls two command line tools: ``cdo`` to describe the input grid, and ``ncrename`` from NCO to write the LPJ-GUESS soil file. Both have to be on ``PATH`` when the tool runs, not only when you install it.

On Levante
==========

A working environment already exists, and using it is quicker than building your own:

.. code-block:: bash

  source /work/ab0246/a270092/software/miniforge3/etc/profile.d/conda.sh
  conda activate ocp-tool2
  module load nco/5.0.6-gcc-11.2.0

Use this miniforge and not the system ``mambaforge`` module. The system Python fails to load the compiled libraries in the environment with a ``GLIBCXX`` version error.

Do not leave out the NCO module. The environment ships ``cdo`` but not ``ncrename``, and the step that needs it (step 12) comes after all the others in a non-paleo run, so you only find out once the rest has finished. See :ref:`trap_ncrename`.

Running with a bare ``python3`` rather than the environment's one tends to pick up a newer matplotlib from ``~/.local``, which breaks the plotting step. That failure looks like a tool bug and is not.

The first import of the heavy packages can sit in Lustre I/O wait for minutes. If ``python -c pass`` returns in a couple of seconds, the environment is fine and the tool is just loading.

Elsewhere
=========

Create an environment from ``environment.yaml`` and add the packages it does not yet list:

.. code-block:: bash

  conda env create -f environment.yaml
  conda activate ocp-tool2
  conda install -c conda-forge cartopy shapely xarray cmocean cdo nco

``environment.yaml`` installs pyfesom2 from its default branch with pip. The pixi configuration uses the ``modernize-packaging`` branch instead, for Python 3.11 and later; try that one if the default branch will not install.

The ``pyproject.toml`` also carries a pixi setup (``pixi install``, then ``pixi run``), which installs the Python dependencies but not ``cdo`` or NCO.

Building these docs
===================

.. code-block:: bash

  pip install -r docs/requirements.txt
  sphinx-build -W -b html docs/source docs/build/html

The docs do not need the package's own dependencies. ``docs/source/conf.py`` mocks ecCodes, NetCDF, cartopy, pyfesom2 and the rest for autodoc, which is also how Read the Docs builds them. If you add a module that imports a new third-party package, add it to ``autodoc_mock_imports`` there or the API page loses that module with a warning.
