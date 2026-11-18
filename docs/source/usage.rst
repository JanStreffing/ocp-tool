*****
Usage
*****

One run of the tool handles one atmosphere grid against one ocean grid. Pick the config in ``configs/`` that is closest to what you want, copy it, change the grid settings, and run it from the top of the repository:

.. code-block:: bash

  cd ocp-tool
  python run_ocp_tool.py configs/TCO95_CORE3.yaml

A TCO95 run takes about two minutes. TCO319 and above take long enough that you want a batch job.

Use ``run_ocp_tool.py``, not ``python -m ocp_tool``. The module entry point is an older copy of the pipeline that stops after step 12: it never attempts the weights and silently skips every paleo step, whatever the config says.

The configuration
=================

Keep the config in ``configs/``. Relative paths in it are resolved against the directory above the one the config sits in, so a config saved somewhere else looks for ``input/`` and writes ``output/`` next to wherever that is. Set ``paths.root_dir`` if you really need it elsewhere.

The keys that decide what you get:

- ``atmosphere.resolution_list``: the spectral truncation, e.g. ``[95]``. Give one value. The loop accepts several, but the output directory and the OASIS grid writer only look at the first, so a second resolution overwrites the first one's files.
- ``atmosphere.truncation_type``: ``cubic-octahedral`` (TCO), ``linear`` (TL) or ``quadratic`` (TQ).
- ``atmosphere.experiment_name``: the four character id of the OpenIFS input files. ``ab45`` means the tool reads ``ICMGGab45INIT`` and ``ICMGGab45INIUA`` from ``input/openifs_input_default/``. It is also written into the GRIB header, see :ref:`awiesm3:paleo` for why that matters.
- ``atmosphere.model_cycle``: ``43r3``, ``48r1`` or ``auto``. The two cycles lay out snow differently. Set it explicitly; the tool then refuses an input file of the other cycle instead of quietly writing the wrong layout. Details at :ref:`awiesm3:model_cycle`.
- ``ocean.grid_name``: the label that goes into the output names, or ``AMIP`` for no ocean. In AMIP mode the mesh, the weights and the runoff steps are skipped.
- ``ocean.mesh_file``: the FESOM2 ``mesh.nc``. If it does not exist, pyfesom2 builds it from the ASCII mesh files in the same directory, see :ref:`trap_mesh_nc`.
- ``ocean.has_ice_cavities``: whether the mesh has ice shelf cavities. It has no default. With it on, cavity nodes come out as land, so the atmosphere sees an ice shelf rather than open water. Keep it in agreement with FESOM2's ``use_cavity`` yourself, because nothing checks the two against each other.
- ``runoff.manual_basin_removal``: closed seas to take out of the runoff routing. ``caspian-sea`` and ``black-sea`` are the only names recognised; anything else is ignored without a message.
- ``lakes.restore_flipped_lakes``: put the input's lake cover back on cells the mesh turned from ocean into land. Without it the Caspian and the Great Lakes lose their FLake lakes and become soil. Off by default only because turning it on changes files that are already in use.
- ``options.output_suffix``: appended to the modified ICMGG and the LPJ-GUESS soil file, ``_v2`` unless set. It exists so that a regenerated file does not overwrite one that running experiments still read under its old name. Set ``""`` for a grid with no published predecessor.
- ``options.generate_rmp``: whether to attempt step 4b at all.
- ``paleo`` and ``ice``: optional sections, off unless ``enabled: true``. Paleo is documented in full at :doc:`awiesm3:paleo`.

Two keys are read and then not used: ``paths.input.fesom_mesh`` and ``ocean.intermediate_resolution``. Changing them does nothing; the mesh comes from ``ocean.mesh_file`` only.

What comes out
==============

Everything goes to ``output/{TCO|TL|TQ}{resolution}_{grid_name}/``:

- ``openifs_input_modified/``: ``ICMGG{expid}INIT_{grid}{suffix}`` with the new land-sea mask, the surface fields of every flipped cell rebuilt from its neighbours, CO2 emissions and bare soil albedo added; ``ICMGG{expid}INIUA`` with 3D CO2 added. Paleo runs add ``ICMSH`` and ``ICMCL`` files.
- ``oasis_mct3_input/``: ``grids.nc``, ``masks.nc``, ``areas.nc``, the CO2 coupling restarts ``rst_co2_ao.nc`` and ``rst_co2_av.nc``, and ``vegin.nc`` for LPJ-GUESS.
- ``runoff_map_modified/``: the runoff mapper input, with the closed seas from ``runoff.manual_basin_removal`` rerouted and the arrival points corrected.
- ``lpj-guess/``: ``slt_{TRUNC}{res}_{grid}{suffix}.nc``, the soil type file.
- ``plots/``: a map of the land points, worth a look before anything else.

Where each of these goes in a runscript is described at :doc:`awiesm3:paleo` (the paths are the same for non-paleo setups) and :doc:`awiesm3:how_to/coupling_oasis`.

The remapping weights
---------------------

Step 4b tries to build the ``rmp_`` files with pyOASIS. It needs the FESOM2 grid in ``grids.nc``, which the tool only writes with ``ocean.write_oasis_grid: true``, and an OASIS build under ``$OASIS_BUILD_PATH``. No production config sets either, so the normal outcome is this, after which the run carries on:

.. code-block:: bash

  Warning: Failed to generate rmp files: "grid 'feom' not found in .../oasis_mct3_input/.wg_link00/grids.nc"

That is expected and not a failure. The ``.wg_link00/`` directory it leaves behind is scratch and can be deleted. Production weights come from the offline weight generation in :doc:`awiesm3:how_to/coupling_oasis`, using the ``tools/submit_weightgen_*.sh`` scripts here. Do not reuse weights from another grid pair: a new truncation or a new mesh is a new grid pair. Moving an existing set to a different FESOM2 core count is a reordering job, not a regeneration, and the same page covers it.

Check that it worked
====================

``OCP-Tool processing complete!`` at the end proves nothing. The tool catches failures in several steps, prints them and carries on (see :ref:`trap_ncrename`), so read the log for ``Error`` and ``Warning`` lines other than the 4b one above.

Then check that the mask came from the mesh you meant. The log prints ``Cavity nodes:`` for a cavity mesh and ``Number of land points:`` after step 3. For TCO95, from runs that are known to be right:

- CORE2: ``Number of land points: 11399``
- CORE3 with cavities: ``Cavity nodes: 10047``, ``Number of land points: 11312``
- CORE3 with the cavities cut out of the mesh (``core3_nocav``): ``Number of land points: 11296``

A count that matches a different mesh than the one you configured usually means a stale ``mesh.nc``. Finally, open the land point plot in ``plots/`` and look at the coasts you care about.

.. _trap_mesh_nc:

A stale mesh.nc is reused
=========================

If ``ocean.mesh_file`` exists it is read as it is. The tool does not compare it with the ASCII mesh files next to it, so a ``mesh.nc`` left over from an earlier version of the mesh, or built with the other cavity setting, gives you the old geography without a warning. Set ``ocean.force_overwrite_griddes: true`` to rebuild it.

The rebuild writes ``mesh.nc`` into the mesh directory itself, so you need write permission there. If you do not have it (a shared pool directory), copy or link the mesh files into a directory of your own and point ``mesh_file`` there. A rebuild that fails stops the run with ``Failed to generate mesh.nc via pyfesom2``.

pyfesom2 warns ``coast information from nod2d.out seems to be corrupt, using diagnosed coast flag instead`` on some meshes, CORE3 among them. The tool takes the coastline from the triangulation edges and not from that column, so the warning does not affect the mask.

.. _trap_ncrename:

The LPJ-GUESS soil file needs NCO
=================================

Step 12 renames the soil type variable with ``ncrename``. Without NCO on ``PATH`` it fails like this, and the run then reports success:

.. code-block:: bash

  sh: ncrename: command not found
  Error: ncrename failed (exit 32512). Is NCO on PATH?

The ``slt_*.nc`` file is left behind with the variable still called ``slt`` instead of ``var43``. LPJ-GUESS falls back to ``slt``, so a run with that file starts, and the mistake surfaces much later if at all. Do not use it. On Levante, ``module load nco/5.0.6-gcc-11.2.0`` and run the tool again.

Run one config at a time per directory
======================================

Step 1 writes the grid description to ``griddes.txt`` in the current directory and reads it straight back. Two runs started from the same directory at the same time can read each other's grid. Run them one after the other, or from different working directories.

Known limitations
=================

- Only one resolution per config, see above.
- ``environment.yaml`` does not list cartopy, shapely, xarray or NCO, which the tool needs. :doc:`installation` gives the extra install line.
- The NEMO path (``ocp_tool.scriptengine_task``, run as a ScriptEngine task for EC-Earth) is a separate pipeline built on ``ocp_tool.grids`` and is not described here.
- ``ocp_tool.dynamic_regen`` reruns the whole pipeline for a moving FESOM2 submesh in ice sheet coupled runs. It is driven by the coupling scripts, and its module docstring in :doc:`api` is the reference.
