Usage
=====

Configuration
-------------

OCP-Tool uses YAML configuration files. Example configs are in the ``configs/`` directory:

- ``TCO95_CORE2.yaml`` - TCO95 atmosphere with CORE2 ocean mesh
- ``TCO319_CORE3.yaml`` - TCO319 atmosphere with CORE3 ocean mesh (with ice cavities)

Running the Tool
----------------

Basic usage::

    python run_ocp_tool.py configs/TCO95_CORE2.yaml

Configuration Options
---------------------

Atmosphere settings::

    atmosphere:
      resolution_list: [95]  # TCO95
      truncation_type: "cubic-octahedral"
      experiment_name: "ab45"  # ICMGG file prefix

Ocean settings::

    ocean:
      grid_name: "CORE2"
      has_ice_cavities: false
      mesh_file: "/path/to/mesh.nc"

Output Structure
----------------

Output is organized by grid combination::

    output/
    └── TCO95_CORE2/
        ├── lpj-guess/
        ├── oasis_mct3_input/
        ├── openifs_input_modified/
        ├── plots/
        └── runoff_map_modified/
