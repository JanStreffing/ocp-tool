[![Documentation Status](https://readthedocs.org/projects/ocp-tool/badge/?version=latest)](https://ocp-tool.readthedocs.io/en/latest/?badge=latest)

# OCP-Tool

OpenIFS Coupling Preparation tool. It fits the OpenIFS land-sea mask to a FESOM2 mesh and writes the files that have to agree with it: the modified OpenIFS initial files, the OASIS3-MCT grids, masks and areas, the runoff mapper input and the LPJ-GUESS soil type file.

Documentation: https://ocp-tool.readthedocs.io

```console
$ conda env create -f environment.yaml
$ conda activate ocp-tool2
$ python run_ocp_tool.py configs/TCO95_CORE3.yaml
```

`environment.yaml` does not yet list everything the tool needs; see the installation page of the documentation before your first run.
