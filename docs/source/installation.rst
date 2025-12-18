Installation
============

Requirements
------------

- Python 3.9+
- NetCDF libraries
- ECCODES for GRIB processing
- CDO (Climate Data Operators)
- NCO (NetCDF Operators)

Quick Install (Recommended)
---------------------------

Clone from GitHub and install dependencies manually::

    git clone https://github.com/JanStreffing/ocp-tool.git
    cd ocp-tool

Then install the required Python packages in your existing environment.
This is the most common approach, especially on HPC systems where
you may already have the dependencies available via modules.

Using conda/mamba
-----------------

Create a dedicated environment with all dependencies::

    git clone https://github.com/JanStreffing/ocp-tool.git
    cd ocp-tool
    conda env create -f environment.yaml
    conda activate ocp-tool
    pip install -e .

Using pixi
----------

Pixi handles all dependencies automatically::

    git clone https://github.com/JanStreffing/ocp-tool.git
    cd ocp-tool
    pixi install
