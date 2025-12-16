Installation
============

Requirements
------------

- Python 3.9+
- NetCDF libraries
- ECCODES for GRIB processing
- CDO (Climate Data Operators)
- NCO (NetCDF Operators)

Quick Install
-------------

Using conda/mamba::

    git clone https://github.com/JanStreffing/ocp-tool.git
    cd ocp-tool
    conda env create -f environment.yaml
    conda activate ocp-tool
    pip install -e .

Using pixi::

    git clone https://github.com/JanStreffing/ocp-tool.git
    cd ocp-tool
    pixi install
