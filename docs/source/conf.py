"""Sphinx configuration for the OCP-Tool documentation."""

import os
import re
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, ROOT)

# -- Project information -----------------------------------------------------

project = 'OCP-Tool'
copyright = '2020-2026, Jan Streffing'
author = 'Jan Streffing'

# Read the version from the package without importing it, so that the version
# is right even if an import fails.
with open(os.path.join(ROOT, 'ocp_tool', '__init__.py')) as f:
    release = re.search(r'^__version__ = "([^"]+)"', f.read(), re.M).group(1)
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Autodoc -----------------------------------------------------------------

# Read the Docs does not install the package. These are imported by ocp_tool
# but are either heavy, need system libraries (ecCodes, NetCDF, PROJ, GEOS,
# MPI) or are not on PyPI in a usable form, so autodoc gets stand-ins.
autodoc_mock_imports = [
    'cartopy',
    'cdo',
    'dask',
    'eccodes',
    'gribapi',
    'matplotlib',
    'mpi4py',
    'netCDF4',
    'pyfesom2',
    'pyoasis',
    'pyproj',
    'scipy',
    'scriptengine',
    'shapely',
    'xarray',
]

autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'show-inheritance': True,
}
autodoc_typehints = 'description'

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'awiesm3': ('https://awi-esm3.readthedocs.io/en/latest/', None),
}

# -- HTML output -------------------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_context = {
    'display_github': True,
    'github_user': 'JanStreffing',
    'github_repo': 'ocp-tool',
    'github_version': 'master',
    'conf_py_path': '/docs/source/',
}
