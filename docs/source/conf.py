"""Sphinx configuration for OCP-Tool documentation."""

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

project = 'OCP-Tool'
copyright = '2024, AWI Climate Dynamics'
author = 'Jan Streffing'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'alabaster'
html_static_path = ['_static']

# Napoleon settings for Google/NumPy style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = True
