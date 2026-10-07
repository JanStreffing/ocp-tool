#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from os import makedirs

from .config import truncation_prefix

def create_outputdirs(config, resolution):
    """Create necessary output directories based on configuration and resolution."""
    # Determine prefix based on truncation type
    prefix = truncation_prefix(config.atmosphere.truncation_type)
    
    output_dir = f"./output/{prefix}{resolution}_{config.ocean.grid_name}/"
    makedirs(f"{output_dir}lpj-guess", exist_ok=True)
    makedirs(f"{output_dir}oasis_mct3_input", exist_ok=True)
    makedirs(f"{output_dir}openifs_input_modified", exist_ok=True)
    makedirs(f"{output_dir}plots", exist_ok=True)
    makedirs(f"{output_dir}runoff_map_modified", exist_ok=True)