#!/bin/bash


# Reading command line arguments
mesh_file=$1
mesh_name=$2
regular_resolution=$3
oifs_icmgg_file=$4
cavity=$5


# Display help message
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
    echo "Usage: $0 <mesh_file> <output_dir> <mesh_resolution> <input_dir> <cavity_option>"
    echo "Example: $0 /work/ab0246/a270092/input/fesom2/PI_ICEv2/mesh.nc PI_ICEv2 r360x181 ../openifs_input_default/ICMGGhf05INIT true"
    exit 0
fi


# Check the number of arguments
if [ "$#" -ne 5 ]; then
    echo "Error: Incorrect number of arguments. Use -h or --help for usage."
    exit 1
fi


# Check if cdo and nco are in PATH
if ! command -v cdo &> /dev/null; then
    echo "cdo is not found in PATH. Loading module..."
    module load cdo
fi
if ! command -v nco &> /dev/null; then
    echo "nco is not found in PATH. Loading module..."
    module load nco
fi

# Cleanup to avoid HDF5 error warning while overwriting
rm -f cell_area.nc cav_nod_mask.nc cav_nod_mask_step1.nc cav_nod_mask_step2.nc cell_area_times_cav_nod_mask.nc weights_cell_area_${regular_resolution}.nc 
rm -f cell_area_times_cav_nod_mask.nc ${mesh_name}_regular.nc  ${mesh_name}_oce.nc  ${mesh_name}_land.nc griddes.txt ${mesh_name}_oifs.nc


# Generate remapping weights from FESOM mesh to regular grid
# If with cavity, ensure cell_area for atm<->oce exchange = fill_value over ice shelf region
if [ "${cavity,,}" = "true" ]; then #forcing lower case
  cdo -selvar,cell_area ${mesh_file} cell_area.nc
  cdo chname,cav_nod_mask,cell_area -selvar,cav_nod_mask ${mesh_file} cav_nod_mask.nc
  
  ncap2 -O -s 'cell_area = 1 - cell_area' cav_nod_mask.nc cav_nod_mask_step1.nc
  ncap2 -O -s 'where(cell_area<0.5) cell_area=-1;' cav_nod_mask_step1.nc cav_nod_mask_step2.nc
  cdo mul cell_area.nc cav_nod_mask_step2.nc cell_area_times_cav_nod_mask.nc

  echo "cdo genycon,${regular_resolution} -selname,cell_area -setgrid,${mesh_file} cell_area_times_cav_nod_mask.nc weights_cell_area_${regular_resolution}.nc"
  cdo genycon,${regular_resolution} -selname,cell_area -setgrid,${mesh_file} cell_area_times_cav_nod_mask.nc weights_cell_area_${regular_resolution}.nc

  # Remap FESOM2 mesh to intermediary regular grid, of arbitrary but high resolution
  echo "cdo -L -remap,${regular_resolution},weights_cell_area_${regular_resolution}.nc -selname,cell_area -setgrid,${mesh_file} cell_area_times_cav_nod_mask.nc ${mesh_name}_regular.nc"
  cdo -L -remap,${regular_resolution},weights_cell_area_${regular_resolution}.nc -selname,cell_area -setgrid,${mesh_file} cell_area_times_cav_nod_mask.nc ${mesh_name}_regular.nc
else
  echo "cdo genycon,${regular_resolution} -selname,cell_area -setgrid,${mesh_file} ${mesh_file} weights_cell_area_${regular_resolution}.nc"
  cdo genycon,${regular_resolution} -selname,cell_area -setgrid,${mesh_file} ${mesh_file} weights_cell_area_${regular_resolution}.nc

  # Remap FESOM2 mesh to intermediary regular grid, of arbitrary but high resolution
  echo "cdo -L -remap,${regular_resolution},weights_cell_area_${regular_resolution}.nc -selname,cell_area -setgrid,${mesh_file} ${mesh_file} ${mesh_name}_regular.nc"
  cdo -L -remap,${regular_resolution},weights_cell_area_${regular_resolution}.nc -selname,cell_area -setgrid,${mesh_file} ${mesh_file} ${mesh_name}_regular.nc
fi


# Create binary field where cell_area is 0 / >0
ncap2 -O -s 'where(cell_area>0.) cell_area=0;' ${mesh_name}_regular.nc ${mesh_name}_oce.nc
cdo setmisstoc,1 ${mesh_name}_oce.nc ${mesh_name}_land.nc

# Remap the binary regular grid mask to OpenIFS grid
cdo griddes $oifs_icmgg_file > griddes.txt
cdo -setgrid,$oifs_icmgg_file -remapdis,griddes.txt ${mesh_name}_land.nc ${mesh_name}_oifs.nc
