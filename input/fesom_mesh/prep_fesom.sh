#!/bin/bash
#exp: ./prep_fesom.sh core2 r3600x1801 ../openifs_input_default/ICMGGhf05INIT

mesh_name=$1
regular_resolution=$2
oifs_icmgg_file=$3

module load cdo nco

cdo genycon,${regular_resolution} -selname,cell_area -setgrid,${mesh_name}_griddes_nodes.nc ${mesh_name}_griddes_nodes.nc weights_cell_area_${regular_resolution}.nc

cdo -L -remap,${regular_resolution},weights_cell_area_${regular_resolution}.nc -selname,cell_area -setgrid,${mesh_name}_griddes_nodes.nc ${mesh_name}_griddes_nodes.nc ${mesh_name}_regular.nc

ncap2 -O -s 'where(cell_area>0.) cell_area=0;' ${mesh_name}_regular.nc ${mesh_name}_oce.nc
cdo setmisstoc,1 ${mesh_name}_oce.nc ${mesh_name}_land.nc

cdo griddes $oifs_icmgg_file > griddes.txt
cdo -setgrid,$oifs_icmgg_file -remapdis,griddes.txt ${mesh_name}_land.nc ${mesh_name}_oifs.nc
