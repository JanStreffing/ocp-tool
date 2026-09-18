#!/bin/bash

set echo
set verbose

if [[ ${host:-`hostname`} == levante* ]]; then
  # On levante: GNU Fortran compiler (gcc@11.2.0) is used for calnoro programme
  NCROOT=/sw/spack-levante/netcdf-fortran-4.5.3-l2ulgp
  module load gcc/11.2.0-gcc-11.2.0
  module load netcdf-fortran/4.5.3-gcc-11.2.0
  gfortran -g -O -I$NCROOT/include -o calnoro calnoro.f90 grid_noro.f90 -L $NCROOT/lib -lnetcdff
elif [[ ${host:-`hostname`} == albedo* ]]; then
  # On albedo: GNU Fortran compiler is used for calnoro programme
  module load gcc/12.1.0
  module load netcdf-fortran/4.5.4-gcc12.1.0 
  NCROOT=/albedo/soft/sw/spack-sw/netcdf-fortran/4.5.4-hdggb4u
  gfortran -g  -O -I$NCROOT/include -o calnoro calnoro.f90 grid_noro.f90  -L $NCROOT/lib -lnetcdff
else
  echo " "
   echo "The system " ${host:-`hostname`} " is not supported."
  echo " "
  exit
fi
