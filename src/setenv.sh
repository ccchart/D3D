#!/bin/bash
###############################################
### load your build environment 	    ###
###############################################
echo "Module Load"

if [ "$1" == "intel23" ]; then
     echo "Loading Intel23 compiled modules"
  
     module load    intel/2023.1.0
     module display intel/2023.1.0
 
     module load    intelmpi/2021.9.0
     module display intelmpi/2021.9.0

      . $SETVARS_VARS_PATH -ofi_internal=1
 
     module load    netcdf/4.9.2_gcc12.2.0
     module display netcdf/4.9.2_gcc12.2.0
  
     module load    petsc/3.19.0_gcc12.2.0
     module display petsc/3.19.0_gcc12.2.0
  
     #module load    metis/5.1.0_intel21.2.0
     #module display metis/5.1.0_intel21.2.0
  
     module load    cmake/3.26.1_intel2023.1.0
     module display cmake/3.26.1_intel2023.1.0

     module load    gcc/12.2.0_gcc12.2.0
     module display gcc/12.2.0_gcc12.2.0
      
     module load    proj/9.2.0_intel2023.1.0
     module display proj/9.2.0_intel2023.1.0

     module load    gdal/3.6.3_intel2023.1.0
     module display gdal/3.6.3_intel2023.1.0

     module load    svn/1.14.2_gcc12.2.0
     module display svn/1.14.2_gcc12.2.0

     module load    patchelf/0.17.2_gcc12.2.0
     module display patchelf/0.17.2_gcc12.2.0
elif [ "$1" == "intel21" ]; then
     echo "Loading Intel21 compiled modules"
  
     module load    intel/21.2.0
     module display intel/21.2.0
 
     module load    intelmpi/21.2.0
     module display intelmpi/21.2.0

      . $SETVARS_VARS_PATH -ofi_internal=1
 
     module load    netcdf/v4.7.4_v4.5.3_intel21.2.0
     module display netcdf/v4.7.4_v4.5.3_intel21.2.0
  
     module load    petsc/3.13.3_intel21.2.0_intelmpi21.2.0_no_mkl
     module display petsc/3.13.3_intel21.2.0_intelmpi21.2.0_no_mkl
  
     module load    metis/5.1.0_intel21.2.0
     module display metis/5.1.0_intel21.2.0
  
     module load    cmake/3.19.3_intel21.2.0 
     module display cmake/3.19.3_intel21.2.0 
             
     # Shapelib is intertangled with the code in third_party_open
     # loading the module is useless
     #module load    shapelib/1.5.0_intel18.0.3
     #module display shapelib/1.5.0_intel18.0.3

     module load    gcc/7.3.0
     module display gcc/7.3.0
      
     module load    proj/7.1.0_gcc7.3.0
     module display proj/7.1.0_gcc7.3.0

     module load    gdal/3.1.2_gcc7.3.0
     module display gdal/3.1.2_gcc7.3.0

     module load    svn/1.9.12serf_gcc7.3.0
     module display svn/1.9.12serf_gcc7.3.0

     module load    patchelf/0.12
     module display patchelf/0.12
else 
     echo "Loading Intel18 compiled modules"
  
     module load    intel/18.0.3
     module display intel/18.0.3
    
     module load    mpich/3.3.2_intel18.0.3
     module display mpich/3.3.2_intel18.0.3
  
     module load    netcdf/v4.7.4_v4.5.3_intel18.0.3
     module display netcdf/v4.7.4_v4.5.3_intel18.0.3
  
     module load    petsc/3.13.3_intel18.0.3_mpich3.3.2
     module display petsc/3.13.3_intel18.0.3_mpich3.3.2
  
     module load    metis/5.1.0_intel18.0.3
     module display metis/5.1.0_intel18.0.3
  
     module load    cmake/3.18.0_intel18.0.3 
     module display cmake/3.18.0_intel18.0.3 
     
     # Shapelib is intertangled with the code in third_party_open
     # loading the module is useless
     #module load    shapelib/1.5.0_intel18.0.3
     #module display shapelib/1.5.0_intel18.0.3

     module load    gcc/7.3.0
     module display gcc/7.3.0
      
     module load    proj/7.1.0_gcc7.3.0
     module display proj/7.1.0_gcc7.3.0

     module load    gdal/3.1.2_gcc7.3.0
     module display gdal/3.1.2_gcc7.3.0

     module load    svn/1.9.12serf_gcc7.3.0
     module display svn/1.9.12serf_gcc7.3.0

     module load    patchelf/0.12
     module display patchelf/0.12
fi


echo "Export environment variables"
if [ "$1" == "intel23" ]; then
     export FC=ifx
     export CXX=icx
     export CC=icc
elif [ "$1" == "intel21" ]; then
     export FC=mpiifort
     export CXX=mpiicpc
     export CC=mpiicc
else
     export FC=mpif90
     export CXX=mpicxx
     export CC=mpicc
fi
echo "FC=$FC"
echo "CXX=$CXX"
echo "CC=$CC"
