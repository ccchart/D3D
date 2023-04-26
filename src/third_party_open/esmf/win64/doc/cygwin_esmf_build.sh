#!/usr/bin/env bash

export ESMF_DIR=/home/westrene/git/esmf
export ESMF_COMM=mpiuni
export ESMF_COMPILER=gfortran
export ESMF_NETCDF=nc-config
export ESMF_PIO=OFF

# Build DLL
make -j4 lib chkdir_apps

# Build ESMF_RegridWeightGen.exe
make -C src/apps/ESMF_RegridWeightGen

# Create dir for binaries/libraries
mkdir esmf_bin

# Copy all binaries/libraries
cp `PATH=lib/libO/Cygwin.${ESMF_COMPILER}.64.${ESMF_COMM}.default/:$PATH ldd apps/apps0/Cygwin.${ESMF_COMPILER}.64.${ESMF_COMM}.default/ESMF_RegridWeightGen.exe | awk '{ print $3 }' | grep -e ^/usr -e ^/home` esmf_bin

# Copy the application itself
cp apps/appsO/Cygwin.${ESMF_COMPILER}.64.${ESMF_COMM}.default/ESMF_RegridWeightGen.exe esmf_bin