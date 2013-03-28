#!/bin/bash


echo "Running testcase 01_standard ..."
cd 01_standard
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d" "<library>flow2d3d_sp"
./run_flow2d3d.sh >screen.log 2>&1


echo "Running testcase 01_standard parallel ..."
./run_flow2d3d_parallel.sh >screen_parallel.log 2>&1
   # Undo changes
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d_sp" "<library>flow2d3d"
cd ..


echo "Running testcase 02_domaindecomposition ..."
cd 02_domaindecomposition
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d" "<library>flow2d3d_sp"
./run_flow2d3d.sh >screen.log 2>&1
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d_sp" "<library>flow2d3d"
cd ..


echo "Running testcase 03_flow-wave ..."
cd 03_flow-wave
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d" "<library>flow2d3d_sp"
./run_flow2d3d.sh >screen.log 2>&1
../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d_sp" "<library>flow2d3d"
cd ..


echo "Running testcase 04_fluidmud ..."
cd 04_fluidmud
../sed_in_file.tcl config_d_hydro_mud.xml "<library>flow2d3d" "<library>flow2d3d_sp"
../sed_in_file.tcl config_d_hydro_sed.xml "<library>flow2d3d" "<library>flow2d3d_sp"
./run_flow2d3d_flm.sh >screen.log 2>&1
../sed_in_file.tcl config_d_hydro_mud.xml "<library>flow2d3d_sp" "<library>flow2d3d"
../sed_in_file.tcl config_d_hydro_sed.xml "<library>flow2d3d_sp" "<library>flow2d3d"
cd ..


echo "Running testcase 05_mormerge ..."
cd 05_mormerge/input
../../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d" "<library>flow2d3d_sp"
cd ../merge
./run_flow2d3d_wave_mormerge.sh >screen.log 2>&1
cd ../input
../../sed_in_file.tcl config_d_hydro.xml "<library>flow2d3d_sp" "<library>flow2d3d"
cd ../..


echo ...finished

