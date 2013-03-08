@ echo off
    rem
    rem This script is an example for running Delft3D-FLOW
    rem Adapt and use it for your own purpose
    rem
    rem adri.mourits@deltares.nl
    rem 01 Mar 2011
    rem 
    rem
    rem This script starts a single-domain Delft3D-FLOW computation online with Delft3D-WAVE on Windows
    rem


    rem
    rem Set the config file and mdw file
    rem 
set argfile=config_d_hydro.xml
set mdwfile=r17.mdw




    rem
    rem Set the directory containing delftflow.exe
    rem
set D3D_HOME=..\..\bin\win32
set flowexedir=%D3D_HOME%\flow2d3d\bin
set waveexedir=%D3D_HOME%\wave\bin
set swanexedir=%D3D_HOME%\swan\bin
set swanbatdir=%D3D_HOME%\swan\scripts

    rem
    rem No adaptions needed below
    rem


    rem Run
set PATH=%flowexedir%;%PATH%
start %flowexedir%\d_hydro.exe %argfile%

set PATH=%waveexedir%;%swanbatdir%;%swanexedir%;%PATH%
%waveexedir%\wave.exe %mdwfile% 1

    rem To prevent the DOS box from disappearing immediately: remove the rem on the following line
rem pause
