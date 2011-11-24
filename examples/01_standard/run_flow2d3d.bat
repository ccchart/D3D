@ echo off
    rem
    rem This script is an example for running Delft3D-FLOW
    rem Adapt and use it for your own purpose
    rem
    rem adri.mourits@deltares.nl
    rem 27 Dec 2010
    rem 
    rem
    rem This script starts a single-domain Delft3D-FLOW computation on Windows
    rem


    rem
    rem Set the config file here
    rem 
set argfile=config_flow2d3d.ini





    rem
    rem Set the directory containing delftflow.exe here
    rem
set D3D_HOME=..\..\bin
set exedir=%D3D_HOME%\w32\flow\bin
set libdir=%D3D_HOME%\w32\lib

    rem
    rem No adaptions needed below
    rem

    rem Set some (environment) parameters
set PATH=%exedir%;%libdir%;%PATH%

    rem Run
%exedir%\deltares_hydro.exe %argfile%


    rem To prevent the DOS box from disappearing immediately: remove the rem on the following line
pause
