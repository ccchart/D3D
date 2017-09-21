@ echo off
    rem
    rem This script is an example for running Delft3D-WAVE
    rem Adapt and use it for your own purpose
    rem
    rem adri.mourits@deltares.nl
    rem 17 sep 2014
    rem 
    rem
    rem This script starts a Delft3D-WAVE calculation on Windows
    rem





    rem
    rem Set the directory containing wave.exe and swan
    rem
set ARCH=win64
set D3D_HOME=..\..\bin
  rem set D3D_HOME=c:\Program Files (x86)\Deltares\Delft3D 4.01.00
set dimrexedir=%D3D_HOME%\%ARCH%\dimr\bin
set waveexedir=%D3D_HOME%\%ARCH%\wave\bin
set swanexedir=%D3D_HOME%\%ARCH%\swan\bin
set swanbatdir=%D3D_HOME%\%ARCH%\swan\scripts
set shareddir=%D3D_HOME%\%ARCH%\shared

    rem
    rem No adaptions needed below
    rem


    rem Run
title Wave simulation
set PATH=%dimrexedir%;%shareddir%;%waveexedir%;%swanbatdir%;%swanexedir%;%PATH%
"%dimrexedir%\dimr.exe" dimr_config.xml
title %CD%

    rem To prevent the DOS box from disappearing immediately: remove the rem on the following line
rem pause
