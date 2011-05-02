@echo off

echo Installing Debug...

rem Example calls:
rem > install.cmd               # Install all dlls in the directory of the executable to be debugged
rem > install.cmd flow2d3d      # Install only project flow2d3d (and its dependencies)

rem 0. defaults:
set project=

rem  The next statement is needed in order for the set commands to work inside the if statement
setlocal enabledelayedexpansion

if [%1] EQU [] (
    rem Install the package/engine specified by the first argument.

    set project=%1
    echo Source          : package/engine !project!
) else (
    rem Install all engines.

    set project=install_all
    echo Source          : all engines
)

rem Change to directory where this batch file resides (necessary when oss-install.cmd is called from outside of oss/trunk/src)
cd %~dp0

call :!project!

goto end

rem  Actual install "routines"


rem ===============
rem === INSTALL_ALL
rem ===============
:install_all
    echo "    installing all open source projects (debug) . . ."

    call :deltares_hydro
    call :flow2d3d
    call :wave
    call :plugin_culvert
    call :plugin_delftflow_trafrm
    call :datsel
    call :kubint
    call :lint
    call :mormerge
    call :vs

goto :endproc



rem ==========================
rem === INSTALL_DELTARES_HYDRO
rem ==========================
:deltares_hydro
    echo "installing deltares_hydro . . ."
    echo "... nothing to be done"
goto :endproc



rem ====================
rem === INSTALL_FLOW2D3D
rem ====================
:flow2d3d
    echo "installing flow2d3d . . ."

    set dest_bin="engines_gpl\deltares_hydro\bin\Debug"

    if not exist !dest_bin!     mkdir !dest_bin!

    copy engines_gpl\flow2d3d\bin\Debug\flow2d3d.dll                                     !dest_bin!
    copy engines_gpl\flow2d3d\bin\Debug\flow2d3d_sp.dll                                  !dest_bin!
    copy third_party_open\DelftOnline\lib\Debug\DelftOnline.dll                          !dest_bin!
    copy third_party_open\DelftOnline\lib\Debug\DelftOnlineJNI.dll                       !dest_bin!
    copy third_party_open\DelftOnline\lib\Debug\JavaLaunch.dll                           !dest_bin!
    copy third_party_open\pthreads\bin\win32\pthreadVCE2.dll                               !dest_bin!
    copy third_party_open\pthreads\bin\win32\pthreadvce.dll                                !dest_bin!
    copy third_party_open\mpich2\bin\*.exe                                                 !dest_bin!
    copy third_party_open\mpich2\lib\*.dll                                                 !dest_bin!
goto :endproc



rem ================
rem === INSTALL_WAVE
rem ================
:wave
    echo "installing wave . . ."
    echo "... nothing to be done"
goto :endproc



rem ==========================
rem === INSTALL_PLUGIN_CULVERT
rem ==========================
:plugin_culvert
    echo "installing plugin_culvert . . ."

    set dest_bin="engines_gpl\deltares_hydro\bin\Debug"

    if not exist !dest_bin!     mkdir !dest_bin!

    copy plugins_lgpl\plugin_culvert\bin\Debug\plugin_culvert.dll                        !dest_bin!
goto :endproc



rem ===================================
rem === INSTALL_PLUGIN_DELFTFLOW_TRAFRM
rem ===================================
:plugin_delftflow_trafrm
    echo "installing plugin_delftflow_trafrm . . ."

    set dest_bin="engines_gpl\deltares_hydro\bin\Debug"

    if not exist !dest_bin!     mkdir !dest_bin!

    copy plugins_lgpl\plugin_delftflow_traform\bin\Debug\plugin_delftflow_traform.dll    !dest_bin!
goto :endproc



rem ==================
rem === INSTALL_DATSEL
rem ==================
:datsel
    echo "installing datsel . . ."
    echo "... nothing to be done"
goto :endproc



rem ==================
rem === INSTALL_KUBINT
rem ==================
:kubint
    echo "installing kubint . . ."
    echo "... nothing to be done"
goto :endproc



rem ================
rem === INSTALL_LINT
rem ================
:lint
    echo "installing lint . . ."
    echo "... nothing to be done"
goto :endproc



rem ====================
rem === INSTALL_MORMERGE
rem ====================
:mormerge
    echo "installing mormerge . . ."
    echo "... nothing to be done"
goto :endproc



rem ==============
rem === INSTALL_VS
rem ==============
:vs
    echo "installing vs . . ."
    echo "... nothing to be done"
goto :endproc






:end
if NOT %ErrorLevel% EQU 0 (
    rem
    rem Only jump to :end when the script is completely finished
    rem 
    exit %ErrorLevel%
)

:endproc
   rem
   rem No exit here
   rem Otherwise the script exits directly at the first missing artefact
