@echo off

rem ######
rem # MAIN
rem ######

setlocal enabledelayedexpansion
set globalErrorLevel=0

echo oss-post_build...

rem Usage:
rem > oss-install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration>
rem > oss-install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> [project] 
rem > oss-install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> [project] ["compiler_redist_dir"]
rem > oss-install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> [project] ["compiler_redist_dir"] ["mkl_redist_dir"]
rem > oss-install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> [project] ["compiler_redist_dir"] ["mkl_redist_dir"] ["mpi_redist_dir"]

rem with:
rem   <install_dir>           : Target directory where all binaries etc. are going to be installed by this script
rem   <build_dir>             : The root build directory where all binaries are build
rem   <configuration>         : Debug or release configuration. Configuration determines the directory structure
rem   [project]               : (optional) project to install. If missing, "everything" is installed
rem   ["compiler_redist_dir"] : (optional) Directory containing compiler specific dll's to be installed
rem   ["mkl_redist_dir"]      : (optional) Directory containing Intel math kernel library specific dll's to be installed
rem   ["mpi_redist_dir"]      : (optional) Directory containing Intel MPI library specific dll's to be installed

rem
rem Example calls:
rem > install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration>                                                    # Install entire solution
rem > install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> dflowfm                                            # Install only project dflowfm (and its dependencies)
rem > install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> dflowfm ["compiler_redist_dir"]                    # "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\ia32\compiler\"                                                                                           
rem                                                                                                                                # Install only project dflowfm (and its dependencies)
rem > install.cmd <install_dir> <build_dir> <checkout_src_root> <configuration> dflowfm ["compiler_redist_dir"] ["mkl_redist_dir"] # "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\ia32\compiler\"  "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries\windows\redist\ia32\mkl\"
rem                                                                                                                                # Install only project dflowfm (and its dependencies including mkl required dlls)

rem Example calls:
rem >  oss-post_build.cmd install_dir Debug            # postbuild all debug
rem >  oss-post_build.cmd install_dir Debug dimr       # postbuild debug dimr

rem  The next statement is needed in order for the set commands to work inside the if statement

set install_dir=%1
set build_dir=%2
set checkout_src_root=%3
set configuration=%4
set project=%5

rem substitute backslashes
set "install_dir=!install_dir:/=\!"
set "build_dir=!build_dir:/=\!"
set "checkout_src_root=!checkout_src_root:/=\!"

if [%5] EQU [] (
    rem Install all engines
    set project=install_all
) 

if [%6] EQU [] (
    set compiler_redist_dir=""
) else (
    set compiler_redist_dir_read=%6
    rem Remove leading and trailing quote (")
    rem These quotes MUST be present in argument number 3, because "compiler_redist_dir" may contain white spaces
    set compiler_redist_dir=!compiler_redist_dir_read:~1,-1!
    set "compiler_redist_dir=!compiler_redist_dir:/=\!"
)

if [%7] EQU [] (
    set mkl_redist_dir=""
) else (
    set mkl_redist_dir_read=%7
    rem Remove leading and trailing quote (")
    rem These quotes MUST be present in argument number 4, because "mkl_redist_dir_read" may contain white spaces
    set mkl_redist_dir=!mkl_redist_dir_read:~1,-1!
    set "mkl_redist_dir=!mkl_redist_dir:/=\!"
)

if [%8] EQU [] (
    set mpi_redist_dir=""
) else (
    set mpi_redist_dir_read=%8
    rem Remove leading and trailing quote (")
    rem These quotes MUST be present in argument number 8, because "mpi_redist_dir_read" may contain white spaces
    set mpi_redist_dir=!mpi_redist_dir_read:~1,-1!
    set "mpi_redist_dir=!mpi_redist_dir:/=\!"
)



echo install_dir         : !install_dir!
echo build_dir           : !build_dir!
echo checkout_src_root   : !checkout_src_root!
echo configuration       : !configuration!
echo project             : !project!
echo compiler_redist_dir : !compiler_redist_dir!
echo mkl_redist_dir      : !mkl_redist_dir!





rem Change to directory tree where this batch file resides (necessary when oss-install.cmd is called from outside of oss/trunk/src)
cd %~dp0\..\..

call :!project!

goto end





rem ############
rem # PROCEDURES
rem ############

rem =============================================================
rem === copyFile takes two arguments: the name of the file to ===
rem === copy to the destiny directory                         ===
rem ===                                                       ===
rem === NOTE: errors will be reported and the script will     ===
rem === with an error code after executing the rest of its    ===
rem === statements                                            ===
rem =============================================================
:copyFile
    set fileName=%~1
    set dest=%~2
    rem
    rem "echo f |" is (only) needed when dest does not exist
    rem and does not harm in other cases
    rem
    echo f | xcopy "%fileName%" %dest% /F /Y
    if NOT !ErrorLevel! EQU 0 (
        echo ERROR: while copying "!fileName!" to "!dest!"
    )
goto :endproc



rem ===============
rem === POSTBUILD_ALL
rem ===============
:install_all
    echo " postbuild all open source projects . . ."

    call :dimr
    call :dimr_dll
    call :dflowfm
    call :dflowfm_dll
    call :dfmoutput
    call :waq_plugin_wasteload
    call :delwaq_dll
    call :delwaq1
    call :delwaq2
    call :waq_run_processes
    call :duprol2delwaq
    call :wave
    call :wave_exe
    call :flow2d3d

goto :endproc



rem =============================================================
rem === makeDir accepts one argument: the name of the         ===
rem === directory it will create if it doesn't already exists ===
rem ===                                                       ===
rem === NOTE: errors will be reported and the script will     ===
rem === return with an error code after executing the rest of ===
rem === its statements                                        ===
rem =============================================================
:makeDir
    set dirName=%~1
    if not exist !dirName! mkdir !dirName!
    if not !ErrorLevel! EQU 0 (
        echo ERROR: while creating directory "!dirName!"
    )
goto :endproc



rem =============================================================
rem === makeAllDirs makes all needed directories              ===
rem =============================================================
:makeAllDirs

    call :makeDir !dest_bin!
    call :makeDir !dest_default!
    call :makeDir !dest_scripts!
    call :makeDir !dest_plugins!
    call :makeDir !dest_share!
    
goto :endproc     



rem =============================================================
rem === copies runtime libraries for dflowfm and dflowfm_dll  ===
rem =============================================================
:copyDflowfmDependentRuntimeLibraries

    echo "copyDflowfmDependentRuntimeLibraries . . ."
    
    call :copyFile "!checkout_src_root!\third_party_open\netcdf\netCDF 4.6.1\bin\*"                                          !dest_bin!  
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\*.dll"                                             !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\bin\*.exe"                                               !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\lib\*.dll"                                               !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\expat\x64\x64\%configuration%\libexpat.dll"                         !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\intelredist\lib\x64\\*.*"                                           !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\\*.dll"                                            !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\lib\\*.dll"                                              !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\Tecplot\lib\x64\\*.dll"                                             !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\xerces-c_3_2.dll"                 !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\gdal300.dll"                      !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\expat.dll"                        !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\libpq.dll"                        !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\sqlite3.dll"                      !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\libmysql.dll"                     !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\spatialite.dll"                   !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\proj.dll"                         !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\proj_6_1.dll"                     !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\openjp2.dll"                      !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\geos_c.dll"                       !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\libxml2.dll"                      !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\iconv.dll"                        !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\geos.dll"                         !dest_bin!
    call :copyFile "!checkout_src_root!\third_party_open\GISInternals\release-1911-x64\bin\freexl.dll"                       !dest_bin!
    rem pthreads depends on MSVCR100.DLL
    call :copyFile "!checkout_src_root!\third_party_open\vcredist\x64\Microsoft.VC100.CRT\msvcr100.dll"                      !dest_bin!
    

    if !compiler_redist_dir!=="" (
        rem Compiler_dir not set
    ) else (
        rem "Compiler_dir:!compiler_redist_dir!"
        set localstring="!compiler_redist_dir!*.dll"
        rem Note the awkward usage of !-characters
        call :copyFile !!localstring! !dest_bin!!
        call :copyFile "!checkout_src_root!\third_party_open\petsc\petsc-3.10.2\lib\x64\Release\libpetsc.dll"               !dest_bin!
        rem is needed for dimr nuget package? please check                                   
        call :copyFile "!checkout_src_root!\third_party_open\petsc\petsc-3.10.2\lib\x64\Release\libpetsc.dll"               !dest_share!
    )

    if !mkl_redist_dir!=="" (
        rem mkl_redist_dir not set
    ) else (
        rem note that for onaApi MKL, the DLL names end in '.1.dll'
        set localstring="!mkl_redist_dir!mkl_core*.dll"
        call :copyFile !!localstring! !dest_bin!
        set localstring="!mkl_redist_dir!mkl_def*.dll"
        call :copyFile !!localstring! !dest_bin!
        set localstring="!mkl_redist_dir!mkl_core*.dll"
        call :copyFile !!localstring! !dest_bin!
        set localstring="!mkl_redist_dir!mkl_avx*.dll"
        call :copyFile !!localstring! !dest_bin!
        rem is needed for dimr nuget package? please check
        call :copyFile !!localstring! !dest_share!
        set localstring="!mkl_redist_dir!mkl_intel_thread*.dll"
        call :copyFile !!localstring! !dest_bin!
        rem is needed for dimr nuget package?  please check
        call :copyFile !!localstring! !dest_share!

        rem if 'mkl_redist_dir' contains 'oneAPI', use the version of petsc built with oneAPI Fortran
        if not "x!mkl_redist_dir:oneAPI=!"=="x!mkl_redist_dir!" (
            call :copyFile "!checkout_src_root!\third_party_open\petsc\petsc-3.10.2\lib\x64\Release-oneAPI\libpetsc.dll"             !dest_bin!            
        ) else (
            call :copyFile "!checkout_src_root!\third_party_open\petsc\petsc-3.10.2\lib\x64\Release\libpetsc.dll"             !dest_bin!
        )
    )    

    if !mpi_redist_dir!=="" (
        rem mpi_redist_dir not set
    ) else (
        set localstring="!mpi_redist_dir!bin\release\impi.dll"
        call :copyFile !!localstring! !dest_bin!
        set localstring="!mpi_redist_dir!libfabric\bin\libfabric.dll"
        call :copyFile !!localstring! !dest_bin!
    )

goto :endproc



rem =============================================================
rem === copies runtime libraries for dimr and dimr_lib        ===
rem =============================================================
:copyDimrDependentRuntimeLibraries

    set destination=%~1
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\*.dll"                                        !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\bin\*.exe"                                          !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\lib\*.dll"                                          !destination!
    call :copyFile "!checkout_src_root!\third_party_open\expat\x64\x64\%configuration%\libexpat.dll"                    !destination!

goto :endproc



rem =============================================================
rem === copies runtime libraries for DWaq                     ===
rem =============================================================
:copyDwaqDependentRuntimeLibraries

    set destination=%~1
    call :copyFile "!checkout_src_root!\third_party_open\netcdf\netCDF 4.6.1\bin\*"                                     !destination!  
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\*.dll"                                        !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\bin\*.exe"                                          !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\lib\*.dll"                                          !destination!

goto :endproc



rem =============================================================
rem === copies runtime libraries for DWaves                   ===
rem =============================================================
:copyDwavesDependentRuntimeLibraries

    set destination=%~1
    call :copyFile "!checkout_src_root!\third_party_open\netcdf\netCDF 4.6.1\bin\*"                                     !destination!  
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\*.dll"                                        !destination!

goto :endproc



rem =============================================================
rem === copies runtime libraries for flow2d3d                 ===
rem =============================================================
:copyFlow2D3DDependentRuntimeLibraries

    set destination=%~1
    call :copyFile "!checkout_src_root!\third_party_open\netcdf\netCDF 4.6.1\bin\*"                                     !destination!  
    call :copyFile "!checkout_src_root!\third_party_open\pthreads\bin\x64\*.dll"                                        !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\bin\*.exe"                                          !destination!
    call :copyFile "!checkout_src_root!\third_party_open\mpich2\x64\lib\*.dll"                                          !destination!
    call :copyFile "!checkout_src_root!\third_party_open\expat\x64\x64\%configuration%\libexpat.dll"                    !destination!

goto :endproc



rem ==========================
rem === POST_BUILD_DFLOWFM_DLL
rem ==========================
:dflowfm_dll

    echo "postbuild dflowfm_dll . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDflowfmDependentRuntimeLibraries
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\dflowfm_lib\!configuration!\dflowfm.*"                                           !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dflowfm\bin"
        set dest_default="!install_dir!\x64\Release\dflowfm\default"
        set dest_scripts="!install_dir!\x64\Release\dflowfm\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDflowfmDependentRuntimeLibraries
        
        rem Temporarily rename dest_bin to share_bin to copy libraries there as well
        set dest_bin=!dest_share!
        call :copyDflowfmDependentRuntimeLibraries
        set dest_bin="!install_dir!\x64\Release\dflowfm\bin"
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\dflowfm_lib\!configuration!\dflowfm.dll"                                           !dest_bin! 
    
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_DFLOWFM-CLI
rem ==========================
:dflowfm-cli
    echo "postbuild dflowfm-cli . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"

        call :makeDir !dest_bin!
        call :copyDflowfmDependentRuntimeLibraries
        call :copyFile "!build_dir!\dflowfm_cli_exe\!configuration!\dflowfm-cli.*"                                            !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dflowfm\bin"
        set dest_default="!install_dir!\x64\Release\dflowfm\default"
        set dest_scripts="!install_dir!\x64\Release\dflowfm\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"

        call :makeAllDirs 
        call :copyDflowfmDependentRuntimeLibraries
        
        call :copyFile "!build_dir!\dflowfm_cli_exe\!configuration!\dflowfm-cli.*"                                           !dest_bin!

        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\bloom.spe"                                             !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\bloominp.d09"                                          !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\proc_def.dat"                                          !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\proc_def.def"                                          !dest_default!
        
        call :copyFile "!checkout_src_root!\engines_gpl\dflowfm\scripts\team-city\run_dflowfm_processes.bat"               !dest_scripts!
        call :copyFile "!checkout_src_root!\engines_gpl\dflowfm\scripts\team-city\run_dflowfm.bat"                         !dest_scripts!
    )
    
goto :endproc



rem =================================
rem === POST_BUILD_DFLOWFM_INTERACTER
rem =================================
:dflowfm_interacter
    echo "postbuild dflowfm_interacter . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"

        call :makeDir !dest_bin!
        call :copyDflowfmDependentRuntimeLibraries
        call :copyFile "!build_dir!\dflowfm\!configuration!\dflowfm.*"                                                   !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dflowfm\bin"
        set dest_default="!install_dir!\x64\Release\dflowfm\default"
        set dest_scripts="!install_dir!\x64\Release\dflowfm\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"

        call :makeAllDirs 
        call :copyDflowfmDependentRuntimeLibraries
        
        call :copyFile "!build_dir!\dflowfm\!configuration!\dflowfm.exe"                                                   !dest_bin!

        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\bloom.spe"                                             !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\bloominp.d09"                                          !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\proc_def.dat"                                          !dest_default!
        call :copyFile "!checkout_src_root!\engines_gpl\waq\default\proc_def.def"                                          !dest_default!
        
        call :copyFile "!checkout_src_root!\engines_gpl\dflowfm\scripts\MSDOS\run_dflowfm_processes.bat"                   !dest_scripts!
        call :copyFile "!checkout_src_root!\engines_gpl\dflowfm\scripts\team-city\run_dflowfm.bat"                         !dest_scripts!
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_DFMOUTPUT
rem ==========================
:dfmoutput
    echo "postbuild dfmoutput . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"

        call :makeDir !dest_bin!
        call :copyDflowfmDependentRuntimeLibraries
        call :copyFile "!build_dir!\dfmoutput\!configuration!\dfmoutput.*"                                                    !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dflowfm\bin"
        set dest_default="!install_dir!\x64\Release\dflowfm\default"
        set dest_scripts="!install_dir!\x64\Release\dflowfm\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"

        call :makeAllDirs 
        call :copyDflowfmDependentRuntimeLibraries
        
        call :copyFile "!build_dir!\dfmoutput\!configuration!\dfmoutput.exe"                                                   !dest_bin!

        call :copyFile "!checkout_src_root!\engines_gpl\dflowfm\scripts\team-city\run_dfmoutput.bat"                       !dest_scripts!
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_DIMR
rem ==========================
:dimr
    echo "postbuild dimr . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyDimrDependentRuntimeLibraries                                                                               !dest_bin!
    call :copyFile "!build_dir!\dimr\!configuration!\dimr.*"                                                              !dest_bin!
    )
    
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"

    set dest_bin="!install_dir!\x64\Release\dimr\bin"
    set dest_default="!install_dir!\x64\Release\dimr\default"
    set dest_scripts="!install_dir!\x64\Release\dimr\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"

    call :makeAllDirs 
    call :copyDimrDependentRuntimeLibraries                                                                             !dest_share!
    call :copyFile "!build_dir!\dimr\!configuration!\dimr.exe"                                                          !dest_bin!

    call :copyFile "!checkout_src_root!\engines_gpl\d_hydro\scripts\create_config_xml.tcl"                              !dest_menu!
    call :copyFile "!checkout_src_root!\engines_gpl\dimr\scripts\generic\win64\*.*"                                     !dest_scripts!

    )

goto :endproc



rem ==========================
rem === POST_BUILD_DIMR_LIB
rem ==========================
:dimr_lib
    echo "postbuild dimr_lib . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyDimrDependentRuntimeLibraries                                                                               !dest_bin!
    call :copyFile "!build_dir!\dimr_lib\!configuration!\dimr_dll.*"                                                      !dest_bin!

    )
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_bin="!install_dir!\x64\Release\dimr\bin"
    set dest_default="!install_dir!\x64\Release\dimr\default"
    set dest_scripts="!install_dir!\x64\Release\dimr\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeAllDirs 
    call :copyDimrDependentRuntimeLibraries                                                                               !dest_share!
    call :copyFile "!build_dir!\dimr_lib\!configuration!\dimr_dll.*"                                                      !dest_bin!
    
    )

goto :endproc



rem ===================================
rem === POST_BUILD_waq_plugin_wasteload
rem ===================================
:waq_plugin_wasteload

    echo "postbuild waq_plugin_wasteload . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                                                            !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_plugin_wasteload\!configuration!\waq_plugin_wasteload.*"                           !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\waq_plugin_wasteload\bin"
        set dest_default="!install_dir!\x64\Release\waq_plugin_wasteload\default"
        set dest_scripts="!install_dir!\x64\Release\waq_plugin_wasteload\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                                             !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_plugin_wasteload\!configuration!\waq_plugin_wasteload.*"                            !dest_bin! 
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_delwaq_dll
rem ==========================
:delwaq_dll

    echo "postbuild delwaq_dll . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                                                          !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_plugin_wasteload\!configuration!\waq_plugin_wasteload.*"                         !dest_bin!
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.*"                                                     !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaq\bin"
        set dest_default="!install_dir!\x64\Release\dwaq\default"
        set dest_scripts="!install_dir!\x64\Release\dwaq\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                        !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.dll"                                 !dest_bin! 

        rem default
        call :copyFile !checkout_src_root!\engines_gpl\waq\default\bloom.spe                           !dest_default!
        call :copyFile !checkout_src_root!\engines_gpl\waq\default\bloominp.d09                        !dest_default!
        call :copyFile !checkout_src_root!\engines_gpl\waq\default\proc_def.dat                        !dest_default!
        call :copyFile !checkout_src_root!\engines_gpl\waq\default\proc_def.def                        !dest_default!

    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_delwaq1
rem ==========================
:delwaq1

    echo "postbuild delwaq1 . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                                       !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_plugin_wasteload\!configuration!\waq_plugin_wasteload.*"      !dest_bin!
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.*"                                  !dest_bin!
        call :copyFile "!build_dir!\delwaq1\!configuration!\delwaq1.*"                                !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaq\bin"
        set dest_default="!install_dir!\x64\Release\dwaq\default"
        set dest_scripts="!install_dir!\x64\Release\dwaq\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                     !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.dll"                              !dest_bin! 
        call :copyFile "!build_dir!\delwaq1\!configuration!\delwaq1.exe"                            !dest_bin! 

        call :copyFile "!checkout_src_root!\engines_gpl\waq\scripts\run_delwaq.bat"                 !dest_scripts!
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_delwaq2
rem ==========================
:delwaq2

    echo "postbuild delwaq2 . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                                        !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_plugin_wasteload\!configuration!\waq_plugin_wasteload.*"       !dest_bin!
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.*"                                   !dest_bin!
        call :copyFile "!build_dir!\delwaq2\!configuration!\delwaq2.*"                                 !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaq\bin"
        set dest_default="!install_dir!\x64\Release\dwaq\default"
        set dest_scripts="!install_dir!\x64\Release\dwaq\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                         !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\delwaq\!configuration!\delwaq.*"                                    !dest_bin! 
        call :copyFile "!build_dir!\delwaq2\!configuration!\delwaq2.*"                                  !dest_bin!

        call :copyFile "!checkout_src_root!\engines_gpl\waq\scripts\run_delwaq.bat"                     !dest_scripts! 
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_waq_run_processes
rem ==========================
:waq_run_processes

    echo "postbuild waq_run_processes . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                                     !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_run_processes\!configuration!\waq_run_processes.*"          !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaq\bin"
        set dest_default="!install_dir!\x64\Release\dwaq\default"
        set dest_scripts="!install_dir!\x64\Release\dwaq\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                         !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\waq_run_processes\!configuration!\waq_run_processes.*"              !dest_bin! 
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_duprol2delwaq
rem ==========================
:duprol2delwaq

    echo "postbuild waq_run_processes . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwaqDependentRuntimeLibraries                                             !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\duprol2delwaq\!configuration!\duprol2delwaq.*"          !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaq\bin"
        set dest_default="!install_dir!\x64\Release\dwaq\default"
        set dest_scripts="!install_dir!\x64\Release\dwaq\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwaqDependentRuntimeLibraries                                                 !dest_share!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\duprol2delwaq\!configuration!\duprol2delwaq.*"              !dest_bin! 
    )
    
goto :endproc



rem ==========================
rem === POST_BUILD_wave
rem ==========================
:wave

    echo "postbuild wave . . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwavesDependentRuntimeLibraries                                                                           !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\wave\!configuration!\wave.*"                                                          !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaves\bin"
        set dest_default="!install_dir!\x64\Release\dwaves\default"
        set dest_scripts="!install_dir!\x64\Release\dwaves\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwavesDependentRuntimeLibraries                                                                            !dest_bin!
        
        rem Temporarily rename dest_bin to share_bin to copy libraries there as well
        set dest_bin=!dest_share!
        call :copyDwaqDependentRuntimeLibraries                                                                              !dest_bin!
        set dest_bin="!install_dir!\x64\Release\dwaves\bin"
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\wave\!configuration!\wave.*"                                                           !dest_bin! 
    )
    
goto :endproc

rem ==========================
rem === POST_BUILD_wave_exe
rem ==========================
:wave_exe

    echo "postbuild wave_exe. . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyDwavesDependentRuntimeLibraries                                                                           !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\wave\!configuration!\wave_exe.*"                                                      !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dwaves\bin"
        set dest_default="!install_dir!\x64\Release\dwaves\default"
        set dest_scripts="!install_dir!\x64\Release\dwaves\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyDwavesDependentRuntimeLibraries                                                                           !dest_bin!
        
        rem Temporarily rename dest_bin to share_bin to copy libraries there as well
        set dest_bin=!dest_share!
        call :copyDwaveDependentRuntimeLibraries                                                                            !dest_bin!
        set dest_bin="!install_dir!\x64\Release\dwaves\bin"
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\wave\!configuration!\wave_exe.*"                                                      !dest_bin! 

        call :copyFile "!checkout_src_root!\engines_gpl\wave\scripts\run_dwaves.bat"                                        !dest_scripts!
    )
    
goto :endproc

rem ==========================
rem === POST_BUILD_flow2d3d
rem ==========================
:flow2d3d

    echo "postbuild flow2d3d. . ."
    
    if "%configuration%" == "Debug" (
    
        echo "Debug postbuild"
        set dest_bin="%install_dir%\x64\Debug"
        
        set dest_bin="!install_dir!\x64\Debug"
        set dest_default="!install_dir!\x64\Debug"
        set dest_scripts="!install_dir!\x64\Debug"
        set dest_plugins="!install_dir!\x64\Debug"
        set dest_share="!install_dir!\x64\Debug"
        
        call :makeDir !dest_bin!   
        call :copyFlow2D3DDependentRuntimeLibraries                                                                             !dest_bin!
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\flow2d3d\!configuration!\flow2d3d.*"                                                      !dest_bin!
    )
    
    if "%configuration%" == "Release" ( 
    
        echo "Release postbuild"

        set dest_bin="!install_dir!\x64\Release\dflow2d3d\bin"
        set dest_default="!install_dir!\x64\Release\dflow2d3d\default"
        set dest_scripts="!install_dir!\x64\Release\dflow2d3d\scripts"
        set dest_plugins="!install_dir!\x64\Release\plugins\bin"
        set dest_share="!install_dir!\x64\Release\share\bin"
        
        call :makeAllDirs   
        call :copyFlow2D3DDependentRuntimeLibraries                                                                             !dest_bin!
        
        rem Temporarily rename dest_bin to share_bin to copy libraries there as well
        set dest_bin=!dest_share!
        call :copyFlow2D3DDependentRuntimeLibraries                                                                             !dest_bin!
        set dest_bin="!install_dir!\x64\Release\dflow2d3d\bin"
        
        rem copy binaries and dll 
        call :copyFile "!build_dir!\flow2d3d\!configuration!\flow2d3d.*"                                                      !dest_bin! 

        call :copyFile "!checkout_src_root!\engines_gpl\flow2d3d\scripts\*.bat"                                                 !dest_scripts!
        call :copyFile "!checkout_src_root!\engines_gpl\flow2d3d\scripts\*.m"                                                   !dest_scripts!

        call :copyFile "!checkout_src_root!\engines_gpl\flow2d3d\default\*"                                                     !dest_default!
    )
    
goto :endproc

rem ===========================
rem === POST_BUILD_GRIDGEOM_DLL
rem ===========================
:gridgeom_dll
    echo "postbuild gridgeom_dll. . ."

    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"
    
    call :makeDir !dest_bin!
    
    call :copyFile "!build_dir!\gridgeom\!configuration!\gridgeom.dll"                                           !dest_bin!
    
    )

    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeDir !dest_share!
    call :copyFile "!build_dir!\gridgeom\!configuration!\gridgeom.dll"                                           !dest_share!
    
    )

goto :endproc


rem ==================================
rem === POST_BUILD_DFLOWFM_KERNEL_TEST
rem ==================================
:dflowfm_kernel_test
    echo "postbuild dflowfm_kernel_test . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyFile "!build_dir!\test_dflowfm_kernel\!configuration!\dflowfm_kernel_test.*"                              !dest_bin!

    )
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_bin="!install_dir!\x64\Release\tests\bin"
    set dest_default="!install_dir!\x64\Release\tests\default"
    set dest_scripts="!install_dir!\x64\Release\tests\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeAllDirs 
    call :copyFile "!build_dir!\test_dflowfm_kernel\!configuration!\dflowfm_kernel_test.*"                                !dest_bin!
    
    )

goto :endproc


rem ===============================
rem === POST_BUILD_WAQ_UTILS_F_TEST
rem ===============================
:waq_utils_f_test
    echo "postbuild waq_utils_f_test . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyFile "!build_dir!\test_waq_utils_f\!configuration!\waq_utils_f_test.*"                              !dest_bin!

    )
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_bin="!install_dir!\x64\Release\tests\bin"
    set dest_default="!install_dir!\x64\Release\tests\default"
    set dest_scripts="!install_dir!\x64\Release\tests\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeAllDirs 
    call :copyFile "!build_dir!\test_waq_utils_f\!configuration!\waq_utils_f_test.*"                                !dest_bin!
    
    )

goto :endproc



rem =============================
rem === POST_BUILD_EC_MODULE_TEST
rem =============================
:ec_module_test
    echo "postbuild ec_module_test . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyFile "!build_dir!\test_ec_module\!configuration!\ec_module_test.*"                              !dest_bin!

    )
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_bin="!install_dir!\x64\Release\tests\bin"
    set dest_default="!install_dir!\x64\Release\tests\default"
    set dest_scripts="!install_dir!\x64\Release\tests\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeAllDirs 
    call :copyFile "!build_dir!\test_ec_module\!configuration!\ec_module_test.*"                                !dest_bin!
    
    )

goto :endproc



rem ===================================
rem === POST_BUILD_TEST_DELTARES_COMMON
rem ===================================
:test_deltares_common
    echo "postbuild test_deltares_common . . ."
    
    if "%configuration%" == "Debug" (
    
    echo "Debug postbuild"
    set dest_bin="!install_dir!\x64\Debug"

    call :makeDir !dest_bin!

    call :copyFile "!build_dir!\test_deltares_common\!configuration!\test_deltares_common.*"                              !dest_bin!

    )
    
    if "%configuration%" == "Release" (
    
    echo "Release postbuild"
    
    set dest_bin="!install_dir!\x64\Release\tests\bin"
    set dest_default="!install_dir!\x64\Release\tests\default"
    set dest_scripts="!install_dir!\x64\Release\tests\scripts"
    set dest_plugins="!install_dir!\x64\Release\plugins\bin"
    set dest_share="!install_dir!\x64\Release\share\bin"
    
    call :makeAllDirs 
    call :copyFile "!build_dir!\test_deltares_common\!configuration!\test_deltares_common.*"                                !dest_bin!
    
    )

goto :endproc









:end
echo oss-post_build.cmd finished
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
