#!/bin/bash

globalErrorLevel=0

# ===============================
# === copyFile: handles error ===
# ===============================
function copyFile () {
    # handle the error
    cp -fp $1 $2
    if [ $? != 0 ]; then
        echo "can't copy \"$1\" to \"$2\"" 1>&2
        globalErrorLevel=1
    fi
}


# ===================
# === INSTALL_ALL ===
# ===================
function install_all () {
    echo "    installing all open source projects . . ."

    mkdir -p $dest_main

    d_hydro
    flow2d3d
    flow2d3d_openda
    wave
    plugin_culvert
    plugin_delftflow_traform
    datsel
    kubint
    lint
    mormerge
    vs

    return
}



# ========================
# === INSTALL_DELFT3D_FLOW
# ========================
function delft3d_flow () {
    echo "    installing delft3d-flow . . ."

    mkdir -p $dest_main

    d_hydro
    flow2d3d
    flow2d3d_openda
    plugin_culvert
    plugin_delftflow_traform
    mormerge
    vs

    return
}



# ==========================
# === INSTALL_D_HYDRO
# ==========================
function d_hydro () {
    echo "installing d_hydro . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile "bin/d_hydro.exe" "$dest_bin"
    copyFile "third_party_open/tclkit/bin/intel/deltares_hydro.exe" "$dest_bin"

    return
}



# ========================
# === INSTALL FLOW2D3D
# ========================
function flow2d3d () {
    echo "installing flow2d3d . . ."
    dest_bin="$dest_main/intel/flow/bin"
    dest_lib="$dest_main/intel/flow/lib"
    dest_default="$dest_main/intel/flow/default"

    mkdir -p $dest_bin
    mkdir -p $dest_lib
    mkdir -p $dest_default

    cp -fp lib/libflow2d3d.so                                        $dest_bin
    cp -fp lib/libflow2d3d_sp.so                                     $dest_bin
    copyFile engines_gpl/flow2d3d/scripts/meteo_old2new.m            $dest_bin
    copyFile third_party_open/DelftOnline/lib/libDelftOnline.so      $dest_lib
    copyFile third_party_open/DelftOnline/lib/libDelftOnlineJNI.so   $dest_lib
    copyFile third_party_open/DelftOnline/lib/libJavaLaunch.so       $dest_lib
    copyFile third_party_open/DelftOnline/lib/libjvm.so              $dest_lib
    copyFile bin/esm_create                                          $dest_bin
    copyFile bin/esm_delete                                          $dest_bin
    copyFile bin/esm_info                                            $dest_bin
    copyFile engines_gpl/flow2d3d/default/*.*                        $dest_default

    return
}



# ===========================
# === INSTALL FLOW2D3D_OPENDA
# ===========================
function flow2d3d_openda () {
    echo "installing flow2d3d_openda . . ."
    dest_bin="$dest_main/intel/flow/bin"
    dest_lib="$dest_main/intel/flow/lib"
    dest_default="$dest_main/intel/flow/default"

    mkdir -p $dest_bin
    mkdir -p $dest_lib
    mkdir -p $dest_default

    cp -fp lib/libflow2d3d_openda.so                                 $dest_bin
    cp -fp lib/libflow2d3d_openda_sp.so                              $dest_bin
    copyFile engines_gpl/flow2d3d/scripts/meteo_old2new.m            $dest_bin
    copyFile third_party_open/DelftOnline/lib/libDelftOnline.so      $dest_lib
    copyFile third_party_open/DelftOnline/lib/libDelftOnlineJNI.so   $dest_lib
    copyFile third_party_open/DelftOnline/lib/libJavaLaunch.so       $dest_lib
    copyFile third_party_open/DelftOnline/lib/libjvm.so              $dest_lib
    copyFile bin/esm_create                                          $dest_bin
    copyFile bin/esm_delete                                          $dest_bin
    copyFile bin/esm_info                                            $dest_bin
    copyFile engines_gpl/flow2d3d/default/*.*                        $dest_default

    return
}



# ========================
# === INSTALL WAVE
# ========================
function wave () {
    echo "installing wave . . ."

    dest_bin="$dest_main/intel/wave/bin"
    dest_default="$dest_main/intel/wave/default"
    dest_lib="$dest_main/intel/lib"

    mkdir -p $dest_bin
    mkdir -p $dest_default
    mkdir -p $dest_lib

    copyFile bin/wave.exe                                        $dest_bin
    copyFile "third_party_open/swan/bin/linux/*.*"               $dest_bin
    copyFile engines_gpl/flow2d3d/default/dioconfig.ini          $dest_default
    copyFile third_party_open/swan/scripts/swan_install.sh       $dest_bin/swan.sh

    return
}



# ==========================
# === INSTALL PLUGIN_CULVERT
# ==========================
function plugin_culvert () {
    echo "installing plugin_culvert . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile lib/libplugin_culvert.so            $dest_bin/plugin_culvert.so

    return
}



# ====================================
# === INSTALL PLUGIN_DELFTFLOW_TRAFORM
# ====================================
function plugin_delftflow_traform () {
    echo "installing plugin_delftflow_traform . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile lib/libplugin_delftflow_traform.so   $dest_bin/plugin_delftflow_traform.so

    return
}



# ==================
# === INSTALL DATSEL
# ==================
function datsel () {
    echo "installing datsel . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile bin/datsel    $dest_bin

    return
}



# ==================
# === INSTALL KUBINT
# ==================
function kubint () {
    echo "installing kubint . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile bin/kubint  $dest_bin

    return
}



# ================
# === INSTALL LINT
# ================
function lint () {
    echo "installing lint . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile bin/lint   $dest_bin

    return
}



# ====================
# === INSTALL MORMERGE
# ====================
function mormerge () {
    echo "installing mormerge . . ."
    dest_bin="$dest_main/intel/flow/bin"

    mkdir -p $dest_bin

    copyFile engines_gpl/flow2d3d/scripts/mormerge.tcl                 $dest_bin
    copyFile bin/mormerge.exe                                          $dest_bin

    return
}



# ==============
# === INSTALL VS
# ==============
function vs () {
    echo "installing vs . . ."
    dest="$dest_main/intel/util"

    mkdir -p $dest

    copyFile bin/vs  $dest

    return
}


# ============
# === MAIN ===
# ============

echo oss-install...

# Example calls:
# > install.cmd <dest directory>              # Install entire solution
# > install.cmd flow2d3d <dest directory>     # Install only project flow2d3d (and its dependencies)

# 0. defaults:
project=
dest_main=


if [ "$2" == '' ]; then
    # Install all engines, assume the first argument is a target directory

    dest_main=$1
    project=install_all
    echo Target directory: $dest_main
    echo Source          : all engines
else
    # Install the package/engine specified by the first argument. The second argument is assumed to be the target directory.

    dest_main=$2
    project=$1
    echo Target directory: $dest_main
    echo Source          : package/engine $project
fi

if [ "$dest_main" == '' ]; then
    echo "ERROR: No target directory specified as argument of oss-install.sh"
    exit 1
fi

# Change to directory tree where this batch file resides (necessary when oss-install.sh is called from outside of oss/trunk/src)
curdir=`pwd`
scriptdirname=`readlink \-f \$0`
scriptdir=`dirname $scriptdirname`
cd $scriptdir/../..

$project

cd $curdir


exit $globalErrorLevel
