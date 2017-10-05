#!/bin/bash
#$ -V
#$ -j yes
#$ -cwd
    #
    # This script runs Delft3D-FLOW on Linux
    # Adapt and use it for your own purpose
    #

function print_usage_info {
    echo "Usage: ${0##*/} [OPTION]..."
    echo "Run a Delft3D-FLOW model on Linux."
    echo
    echo "Options:"
    echo "-h, --help"
    echo "       print this help message and exit"
    echo "<filename>"
    echo "       Delft3D-FLOW configuration filename, default config_d_hydro.xml"
    exit 1
}


# ============
# === MAIN ===
# ============

#
## Defaults
configfile=config_d_hydro.xml
D3D_HOME=
ulimit -s unlimited


#
## Start processing command line options:

while [[ $# -ge 1 ]]
do
key="$1"
shift

case $key in
    -h|--help)
    print_usage_info
    ;;
    -m|--masterfile)
    configfile="$1"
    shift
    ;;
    --D3D_HOME)
    D3D_HOME="$1"
    shift
    ;;
	--NNODES)
    NNODES="$1"
    shift
    ;;
    *)
    configfile="$*"
    break
    ;;
esac
done


if [ ! -f $configfile ]; then
    echo "ERROR: configfile $configfile does not exist"
    print_usage_info
fi


workdir=`pwd`

if [ -z "${D3D_HOME}" ]; then
    scriptdirname=`readlink \-f \$0`
    scriptdir=`dirname $scriptdirname`
    D3D_HOME=$scriptdir/../..
else
    # D3D_HOME is passed through via argument --D3D_HOME
    # Commonly its value is "/some/path/lnx64/scripts/../.."
    # Remove "/../.." at the end of the string
    scriptdir=${D3D_HOME%"/../.."}
fi
if [ ! -d $D3D_HOME ]; then
    echo "ERROR: directory $D3D_HOME does not exist"
    print_usage_info
fi
export D3D_HOME
 
    # find ARCH from scriptdir path
pth=( $( echo $scriptdir | tr "/" "\n" ) )
a=${#pth[@]}-2
export ARCH=${pth[a]}

echo "    Configfile       : $configfile"
echo "    D3D_HOME         : $D3D_HOME"
echo "    ARCH             : $ARCH"
echo "    Working directory: $workdir"
echo 

    #
    # Set the directories containing the binaries
    #

flow2d3dexedir=$D3D_HOME/$ARCH/dflow2d3d/bin
shareddir=$D3D_HOME/$ARCH/shared


    #
    # No adaptions needed below
    #

    # Run
export LD_LIBRARY_PATH=$shareddir:$flow2d3dexedir


    echo "executing:"
    echo "$flow2d3dexedir/d_hydro.exe $configfile"
    echo 
    $flow2d3dexedir/d_hydro.exe $configfile



    # Wait until all child processes are finished
wait

