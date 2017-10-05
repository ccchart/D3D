#!/bin/bash
#$ -V
#$ -j yes
#$ -cwd
    #
    # This script runs Delft3D-WAVE on Linux
    # Adapt and use it for your own purpose
    #

function print_usage_info {
    echo "Usage: ${0##*/} <mdw-file> [OPTION]..."
    echo "Run a Delft3D-WAVE model on Linux."
    echo
    echo "<mdw-file>"
    echo "        Name of mdw-file"

    echo "Options:"
    echo "-h, --help"
    echo "       print this help message and exit"
    exit 1
}


# ============
# === MAIN ===
# ============

#
## Defaults
mdwfile=
D3D_HOME=
ulimit -s unlimited


#
## Start processing command line options:

mdwfile=$1
shift
while [[ $# -ge 1 ]]
do
key="$1"
shift

case $key in
    -h|--help)
    print_usage_info
    ;;
    --D3D_HOME)
    D3D_HOME="$1"
    shift
    ;;
    --NNODES)
    NNODES="$1"
    shift
    ;;
esac
done


if [ ! -f $mdwfile ]; then
    echo "ERROR: mdwfile $mdwfile does not exist"
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

echo "    mdw-file         : $mdwfile"
echo "    D3D_HOME         : $D3D_HOME"
echo "    ARCH             : $ARCH"
echo "    Working directory: $workdir"
echo 

    #
    # Set the directories containing the binaries
    #

swanexedir=$D3D_HOME/$ARCH/swan/bin
swanbatdir=$D3D_HOME/$ARCH/swan/scripts
shareddir=$D3D_HOME/$ARCH/shared
waveexedir=$D3D_HOME/$ARCH/dwaves/bin


    #
    # No adaptions needed below
    #

    # Run
export LD_LIBRARY_PATH=$swanbatdir:$swanexedir:$waveexedir
export PATH=$swanbatdir:$PATH
    echo "executing:"
    echo "$waveexedir/wave.exe $mdwfile 0"
    echo 
$waveexedir/wave.exe $mdwfile 0 


    # Wait until all child processes are finished
wait

