#!/bin/bash
#$ -V
#$ -j yes
#$ -cwd
    #
    # This script runs dimr on Linux
    # Adapt and use it for your own purpose
    #
    # Usage example:
    # run_dimr.sh
    #
    # adri.mourits@deltares.nl
    # 15 Feb 2017
    #

function print_usage_info {
    echo "Usage: ${0##*/} [OPTION]..."
    echo "Run a dimr model on Linux."
    echo
    echo "Options:"
    echo "-c, --corespernode <M>"
    echo "       number of cores per node, default $corespernodedefault"
    echo "-d, --debug 0xFFFFFFFF"
    echo "       maximum debug output"
    echo "-h, --help"
    echo "       print this help message and exit"
    echo "-m, --masterfile <filename>"
    echo "       dimr configuration filename, default dimr_config.xml"
    echo "    --D3D_HOME <path>"
    echo "       path to binaries and scripts"
    echo "    --NSLOTS <N>"
    echo "       number of slots=NHOSTS*CoresPerNode, default 1 (not parallel)"
    exit 1
}


# ============
# === MAIN ===
# ============

#
## Defaults
corespernodedefault=1
corespernode=1
debuglevel=0
configfile=dimr_config.xml
D3D_HOME=
runscript_extraopts=
NHOSTS=1


#
## Start processing command line options:

while [[ $# -ge 1 ]]
do
key="$1"
shift

case $key in
    -c|--corespernode)
    corespernode=$1
    shift
    ;;
    -d|--debug)
    debuglevel="$1"
    shift
    ;;
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
    --)
    echo "-- sign detected, remained options are going to be passed to dimr"
    runscript_extraopts="$runscript_extraopts $*"
    break       # exit loop, stop shifting, all remaining arguments without dashes handled below
    ;;
    -*)
    echo "option ${key} seems dedicated for dimr, therefore passing it and the following ones to the dimr"
    runscript_extraopts="$key $*"
    break       # exit loop, $key+all remaining options to dflowfm executable
    ;;
esac
done


if [ ! -f $configfile ]; then
    echo "ERROR: configfile $configfile does not exist"
    print_usage_info
fi



export OMP_NUM_THREADS=1
export NSLOTS=`expr $NHOSTS \* $corespernode` 


workdir=`pwd`



if [ -z "${D3D_HOME}" ]; then
    scriptdirname=`readlink \-f \$0`
    scriptdir=`dirname $scriptdirname`
    export D3D_HOME=$scriptdir/../..
fi
if [ ! -d $D3D_HOME ]; then
    echo "ERROR: directory $D3D_HOME does not exist"
    print_usage_info
fi
export D3D_HOME


echo "    Configfile       : $configfile"
echo "    D3D_HOME         : $D3D_HOME"
echo "    Working directory: $workdir"
echo "    Number of slots  : $NSLOTS"
echo 

    #
    # Set the directories containing the binaries
    #
export ARCH=lnx64

delwaqexedir=$D3D_HOME/$ARCH/delwaq/bin
dflowfmexedir=$D3D_HOME/$ARCH/dflowfm/lib
dimrexedir=$D3D_HOME/$ARCH/dimr/bin
esmfexedir=$D3D_HOME/$ARCH/esmf/bin
esmfbatdir=$D3D_HOME/$ARCH/esmf/scripts
flow1d2dexedir=$D3D_HOME/$ARCH/flow1d2d/bin
rrexedir=$D3D_HOME/$ARCH/rr/bin
rtcexedir=$D3D_HOME/$ARCH/rtctools/bin
swanexedir=$D3D_HOME/$ARCH/swan/bin
swanbatdir=$D3D_HOME/$ARCH/swan/scripts
thirdpartydir=$D3D_HOME/$ARCH/third_party
waveexedir=$D3D_HOME/$ARCH/wave/bin


    #
    # No adaptions needed below
    #

    # Run
export LD_LIBRARY_PATH=$dimrexedir:$dflowfmexedir:$flow1d2dexedir:$delwaqexedir:$rtcexedir:$rrexedir:$waveexedir:$swanbatdir:$swanexedir:$esmfbatdir:$esmfexedir:$thirdpartydir:$LD_LIBRARY_PATH
export PATH=$swanbatdir:$esmfbatdir:$PATH
export LD_PRELOAD=$thirdpartydir/libmkl_core.so

    echo "=== LD_LIBRARY_PATH ======================================"
    echo $LD_LIBRARY_PATH
    echo "=============================================="

    echo "=== ldd ======================================"
    ldd $rtcexedir/libRTCTools_BMI.so
    echo "=============================================="


if [ $NSLOTS -eq 1 ]; then
    echo "executing:"
    echo "$dimrexedir/dimr.exe $configfile -d $debuglevel"
    echo 
    $dimrexedir/dimr.exe $configfile -d $debuglevel
else
    #
    # Create machinefile using $PE_HOSTFILE
    if [ -n $corespernode ]; then
        if [ -e $(pwd)/machinefile ]; then
            rm -f machinefile
        fi
        for (( i = 1 ; i <= $corespernode; i++ )); do
            awk '{print $1":"1}' $PE_HOSTFILE >> $(pwd)/machinefile
        done
    else
       awk '{print $1":"2}' $PE_HOSTFILE > $(pwd)/machinefile
    fi
    echo Contents of machinefile:
    cat $(pwd)/machinefile
    echo ----------------------------------------------------------------------

    #if [ $NHOSTS -eq 1 ]; then
    #    echo "Starting mpd..."
    #    mpd &
    #fi

    node_number=$NSLOTS
    while [ $node_number -ge 1 ]; do
       node_number=`expr $node_number - 1`
       ln -s /dev/null log$node_number.irlog
    done
    
    echo "/opt/mpich2/1.4.1_intel14.0.3/bin/mpiexec -np $NSLOTS $dimrexedir/dimr.exe $configfile -d $debuglevel"
    /opt/mpich2/1.4.1_intel14.0.3/bin/mpiexec -np $NSLOTS $dimrexedir/dimr.exe $configfile -d $debuglevel
    rm -f log*.irlog
fi



    # Wait until all child processes are finished
wait

