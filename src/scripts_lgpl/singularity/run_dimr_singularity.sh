#!/bin/bash

function print_usage_info {
    echo "Usage: ${0##*/} dimr_config_file [run_dimr.sh OPTIONS]..."
    echo "       ${0##*/} [-p | --parentlevel] 3 dimr_config_file [run_dimr.sh OPTIONS]..."
    echo "Runs the run_dimr.sh command of a Singularity container by wrapping and passing additional arguments."
    echo
    echo "Options:"
    echo "-h, --help"
    echo "       Print this help message and exit"
    echo "-p, --parentlevel"
    echo "       A numeric value which specifies the amount of levels to navigate to the parent directory to mount"
    echo "       (Default value: 1)"
    exit 1
}

# Variables
parent_level=1
executable=run_dimr.sh
executable_extraopts=
dimr_config_file=
container_libdir=/opt/delft3dfm_latest/lnx64/bin/ # The directory WITHIN the container that contains all the executables

# Main
if [[ $# -eq 0 ]]; then
    print_usage_info
fi

# Parse the first argument of the script
while [[ $# -ge 1 ]]
do
    key="$1"
    shift
    case $key in
         -h|--help)
        print_usage_info
        ;;
        -p|--parentlevel)
        parent_level=$1
        shift
        ;;
        *)
        dimr_config_file=$key # Assume that the first unknown argument is the dimr_config_file
        executable_extraopts=$* # Parse the remaining arguments and pass it as additional arguments to the executable as extra options
        break
        ;;
    esac
done

# Check if the dimr_config_file exists
if [ ! -z "$dimr_config_file" ] && [ ! -f "$dimr_config_file" ]; then
    echo "ERROR: Dimr config file must be set and existent."
    exit 2 # Exit code 2 to indicate that no such file is present
fi

# Retrieve the script directory
scriptdirname=`readlink \-f \$0`
scriptdir=`dirname $scriptdirname`

# Scan script directory for sif containers
shopt -s nullglob
container_file_paths=($(ls $scriptdir/*.sif))
container_file_path=

if [ ${#container_file_paths[@]} -eq 1 ]; then
    container_file_path=${container_file_paths[0]}
else
    echo "ERROR: Directory must contain one and only one *.sif file."
    exit 2 # Exit code 2 to indicate that no such file is present
fi
shopt -u nullglob

# Avoid referring to mpi libraries outside the container by feeding it a restricted PATH (and the --cleanenv flag)
restricted_path=/usr/lib64/mpich/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin 

# Set container properties
# Note that the working directory is set to a custom mounting directory
# for the container runtime environment. This mounting is to prevent 
# clashes with the internal opt directory of the container
current_working_dir=$(pwd)
mountdir=/mnt/data

working_dir=
container_working_dir=

if [[ $parent_level -ge 1 ]]; then 
    let target_workingdir_level=$parent_level+1
    working_dir=$(echo $current_working_dir | rev | cut -d'/' -f$target_workingdir_level- | rev) # Returns the desired mounting parent directory 
    container_working_dir=$mountdir/$(echo $current_working_dir | rev | cut -d'/' -f-$parent_level | rev) # Extract the directories that are traversed
elif [[ $parent_level -eq 0 ]]; then
    # Parent directory is equal to the current directory and as such there is no need to traverse the directory structure
    # while mounting
    working_dir=$current_working_dir
    container_working_dir=$mountdir
else 
    echo "Invalid parent level setting: value must be greater or equal to 0"
    exit 1
fi

echo "Executing Singularity container with:"
echo "Dimr config file                  :   $dimr_config_file"
echo "Container file                    :   $container_file_path"
echo "Current working directory         :   $current_working_dir"
echo "Mounting source directory         :   $working_dir"
echo "Mounting target directory         :   $mountdir"
echo "Container working directory       :   $container_working_dir"
echo "Executable                        :   $executable"
echo "Extra executable flags            :   $executable_extraopts"
echo "Restricted PATH                   :   $restricted_path"

singularity exec --cleanenv --bind $working_dir:$mountdir --pwd $container_working_dir --env PATH=$restricted_path $container_file_path $container_libdir/$executable -m $dimr_config_file $executable_extraopts
