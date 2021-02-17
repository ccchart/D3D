#!/bin/bash

# set NHOSTS manually here:
NHOSTS=$(cat ../../../nhosts)
export NSLOTS=$NHOSTS




# user defined settings
# export windownsNode=CHATWDELT01
export FF2NFdir=coupleFiles

#fssid=$(cat ../../../fssid)

# ------------- Processing -----------------------
# export paths
export rundir=$(pwd)
export reldir=${PWD##*/}

# copy COSUMO_template_settings.xml to COSUMOsettings.xml
cp COSUMO_template_settings.xml COSUMOsettings.xml

# replace keywords COSUMOsettings.xml
export find1=%RUNFOLDERNAME%
export replace1=$rundir
# export find2=%WINDOWSNODE%
# export replace2=$windownsNode
export find3=%FF2NFDIR%
export replace3=$FF2NFdir

sed -i 's='$find1'='$rundir'=g' COSUMOsettings.xml
# sed -i "s/$find2/$replace2/" COSUMOsettings.xml
sed -i "s/$find3/$replace3/" COSUMOsettings.xml


# mount directory
#/opt/cjlanham-wa/mount_smb.sh $windownsNode
#mkdir -p /apps/deltares/nfsShare/CosumoShare/$fssid
#cp -rp $rundir/* /apps/deltares/nfsShare/CosumoShare/$fssid



    #
    # Specify the config file to be used here
argfile=config_d_hydro.xml


    #
    # Set the directory containing d_hydro.exe here
export ARCH=lnx64
export D3D_HOME=../bin
exedir=$D3D_HOME/$ARCH/flow2d3d/bin

 
    #
    # No adaptions needed below
    #

    #
    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:/usr/lib64/mpich-3.2/lib:$LD_LIBRARY_PATH 
export PATH="/usr/lib64/mpich-3.2/bin:${PATH}"
    ### limit MpiCH to running on ports 50100 - 50500
export MPIEXEC_PORT_RANGE=50100:50500
 
    ### The file "machinefile" is assumed to be created manually and should
    ### look like:
    ### machinename1:4
    ### machinename2:4
    ### machinename3:4

echo Contents of machinefile:
cat $(pwd)/../../../machinefile
echo ----------------------------------------------------------------------






# link mpich debug rubbish to /dev/null
node_number=$NSLOTS
while [ $node_number -ge 1 ]; do
   node_number=`expr $node_number - 1`
   ln -s /dev/null log$node_number.irlog
done


    #
    # Run
/usr/lib64/mpich-3.2/bin/mpirun -f $(pwd)/../../../machinefile -np $NSLOTS $exedir/d_hydro.exe $argfile 1>screen.log 2>&1

    #
    # Copy Cosumo/BFSCH files to the working directory, replacing existing files  when they already exist
#cp -rfp /apps/deltares/nfsShare/CosumoShare/COSUMO/ . 


    #
    # Clean up, finish MPICH2 network
rm -f log*.irlog
#rm -rf /apps/deltares/nfsShare/CosumoShare/$fssid

