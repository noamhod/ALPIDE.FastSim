#! /bin/bash
#PBS -m n
#PBS -l walltime=12:00:00,vmem=9gb

### job arguments
eventNeeded=${parname1}
seedSet=${parname2}
processName=${parname3}
echo "You are here: $PWD"
export BASEPATH=/storage/agrp/arkas/LUXETrackerFastSim/ALPIDE.FastSim/Root
cd $BASEPATH
echo "Now you are here: $PWD"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS
cd $ATLAS_LOCAL_ROOT_BASE/user
echo "Now you are here: $PWD"
source atlasLocalSetup.sh
### this sets root, but it is slow as it has a huge library path to search for
# lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
## set the views, this should be faster approach
lsetup "views LCG_101 x86_64-centos7-gcc8-opt"
export STORAGEDIR=/storage/agrp/arkas/LUXETrackerFastSim/ALPIDE.FastSim
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$STORAGEDIR/src:$STORAGEDIR/Root
echo "STORAGEDIR=: $STORAGEDIR"
cd $BASEPATH
echo "Now you are here: $PWD"
ls -lrth
#./digi -proc=glaser_bkg -evnt=$1 -seed=$2
./digi -proc=${processName} -path=${STORAGEDIR}/data/root/raw/${processName} -evnt=${eventNeeded} -seed=${seedSet}
./recoNew -proc=${processName} -dobg=0 -evnt=${eventNeeded} -seed=${seedSet}
