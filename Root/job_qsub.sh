echo "You are here: $PWD"
export BASEBPTH=/srv01/agrp/nhod/ALPIDE.FastSim/Root
cd $BASEBPTH
echo "Now you are here: $PWD"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
cd $ATLAS_LOCAL_ROOT_BASE/user
echo "Now you are here: $PWD"
source atlasLocalSetup.sh
setupATLAS
lsetup "root 6.16.00-x86_64-centos7-gcc8-opt"
export STORAGEDIR=/storage/agrp/nhod/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
echo "STORAGEDIR=: $STORAGEDIR"
cd $BASEBPTH
echo "Now you are here: $PWD"
ls -lrth
./digi -proc=bppp_bkg -evnt=$1
./digi -proc=bppp     -evnt=$1
./reco -proc=bppp     -evnt=$1
