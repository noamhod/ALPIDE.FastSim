#! /bin/bash


export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS
### this sets root, but it is slow as it has a huge library path to search for
# lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
## set the views, this should be faster approach
lsetup "views LCG_101 x86_64-centos7-gcc8-opt"
### the following python installation is broken, use the default python3 from LCG setup above
#export PYTHONPATH=$PYTHONPATH:/usr/local/anaconda3/lib/python3.7
#export PATH=/usr/local/anaconda3/bin:$PATH
#alias python="python3.7"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src:$PWD/Root
export STORAGEDIR=$PWD