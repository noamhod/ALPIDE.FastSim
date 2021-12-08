#! /bin/bash


export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'
setupATLAS
lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/src:$PWD/Root
export STORAGEDIR=$PWD

