#source /Users/hod/ROOT/root6.18/bin/thisroot.sh
#source /Users/hod/ROOT/root6.20/bin/thisroot.sh

#source /cvmfs/clicdp.cern.ch/compilers/gcc/8.1.0/x86_64-centos7/setup.sh
#source /cvmfs/clicdp.cern.ch/software/ROOT/6.12.06/x86_64-centos7-gcc8-opt/bin/thisroot.sh


export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

lsetup "root 6.20.06-x86_64-centos7-gcc8-opt"


alias python="/usr/bin/python2.7"
export LC_CTYPE="en_US.UTF-8"
