#! /bin/bash

echo "You are here: $PWD"
cd .. && source setupCluster.sh && cd Root/

# ./digi -proc=glaser_bkg -evnt=$1 -seed=$2
./digi -proc=glaser -path=${STORAGEDIR}/data/root/raw/glaser -evnt=${1} -seed=${2}
./recoNew -proc=glaser     -evnt=$1 -seed=$2
