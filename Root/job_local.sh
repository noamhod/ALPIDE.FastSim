echo "You are here: $PWD"
cd ../
source setup2.sh
cd -
./digi -proc=bppp_bkg -evnt=$1 -seed=$2
./digi -proc=bppp     -evnt=$1 -seed=$2
./reco -proc=bppp     -evnt=$1 -seed=$2
