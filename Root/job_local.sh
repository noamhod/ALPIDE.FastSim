echo "You are here: $PWD"
cd ../
source setupROOT6.binaries.python2.sh
cd -
./digi -proc=bppp_bkg -evnt=$1
./digi -proc=bppp     -evnt=$1
./reco -proc=bppp     -evnt=$1
