echo hostname is: $HOSTNAME
if [[ $HOSTNAME == *"wipp"* ]]; then
  cd ALPIDE.FastSim/python/
fi
echo you are here: $PWD
export PYTHONPATH=$PYTHONPATH:/usr/local/anaconda3/lib/python3.7
export PATH=/usr/local/anaconda3/bin:$PATH
export ROOTSYS=/usr/local/anaconda3/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/anaconda3/lib
alias python="python3.7"
export STORAGEDIR=/storage/agrp/nhod/
echo storage dir is set: $STORAGEDIR
echo ROOTSYS is $ROOTSYS
echo PYTHONPATH is $PYTHONPATH
echo LD_LIBRARY_PATH is $LD_LIBRARY_PATH
which python
which root

echo submitting stdhep2root:
python3.7 stdhep2root.py   -p $1 -d $2

echo submitting truthanalysis:
python3.7 truthanalysis.py -p $1 -d $3