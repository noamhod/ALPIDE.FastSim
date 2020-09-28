echo hostname is: $HOSTNAME
if [[ $HOSTNAME == *"weizmann"* ]]; then
  cd ALPIDE.FastSim/python/
fi
echo you are here: $PWD
export PYTHONPATH=$PYTHONPATH:/usr/local/anaconda3/lib/python3.7
export PATH=/usr/local/anaconda3/bin:$PATH
alias python="python3.7"
export STORAGEDIR=/storage/agrp/nhod/
echo storage dir is set: $STORAGEDIR
echo ROOTSYS=$ROOTSYS
echo PYTHONPATH=$PYTHONPATH
which python

echo submitting stdhep2root:
python3.7 stdhep2root.py   -p $1 -d $2

echo submitting truthanalysis:
python3.7 truthanalysis.py -p $1 -d $3
