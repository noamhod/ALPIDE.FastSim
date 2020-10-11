echo hostname is: $HOSTNAME
if [[ $HOSTNAME == *"wipp"* ]]; then
  cd ALPIDE.FastSim/python/
fi
echo you are here: $PWD
export PYTHONPATH=$PYTHONPATH:/usr/local/anaconda3/lib/python3.7
export PATH=/usr/local/anaconda3/bin:$PATH
export ROOTSYS=/usr/local/anaconda3/bin
ln -s /usr/local/anaconda3/pkgs/zstd-1.4.4-h3b9ef0a_1/lib/libzstd.so.1 libzstd.so.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/anaconda3/lib:$PWD
alias python="python3.7"
export STORAGEDIR=/storage/agrp/nhod/
echo storage dir is set: $STORAGEDIR
echo ROOTSYS is $ROOTSYS
echo PYTHONPATH is $PYTHONPATH
echo LD_LIBRARY_PATH is $LD_LIBRARY_PATH
which python
which root

if [[ $4 == *"y"* ]]; then
	echo submitting stdhep2root:
	python3.7 stdhep2root.py -p $1 -d $2 -g 1
else 
	echo "NOT submitting stdhep2root !!"
fi 

echo submitting truthanalysis:
python3.7 truthanalysis.py -p $1 -d $3
