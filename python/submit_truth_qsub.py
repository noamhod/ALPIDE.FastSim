#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
import queue
import threading
import multiprocessing
from random import seed, randint
import argparse
parser = argparse.ArgumentParser(description='submit_trut_qsub.py...')
parser.add_argument('-s', metavar='do stdhep2root step?', required=True,  help='do stdhep2root step? [y/n]')
argus = parser.parse_args()
dostdhep = argus.s

basepath_stdhep = os.path.expandvars('$STORAGEDIR/data/stdhep/')
basepath_root   = os.path.expandvars('$STORAGEDIR/output/root/raw/')

jeti40elaser = "trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/"
jeti40glaser = "bppp/IPstrong_V1.1.00/JETI40/g_laser/16.5GeV/"
phase2elaser = "trident/IPstrong_V1.1.00/phaseII/e_laser/16.5GeV/"
phase2glaser = "bppp/IPstrong_V1.1.00/phaseII/g_laser/16.5GeV/"

signals = {}
#signals.update( {jeti40elaser : ["w0_3000nm","w0_3500nm","w0_4000nm","w0_4500nm","w0_5000nm","w0_6500nm","w0_8000nm","w0_10000nm","w0_13000nm","w0_15000nm","w0_20000nm","w0_50000nm","w0_100000nm"]} )
signals.update( {jeti40glaser : ["w0_3000nm","w0_3500nm","w0_4000nm","w0_4500nm","w0_5000nm","w0_6500nm","w0_8000nm"]} )
#signals.update( {phase2elaser : ["w0_11000nm", "w0_12000nm", "w0_16000nm", "w0_20000nm"]} )
signals.update( {phase2glaser : ["w0_9000nm", "w0_10000nm", "w0_11000nm", "w0_12000nm","w0_16000nm","w0_20000nm"]} )

## run!
for signalpath,spotsizes in signals.items():
   for spotsize in spotsizes:
      relpath  = signalpath+"/"+spotsize
      fullpath = basepath_stdhep+relpath
      logname = relpath.replace("/","_").replace("__","_")
      proc = "trident" if("trident" in relpath) else "bppp"
      p = subprocess.Popen("mkdir -p "+basepath_root+relpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      
      qsub = 'qsub -F "'+proc+' '+fullpath+' '+relpath+' '+dostdhep+'" job_truth.sh -q N  -o $STORAGEDIR/logs/log_'+logname+'.out -e $STORAGEDIR/logs/log_'+logname+'.err'
      print(qsub)
      p = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      print("out=",str(out).replace("\n",""))
      if(err!=""): print("err=",str(err).replace("\n",""))

print("Check with:        qstat -u nhod")
print("List logs with:    ls -lrth $STORAGEDIR/logs/")
print("kill all jobs:     qselect -u nhod | xargs qdel")
