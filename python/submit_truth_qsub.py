#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
import queue
import threading
import multiprocessing
from random import seed, randint

basepath_stdhep = os.path.expandvars('$STORAGEDIR/data/stdhep/')
basepath_root   = os.path.expandvars('$STORAGEDIR/output/root/raw/')

signals = {}
signals.update( {"trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/" :["w0_3000nm","w0_3500nm","w0_4000nm","w0_4500nm","w0_5000nm","w0_8000nm","w0_20000nm","w0_50000nm","w0_100000nm"]} )
signals.update( {"bppp/IPstrong_V1.1.00/JETI40/g_laser/16.5GeV/"    :["w0_3000nm","w0_3500nm","w0_4000nm","w0_4500nm","w0_5000nm","w0_8000nm"]} )
signals.update( {"trident/IPstrong_V1.1.00/phaseII/e_laser/16.5GeV/":["w0_8000nm","w0_9000nm","w0_10000nm","w0_11000nm","w0_12000nm"]} )
signals.update( {"bppp/IPstrong_V1.1.00/phaseII/g_laser/16.5GeV/"   :["w0_8000nm","w0_9000nm","w0_10000nm","w0_11000nm","w0_12000nm"]} )

## run!
for signalpath,spotsizes in signals.items():
   for spotsize in spotsizes:
      relpath  = signalpath+"/"+spotsize
      fullpath = basepath_stdhep+relpath
      logname = relpath.replace("/","_")
      proc = "trident" if("trident" in relpath) else "bppp"
      p = subprocess.Popen("mkdir -p "+basepath_root+relpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      
      qsub = 'qsub -F "'+proc+' '+fullpath+' '+relpath+'" job_truth.sh -q N  -o $STORAGEDIR/logs/log_'+logname+'.out -e $STORAGEDIR/logs/log_'+logname+'.err'
      print(qsub)
      p = subprocess.Popen(qsub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      print("out=",out.replace("\n",""))
      if(err!=""): print("err=",err.replace("\n",""))

print("Check with:        qstat -u user")
print("List logs with:    ls -lrth $STORAGEDIR/logs/")
print("Check output with: python check_submission.py")