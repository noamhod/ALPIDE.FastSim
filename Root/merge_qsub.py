#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call

nevents = 980

for ievnt in range(nevents):
   sievnt = str(ievnt)
   fOname = '$STORAGEDIR/logs/log_'+sievnt+'.out'
   fEname = '$STORAGEDIR/logs/log_'+sievnt+'.err'
   with open(fOname,'r') as f:
      lines = f.read().splitlines()
      last_line = lines[-1]
      if("Done 0 out of 1 --> CPUav=" not in last_line):
         print("ERROR in:",fname)
         quit()

command = 'hadd -f $STORAGEDIR/rec/rec_bppp.root  $STORAGEDIR/rec/*.root'
print(command)
#p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#out, err = p.communicate()
#print "out=",out
#if(err!=""): print "err=",err

print("Check with: `qstat -u user`")
print("List logs: `ls -lrth $STORAGEDIR/logs/`")
