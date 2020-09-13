#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call

#ievents = []
#nevents = 980
#for ievnt in range(nevents): ievents.append(ievnt)

ievents = [99, 199, 299, 399, 499, 599, 699, 799, 899]

for ievnt in ievents:
   sievnt = str(ievnt)
   p = subprocess.Popen("mkdir -p $STORAGEDIR/logs", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("rm -f $STORAGEDIR/logs/*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   command = 'qsub -F "'+sievnt+'" job.sh -q N  -o $STORAGEDIR/logs/log_'+sievnt+'.out -e $STORAGEDIR/logs/log_'+sievnt+'.err'
   print command
   p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   print "out=",out.replace("\n","")
   if(err!=""): print "err=",err.replace("\n","")

print "Check with: `qstat -u user`"
print "List logs: `ls -lrth $STORAGEDIR/logs/`"
