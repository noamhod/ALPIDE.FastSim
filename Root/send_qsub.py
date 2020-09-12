#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call

nevents = 10

for ievnt in range(nevents):
   sievnt = str(ievnt)
   # command = "source job.sh "+str(ievnt)
   command = 'qsub job.sh -- '+sievnt+'  -q N  -o $STORAGEDIR/log_'+sievnt+'.out -e $STORAGEDIR/log_'+sievnt+'.err'
   print command
   p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
