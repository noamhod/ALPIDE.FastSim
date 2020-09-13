#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
import argparse
parser = argparse.ArgumentParser(description='send_qsub.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-c', metavar='clean',   required=True,  help='y/n')
parser.add_argument('-r', metavar='resubmit',required=False, help='"[1,2,6,99,134,...]"')
argus = parser.parse_args()
proc  = argus.p
clean = (argus.c=="y")
resubmit = (argus.r!=None)

ievents = []
if(not resubmit):
   nevents = 980
   for ievnt in range(nevents): ievents.append(ievnt)
else:
   sievents = (argus.r).replace("[","").replace("]","").split(",")
   for sievent in sievents:
      ievents.append(int(sievent))


p = subprocess.Popen("mkdir -p $STORAGEDIR/logs", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p $STORAGEDIR/data/root/dig", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p $STORAGEDIR/data/root/rec", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()

if(clean):
   p = subprocess.Popen("rm -f $STORAGEDIR/logs/*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("rm -f $STORAGEDIR/data/root/dig/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("rm -f $STORAGEDIR/data/root/rec/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()

for ievnt in ievents:
   sievnt = str(ievnt)
   command = 'qsub -F "'+sievnt+'" job.sh -q N  -o $STORAGEDIR/logs/log_'+sievnt+'.out -e $STORAGEDIR/logs/log_'+sievnt+'.err'
   print command
   p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   print("out=",out.replace("\n",""))
   if(err!=""): print("err=",err.replace("\n",""))


print("Check with:        qstat -u user")
print("List logs with:    ls -lrth $STORAGEDIR/logs/")
print("Check output with: python check_qsub.py")