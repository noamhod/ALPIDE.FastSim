#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
from random import seed, randint
import argparse
parser = argparse.ArgumentParser(description='send_qsub.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-c', metavar='clean',   required=True,  help='y/n')
parser.add_argument('-n', metavar='nevents', required=True,  help='980')
parser.add_argument('-r', metavar='resubmit',required=False, help='"[1,2,6,99,134,...]"')
argus = parser.parse_args()
proc  = argus.p
clean = (argus.c=="y")
nevents = int(argus.n)
resubmit = (argus.r!=None)


## build the list from scratch / resubmission subset
ievents = []
if(not resubmit):
   for ievnt in range(nevents): ievents.append(ievnt)
else:
   sievents = (argus.r).replace("[","").replace("]","").split(",")
   for sievent in sievents:
      ievents.append(int(sievent))


## make the directories if these do not exist
p = subprocess.Popen("mkdir -p $STORAGEDIR/logs", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p $STORAGEDIR/data/root/dig", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p $STORAGEDIR/data/root/rec", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()


## clean all dirs?
if(clean):
   reallyclean = raw_input("Are you sure you really want to clean ??? [yes/no]")
   print("OK, your choice is: "+reallyclean)
   if(reallyclean!="yes" and reallyclean!="no"):
      print("cannot understand your choice: "+reallyclean)
      quit()
   if(reallyclean=="yes"):
      p = subprocess.Popen("rm -f $STORAGEDIR/logs/*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      p = subprocess.Popen("rm -f $STORAGEDIR/data/root/dig/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      p = subprocess.Popen("rm -f $STORAGEDIR/data/root/rec/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()


## send a pilot job to make all the dictionaries
if(not os.path.exists('AutoDict_vector_vector_int____cxx.so')):
   print('Dictionaries should be generated on the first call')
   print('Please first manually execute: ./job_local.sh 0 0')
   print('and then re-run this script as usual')
   quit()

# seed random number generator
seed(1)   
   
## run!
for ievnt in ievents:
   sievnt = str(ievnt)
   randseed = str(randint(0,1000000))
   command = 'qsub -F "'+sievnt+' '+randseed+'" job_qsub.sh -q N  -o $STORAGEDIR/logs/log_'+sievnt+'.out -e $STORAGEDIR/logs/log_'+sievnt+'.err'
   print command
   p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   print("out=",out.replace("\n",""))
   if(err!=""): print("err=",err.replace("\n",""))


print("Check with:        qstat -u user")
print("List logs with:    ls -lrth $STORAGEDIR/logs/")
print("Check output with: python check_submission.py")
