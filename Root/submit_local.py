#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
import Queue
import threading
import multiprocessing
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
   p = subprocess.Popen("rm -f $STORAGEDIR/logs/*", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("rm -f $STORAGEDIR/data/root/dig/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()
   p = subprocess.Popen("rm -f $STORAGEDIR/data/root/rec/*.root", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   out, err = p.communicate()


## send a pilot job to make all the dictionaries
p = subprocess.Popen('./job_local.sh 0', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()


## run!
q = Queue.Queue()
for ievnt in ievents:
   sievnt = str(ievnt)
   command = './job_local.sh '+sievnt
   print(command)
   q.put(command)

def worker():
   while True:
      item = q.get()
      # execute a task: call a shell program and wait until it completes
      command = str(item)
      evntid = command.split(" ")[1]
      evntid.split
      path = os.path.expandvars('$STORAGEDIR/logs')
      p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      f = open(path+"/log_"+evntid+".out", 'w')
      f.write(out)
      f.close()
      f = open(path+"/log_"+evntid+".err", 'w')
      f.write(err)
      f.close()
      q.task_done()

cpus=multiprocessing.cpu_count() #detect number of cores
print("Creating %d threads" % cpus)
for i in range(cpus):
   t = threading.Thread(target=worker)
   t.daemon = True
   t.start()

q.join() # block until all tasks are done




print("List logs with:    ls -lrth $STORAGEDIR/logs/")
print("Check output with: python check_submission.py")