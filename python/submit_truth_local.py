#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call
import queue
import threading
import multiprocessing
import argparse
from random import seed, randint
parser = argparse.ArgumentParser(description='submit_local.py...')
parser.add_argument('-p', metavar='prefix', required=True,  help='e.g. trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/')
argus = parser.parse_args()
prefix  = argus.p

spotsizes = ["w0_3000nm","w0_3500nm","w0_4000nm","w0_4500nm","w0_5000nm","w0_8000nm","w0_20000nm","w0_50000nm","w0_100000nm"]

## run!
q = queue.Queue()
for spotsize in spotsizes:
   relpath = prefix+"/"+spotsize
   proc = "trident" if("trident" in relpath) else "bppp"
   command = "/usr/local/bin/python3 truthanalysis.py -p "+proc+" -d "+relpath
   print(command)
   q.put(command)

def worker():
   while True:
      item = q.get()
      # execute a task: call a shell program and wait until it completes
      command = str(item)
      # print(command)
      relpath = command.split(" ")[5]
      basepath = os.path.expandvars('$STORAGEDIR/output/root/raw/')
      p = subprocess.Popen("mkdir -p "+basepath+relpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      echo1 = "echo you are here: $PWD; "
      setup = 'cd /usr/local/bin/; source thisroot.sh; cd -;\
               export PYTHONPATH=$PYTHONPATH:/usr/local/Cellar/root/6.18.04_1/lib/root;\
               export STORAGEDIR=/Users/hod/GitHub/ALPIDE.FastSim;'
      echo2 = 'echo storage dir is set: $STORAGEDIR;\
               echo ROOTSYS=$ROOTSYS;\
               echo PYTHONPATH=$PYTHONPATH;\
               which python; '
      fullcommand = echo1+setup+echo2+command
      p = subprocess.Popen(fullcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      
      f = open(basepath+relpath+"/log.out", 'w')
      f.write(str(out))
      f.close()
      print("created log file in: "+basepath+relpath+"/log.out")
      f = open(basepath+relpath+"/log.err", 'w')
      f.write(str(err))
      f.close()
      print("created err file in: "+basepath+relpath+"/log.err")
      q.task_done()

cpus=multiprocessing.cpu_count() #detect number of cores
print("Creating %d threads" % cpus)
for i in range(cpus):
   t = threading.Thread(target=worker)
   t.daemon = True
   t.start()

q.join() # block until all tasks are done