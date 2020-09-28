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
q = queue.Queue()
for signalpath,spotsizes in signals.items():
   for spotsize in spotsizes:
      relpath  = signalpath+"/"+spotsize
      fullpath = basepath_stdhep+relpath
      proc = "trident" if("trident" in relpath) else "bppp"
      command1 = "python3.7 stdhep2root.py   -p "+proc+" -d "+fullpath
      command2 = "python3.7 truthanalysis.py -p "+proc+" -d "+relpath
      # command = command1+"; "+command2
      command = "echo "+command1+"; "+command2
      print(command)
      q.put(command)

def worker():
   while True:
      item = q.get()
      # execute a task: call a shell program and wait until it completes
      command = str(item)
      # print(command)
      relpath = command.split(" ")[-1]
      p = subprocess.Popen("mkdir -p "+basepath_root+relpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      echo1 = "echo you are here: $PWD; "
      setup = 'export PYTHONPATH=$PYTHONPATH:/usr/local/anaconda3/lib/python3.7\
               export PATH=/usr/local/anaconda3/bin:$PATH\
               alias python="python3.7"\
               export STORAGEDIR=/storage/agrp/nhod/'
      echo2 = 'echo storage dir is set: $STORAGEDIR;\
               echo ROOTSYS=$ROOTSYS;\
               echo PYTHONPATH=$PYTHONPATH;\
               which python; '
      fullcommand = echo1+setup+echo2+command
      p = subprocess.Popen(fullcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      out, err = p.communicate()
      
      f = open(basepath_root+relpath+"/log.out", 'w')
      f.write(str(out))
      f.close()
      print("created log file in: "+basepath_root+relpath+"/log.out")
      f = open(basepath_root+relpath+"/log.err", 'w')
      f.write(str(err))
      f.close()
      print("created err file in: "+basepath_root+relpath+"/log.err")
      q.task_done()

cpus=multiprocessing.cpu_count() #detect number of cores
print("Creating %d threads" % cpus)
for i in range(cpus):
   t = threading.Thread(target=worker)
   t.daemon = True
   t.start()

q.join() # block until all tasks are done