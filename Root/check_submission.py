#!/usr/bin/env python
import os, sys
import subprocess
from subprocess import call

nevents = 980

def FormatEventID(evnt):
   sevnt = ""
   if(evnt<10):                         sevnt = "000000"+str(evnt)
   if(evnt>=10 and evnt<100):           sevnt = "00000"+str(evnt)
   if(evnt>=100 and evnt<1000):         sevnt = "0000"+str(evnt)
   if(evnt>=1000 and evnt<10000):       sevnt = "000"+str(evnt)
   if(evnt>=10000 and evnt<100000):     sevnt = "00"+str(evnt)
   if(evnt>=100000 and evnt<1000000):   sevnt = "0"+str(evnt)
   if(evnt>=1000000 and evnt<10000000): sevnt = str(evnt) ## assume no more than 9,999,999 events...
   return sevnt

problems = []
for ievnt in range(nevents):
   sievnt = FormatEventID(ievnt)
   path = os.path.expandvars('$STORAGEDIR/data/root/rec/')
   fname = path+'rec_bppp_'+sievnt+'.root'
   if(os.path.isfile(fname)):
      if(os.stat(fname).st_size<10000):
         print("ERROR - file too small:",fname)
         problems.append(ievnt)
   else:
      print("ERROR - file do not exist:",fname)
      problems.append(ievnt)

if(len(problems)>0):
   print("----------\nPlease resubmit with this list:")
   print(problems)
else:
   print("No problems found!")
   print("Can proceed to merge the output as:")
   print('   hadd -f $STORAGEDIR/data/root/rec_bppp.root  $STORAGEDIR/data/root/rec/rec_bppp_*.root')
