#!/usr/bin/python
import os
import math
import subprocess
from array import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TLorentzVector
import glob
import subprocess
from subprocess import call
import argparse
parser = argparse.ArgumentParser(description='stdhep2root.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-d', metavar='fullpath',required=True,  help='the full path to the stdhep files')
parser.add_argument('-g', metavar='photons?',required=False,help='photons only? [default=0/1]')
argus  = parser.parse_args()
proc   = argus.p
path   = argus.d
photon = (argus.g=="1")
isTom  = ("ptarmigan" in path)

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")

meMeV = 0.5109989461 ## MeV
meGeV = meMeV/1000.

def beta2(betax,betay,betaz):
   return betax*betax+betay*betay+betaz*betaz

def gamma(betax,betay,betaz):
   return 1./math.sqrt(1.-beta2(betax,betay,betaz))

def elmomentum(bx,by,bz,E0):
   b2 = beta2(bx,by,bz)
   gm = gamma(bx,by,bz)
   E = gm*meGeV
   px = bx*E
   py = by*E
   pz = bz*E
   if((E0-E)/E0>0.00001): print("E=%g, E0=%g" % (E,E0))
   return px,py,pz,E

def gmmomentum(bx,by,bz,E):
   px = bx*E
   py = by*E
   pz = bz*E
   return px,py,pz,E

def readparticles(name,mpids=[]):
   print("reading: ",name)
   filtermpids = (len(mpids)>0)
   particles = {}
   i = 0
   fIn = open(name,'r')
   for line in fIn:
      if(line.startswith("#")): continue
      words = line.split()
      words = [word for word in words if(len(word.split())!=0)]
      if( not photon and abs(int(words[7]))!=11 ): continue
      # if(     photon and abs(int(words[7]))!=22 ): continue
      ## Tom: E (GeV) | x (micron) | y (micron) | z (micron) | beta_x | beta_y | beta_z | PDG_NUM | MP_Wgt | MP_ID | t (um/c) | xi
      E     = float(words[0])
      vx    = float(words[1])
      vy    = float(words[2])
      vz    = float(words[3])
      bx    = float(words[4])
      by    = float(words[5])
      bz    = float(words[6])
      pdgId = int(  words[7])
      wgt   = float(words[8])
      mpid  = int(  words[9])
      time  = float(words[10]) if(isTom) else -1
      xi    = float(words[11]) if(isTom) else -1
      if(bx*bx+by*by+bz*bz>1.): continue
      ### filter out mpids which are not in the stdhep array
      if(filtermpids and mpid not in mpids): continue
      particles.update( {mpid:{"E":E, "vx":vx, "vy":vy, "vz":vz, "betax":bx, "betay":by, "betaz":bz, "pdgId":pdgId, "wgt":wgt, "time":time, "xi":xi}} )
   fIn.close()
   # print(particles)
   print("got %g particles" % len(particles))
   return particles

# def readstdhep(name):
#    print("reading: ",name)
#    stdhep = []
#    mpids = []
#    index = 0
#    nblock = 0
#    fIn = open(name,'r')
#    for line in fIn:
#       if(line.startswith("#")): continue
#       words = line.split()
#       words = [word for word in words if(len(word.split())!=0)]
#       if(len(words)==2):
#          nblock = 0
#          continue
#       #status  pdgid mom1 mom2 child1 child2 px py pz E m vx vy vz t xi mp_weight, mp_id
#       status = int(words[0])
#       pdgid = int(words[1])
#       px = float(words[6])
#       py = float(words[7])
#       pz = float(words[8])
#       E = float(words[9])
#       m = float(words[10])
#       vx = float(words[11])
#       vy = float(words[12])
#       vz = float(words[13])
#       t = float(words[14])
#       xi = float(words[15])
#       wgt = float(words[16])
#       mpid = int(words[17])
#       mpids.append(mpid)
#       parents  = [index-nblock,index-nblock]
#       children = [index+1,index+2] if(nblock==0) else [index,index]
#       stdhep.append( {"index":index, "status":status, "pdgId":pdgid, "parents":parents, "children":children,
#                       "px":px, "py":py, "pz":pz, "E":E, "m":m, "vx":vx, "vy":vy, "vz":vz, "t":t, "xi":xi, "wgt":wgt, "mpid":mpid} )
#       nblock  += 1
#       index += 1
#    fIn.close()
#    print("got %g stdhep" % len(stdhep))
#    return stdhep,mpids
   
##############################


process = proc
targetdir = path.replace("stdhep","root/raw")
p = subprocess.Popen("mkdir -p "+targetdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()


tf = TFile( targetdir+'/raw_'+process+'.root', 'recreate' )

tt_out = TTree( 'tt','tt' )
vx_out    = ROOT.std.vector( float )()
vy_out    = ROOT.std.vector( float )()
vz_out    = ROOT.std.vector( float )()
px_out    = ROOT.std.vector( float )()
py_out    = ROOT.std.vector( float )()
pz_out    = ROOT.std.vector( float )()
E_out     = ROOT.std.vector( float )()
pdgId_out = ROOT.std.vector( int )()
mpid_out  = ROOT.std.vector( int )()
wgt_out   = ROOT.std.vector( float )()
time_out  = ROOT.std.vector( float )()
xi_out    = ROOT.std.vector( float )()
tt_out.Branch('vx', vx_out)
tt_out.Branch('vy', vy_out)
tt_out.Branch('vz', vz_out)
tt_out.Branch('px', px_out)
tt_out.Branch('py', py_out)
tt_out.Branch('pz', pz_out)
tt_out.Branch('E',  E_out)
tt_out.Branch('pdgId',pdgId_out)
tt_out.Branch('mpid',mpid_out)
tt_out.Branch('time',time_out)
tt_out.Branch('xi',xi_out)

# tt_in = TTree( 'tt_in','tt_in' )
# vx_in    = ROOT.std.vector( float )()
# vy_in    = ROOT.std.vector( float )()
# vz_in    = ROOT.std.vector( float )()
# px_in    = ROOT.std.vector( float )()
# py_in    = ROOT.std.vector( float )()
# pz_in    = ROOT.std.vector( float )()
# E_in     = ROOT.std.vector( float )()
# pdgId_in = ROOT.std.vector( int )()
# mpid_in  = ROOT.std.vector( int )()
# wgt_in   = ROOT.std.vector( float )()
# tt_in.Branch('vx', vx_in)
# tt_in.Branch('vy', vy_in)
# tt_in.Branch('vz', vz_in)
# tt_in.Branch('px', px_in)
# tt_in.Branch('py', py_in)
# tt_in.Branch('pz', pz_in)
# tt_in.Branch('E',  E_in)
# tt_in.Branch('pdgId',pdgId_in)
# tt_in.Branch('mpid',mpid_in)
# tt_in.Branch('wgt',wgt_in)

'''
tt_stdhep = TTree( 'stdhep','stdhep' )
t_stdhep        = ROOT.std.vector( float )()
vx_stdhep       = ROOT.std.vector( float )()
vy_stdhep       = ROOT.std.vector( float )()
vz_stdhep       = ROOT.std.vector( float )()
px_stdhep       = ROOT.std.vector( float )()
py_stdhep       = ROOT.std.vector( float )()
pz_stdhep       = ROOT.std.vector( float )()
E_stdhep        = ROOT.std.vector( float )()
xi_stdhep       = ROOT.std.vector( float )()
status_stdhep   = ROOT.std.vector( int )()
pdgId_stdhep    = ROOT.std.vector( int )()
mpid_stdhep     = ROOT.std.vector( int )()
wgt_stdhep      = ROOT.std.vector( float )()
parents_stdhep  = ROOT.std.vector( ROOT.std.vector( int ) )()
children_stdhep = ROOT.std.vector( ROOT.std.vector( int ) )()
tt_stdhep.Branch('t',  t_stdhep)
tt_stdhep.Branch('vx', vx_stdhep)
tt_stdhep.Branch('vy', vy_stdhep)
tt_stdhep.Branch('vz', vz_stdhep)
tt_stdhep.Branch('px', px_stdhep)
tt_stdhep.Branch('py', py_stdhep)
tt_stdhep.Branch('pz', pz_stdhep)
tt_stdhep.Branch('E',  E_stdhep)
tt_stdhep.Branch('xi', xi_stdhep)
tt_stdhep.Branch('status',status_stdhep)
tt_stdhep.Branch('pdgId',pdgId_stdhep)
tt_stdhep.Branch('mpid',mpid_stdhep)
tt_stdhep.Branch('wgt',wgt_stdhep)
tt_stdhep.Branch('parents',parents_stdhep)
tt_stdhep.Branch('children',children_stdhep)
'''


fIns = glob.glob(path+"/*.out")
for name in fIns:
   ### clear output tree branches
   mpid_out.clear()
   pdgId_out.clear()
   wgt_out.clear()
   vx_out.clear()
   vy_out.clear()
   vz_out.clear()
   px_out.clear()
   py_out.clear()
   pz_out.clear()
   E_out.clear()
   xi_out.clear()
   time_out.clear()
   
   '''
   ### clear input tree branches
   mpid_in.clear()
   pdgId_in.clear()
   wgt_in.clear()
   vx_in.clear()
   vy_in.clear()
   vz_in.clear()
   px_in.clear()
   py_in.clear()
   pz_in.clear()
   E_in.clear()

   ### clear stdhep tree branches
   t_stdhep.clear()
   vx_stdhep.clear()
   vy_stdhep.clear()
   vz_stdhep.clear()
   px_stdhep.clear()
   py_stdhep.clear()
   pz_stdhep.clear()
   E_stdhep.clear()
   xi_stdhep.clear()
   wgt_stdhep.clear()
   status_stdhep.clear()
   pdgId_stdhep.clear()
   mpid_stdhep.clear()
   for n in range(parents_stdhep.size()): parents_stdhep[n].clear()
   parents_stdhep.clear()
   for n in range(children_stdhep.size()): children_stdhep[n].clear()
   children_stdhep.clear()
   '''

   ### read files
   particles    = readparticles(name)
   '''
   stdhep,mpids = readstdhep(name.replace("particles.out","events.stdhep"))
   inparticles  = readparticles(name.replace(".out",".in"),mpids)
   # print(stdhep)
   '''

   ### loop over all output particles
   for MP_ID,particle in particles.items():
      betax  = particle["betax"]
      betay  = particle["betay"]
      betaz  = particle["betaz"]
      Energy = particle["E"]
      px0,py0,pz0,E0 = elmomentum(betax,betay,betaz,Energy) if(not photon) else gmmomentum(betax,betay,betaz,Energy)
      wgt0 = particle["wgt"]
      pdgId0 = particle["pdgId"]
      vx0 = particle["vx"]*1.e-4 ## um to cm
      vy0 = particle["vy"]*1.e-4 ## um to cm
      vz0 = particle["vz"]*1.e-4 ## um to cm
      xi0 = particle["xi"]
      t0  = particle["time"]
      mpid_out.push_back(MP_ID)
      wgt_out.push_back(wgt0)  
      pdgId_out.push_back(pdgId0)  
      vx_out.push_back(vx0)
      vy_out.push_back(vy0)
      vz_out.push_back(vz0)
      px_out.push_back(px0)
      py_out.push_back(py0)
      pz_out.push_back(pz0)
      E_out.push_back(E0)
      time_out.push_back(t0)
      xi_out.push_back(xi0)
   tt_out.Fill()
   
   '''
   ### loop over all input particles
   for MP_ID,particle in inparticles.items():
      betax  = particle["betax"]
      betay  = particle["betay"]
      betaz  = particle["betaz"]
      Energy = particle["E"]
      px0,py0,pz0,E0 = elmomentum(betax,betay,betaz,Energy) if(not photon) else gmmomentum(betax,betay,betaz,Energy)
      wgt0 = particle["wgt"]
      pdgId0 = particle["pdgId"]
      vx0 = particle["vx"]*1.e-4 ## um to cm
      vy0 = particle["vy"]*1.e-4 ## um to cm
      vz0 = particle["vz"]*1.e-4 ## um to cm
      mpid_in.push_back(MP_ID)
      wgt_in.push_back(wgt0)  
      pdgId_in.push_back(pdgId0)  
      vx_in.push_back(vx0)
      vy_in.push_back(vy0)
      vz_in.push_back(vz0)
      px_in.push_back(px0)
      py_in.push_back(py0)
      pz_in.push_back(pz0)
      E_in.push_back(E0)
   tt_in.Fill()
   '''
   
   '''
   ### loop over all output particles
   for stdhep in stdhep:
      status0 = stdhep["status"]
      wgt0 = stdhep["wgt"]
      pdgId0 = stdhep["pdgId"]
      mpid0 = stdhep["mpid"]
      xi0 = stdhep["xi"]
      vx0 = stdhep["vx"]*1.e-4 ## um to cm
      vy0 = stdhep["vy"]*1.e-4 ## um to cm
      vz0 = stdhep["vz"]*1.e-4 ## um to cm
      px0 = stdhep["px"]
      py0 = stdhep["py"]
      pz0 = stdhep["pz"]
      E0 = stdhep["E"]
      t0 = stdhep["t"]
      status_stdhep.push_back(status0)
      mpid_stdhep.push_back(mpid0)
      xi_stdhep.push_back(xi0)
      wgt_stdhep.push_back(wgt0)  
      pdgId_stdhep.push_back(pdgId0)  
      vx_stdhep.push_back(vx0)
      vy_stdhep.push_back(vy0)
      vz_stdhep.push_back(vz0)
      px_stdhep.push_back(px0)
      py_stdhep.push_back(py0)
      pz_stdhep.push_back(pz0)
      E_stdhep.push_back(E0)
      t_stdhep.push_back(t0)
      parents_stdhep.push_back( ROOT.std.vector( int )() )
      children_stdhep.push_back( ROOT.std.vector( int )() )
      nparticles = mpid_stdhep.size()-1
      for parent in stdhep["parents"]:  parents_stdhep[nparticles].push_back(parent)
      for child  in stdhep["children"]: children_stdhep[nparticles].push_back(child)
   tt_stdhep.Fill()
   '''

tt_out.Write()
'''
tt_in.Write()
tt_stdhep.Write()
'''
tf.Write()
tf.Write()
tf.Close()
