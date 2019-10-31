#!/usr/bin/python
import os
import math
import subprocess
from array import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TLorentzVector
import glob
import argparse
parser = argparse.ArgumentParser(description='stdhep2root.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-e', metavar='energy', required=False,  help='beam energy')
parser.add_argument('-g', metavar='photons?', required=False,help='photons only? [default=0/1]')
argus  = parser.parse_args()
proc   = argus.p
ebeam  = argus.e
photon = (argus.g=="1")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)

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

def readparticles(name):
   particles = []
   i = 0
   fIn = open(name,'r')
   for line in fIn:
      if(line.startswith("#")): continue
      words = line.split()
      words = [word for word in words if(len(word.split())!=0)]
      if( not photon and abs(int(words[7]))!=11 ): continue
      if(     photon and abs(int(words[7]))!=22 ): continue
      particles.append( {"E":float(words[0]), "vx":float(words[1]), "vy":float(words[2]), "vz":float(words[3]), "betax":float(words[4]), "betay":float(words[5]), "betaz":float(words[6]), "pdgId":int(words[7]), "wgt":float(words[8])} )
      i += 1
   fIn.close()
   # print(particles)
   print("got %g particles" % len(particles))
   return particles



process = proc
beamenergy = ebeam
tf = TFile( '../data/root/raw_'+process+'.root', 'recreate' ) if(not photon) else TFile( '../data/root/raw_photons_'+process+'.root', 'recreate' )
tt = TTree( 'tt','tt' )
vx    = ROOT.std.vector( float )()
vy    = ROOT.std.vector( float )()
vz    = ROOT.std.vector( float )()
px    = ROOT.std.vector( float )()
py    = ROOT.std.vector( float )()
pz    = ROOT.std.vector( float )()
E     = ROOT.std.vector( float )()
pdgId = ROOT.std.vector( int )()
wgt   = ROOT.std.vector( float )()
tt.Branch('vx', vx)
tt.Branch('vy', vy)
tt.Branch('vz', vz)
tt.Branch('px', px)
tt.Branch('py', py)
tt.Branch('pz', pz)
tt.Branch('E',  E)
tt.Branch('pdgId',pdgId)
tt.Branch('wgt',wgt)

fIns = glob.glob("../data/stdhep/"+process+"/*.out")
for name in fIns:
   ### clear
   pdgId.clear()
   wgt.clear()
   vx.clear()
   vy.clear()
   vz.clear()
   px.clear()
   py.clear()
   pz.clear()
   E.clear()
   ### read file
   particles = readparticles(name)
   for ip in range(len(particles)):
      betax  = particles[ip]["betax"]
      betay  = particles[ip]["betay"]
      betaz  = particles[ip]["betaz"]
      Energy = particles[ip]["E"]
      px0,py0,pz0,E0 = elmomentum(betax,betay,betaz,Energy) if(not photon) else gmmomentum(betax,betay,betaz,Energy)
      wgt0 = particles[ip]["wgt"]
      pdgId0 = particles[ip]["pdgId"]
      vx0 = particles[ip]["vx"]*1.e-4 ## um to cm
      vy0 = particles[ip]["vy"]*1.e-4 ## um to cm
      vz0 = particles[ip]["vz"]*1.e-4 ## um to cm
      wgt.push_back(wgt0)  
      pdgId.push_back(pdgId0)  
      vx.push_back(vx0)
      vy.push_back(vy0)
      vz.push_back(vz0)
      px.push_back(px0)
      py.push_back(py0)
      pz.push_back(pz0)
      E.push_back(E0)
   # for i in range(E.size()):
      # print("i=%g, E=%g" % (i,E[i]))
   tt.Fill()
tf.Write()
tf.Write()
tf.Close()



