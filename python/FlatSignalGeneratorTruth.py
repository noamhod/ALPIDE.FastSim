#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import config as cfg
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='BackgroundGenerator.py...')
parser.add_argument('-proc', metavar='proc',  required=True,  help='process')
parser.add_argument('-nevnt', metavar='events',  required=True,  help='N events')
parser.add_argument('-ntrk', metavar='background tracks',  required=True,  help='N background tracks')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.proc
Nevt  = int(argus.nevnt)
Ntracks = int(argus.ntrk) ## tracks to generate

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
storage = os.path.expandvars("$STORAGEDIR")

#############################################

### electron mass:
me = 0.51099895/1000. ### GeV
me2 = me*me
cm2m = 1.e-2
cm2um = 1.e4
um2cm = 1.e-4

Emin = 1
Emax = 16.5

tf = TFile( storage+"/data/root/raw/flat/raw_"+proc+".root", 'recreate' )
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

n=0 ### init n
for n in range(Nevt):
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
   
   rnd = TRandom()
   rnd.SetSeed()
   
   ### generate tracks
   ntrks = 0
   while(ntrks<Ntracks):      
      ntrks+=1
      vx0 = 0 ## should be random
      vy0 = 0 ## should be random
      vz0 = 0 ## should be random
      Px = 0  ## should be random
      Py = 0  ## should be random
      Pz = rnd.Uniform(Emin,Emax)
      E0  = math.sqrt(Px*Px+Py*Py+Pz*Pz+me2)
      p4 = TLorentzVector()
      p4.SetPxPyPzE(Px,Py,Pz,E0)
      wgt0 = 1
      pdgId0 = 11 if(rnd.Uniform()>0.5) else -11
      
      ### fill output vectors
      wgt.push_back(wgt0)  
      pdgId.push_back(pdgId0)  
      vx.push_back(vx0)
      vy.push_back(vy0)
      vz.push_back(vz0)
      px.push_back(p4.Px())
      py.push_back(p4.Py())
      pz.push_back(p4.Pz())
      E.push_back(p4.E())
      
   if(n%100==0): print("done %g out of %g" % (n,Nevt))
   tt.Fill()
   
tf.Write()
tf.Write()
tf.Close()




