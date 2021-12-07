#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import config as cfg
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom,TRandom3, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='BackgroundGenerator.py...')
parser.add_argument('-proc', metavar='proc',  required=True,  help='process')
parser.add_argument('-nevnt', metavar='events',  required=True,  help='N events')
parser.add_argument('-ntrk', metavar='background tracks',  required=True,  help='N background tracks')
parser.add_argument('-refsig', metavar='reference signal',  required=True,  help='reference signal path')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.proc
Nevt  = int(argus.nevnt)
Ntracks = int(argus.ntrk) ## tracks to generate
Refsig = argus.refsig ## tracks to generate

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


h_Px = TH1D("h_Px",";Px;Tracks",200,-0.002,0.002)   if(proc=="glaser") else TH1D("h_Px",";Px;Tracks",200,-0.006,0.006)
h_Py = TH1D("h_Py",";Py;Tracks",200,-0.002,0.002)   if(proc=="glaser") else TH1D("h_Py",";Py;Tracks",200,-0.006,0.006)
h_vx = TH1D("h_vx",";vx;Tracks",200,-0.008,0.008)   if(proc=="glaser") else TH1D("h_vx",";vx;Tracks",200,-0.006,0.006)
h_vy = TH1D("h_vy",";vy;Tracks",200,-0.0025,0.0025) if(proc=="glaser") else TH1D("h_vy",";vy;Tracks",200,-0.0015,0.0015)
h_vz = TH1D("h_vz",";vz;Tracks",200,-0.0015,0.0015) if(proc=="glaser") else TH1D("h_vz",";vz;Tracks",200,-0.0005,0.0005)
tfreal = TFile( storage+"/data/root/raw/"+Refsig+"/raw_"+proc+".root", "READ")
ttreal = tfreal.Get("tt")
print("getting events from",storage+"/data/root/raw/"+Refsig+"/raw_"+proc+".root")
m = 0
for event in ttreal:
   if(m%2==0 and m>0): print("done %g out of %g" % (m,ttreal.GetEntries()))
   for j in range(event.px.size()):
      if(event.pdgId[j]==11 and proc=="elaser"): continue
      h_Px.Fill(event.px[j])
      h_Py.Fill(event.py[j])
      h_vx.Fill(event.vx[j])
      h_vy.Fill(event.vy[j])
      h_vz.Fill(event.vz[j])
   m+=1
print("done")

tf    = TFile( storage+"/data/root/raw/flat/raw_"+proc+".root", 'recreate' )
tt    = TTree( 'tt','tt' )
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
   
   gRandom = TRandom3(0)
   gRandom.SetSeed(0)

   ### generate tracks
   ntrks = 0
   while(ntrks<Ntracks):      
      ntrks+=1
      vx0  = h_vx.GetRandom()  ## should be random
      vy0  = h_vy.GetRandom()  ## should be random
      vz0  = h_vz.GetRandom()  ## should be random
      Px   = h_Px.GetRandom()  ## should be random
      Py   = h_Py.GetRandom()  ## should be random

      # Px = -999
      # Py = -999
      # ievt = rnd.Integer(nttreal)
      # while ttreal.GetEntry(ievt):
      #    ntrkreal = rnd.Integer(ttreal.px.size())
      #    Px = ttreal.px[ntrkreal]
      #    Py = ttreal.py[ntrkreal]
      #    vx0 = ttreal.vx[ntrkreal]
      #    vy0 = ttreal.vy[ntrkreal]
      #    vz0 = ttreal.vz[ntrkreal]
      #    break

      
      
      Pz   = rnd.Uniform(Emin,Emax)
      E0   = math.sqrt(Px*Px+Py*Py+Pz*Pz+me2)
      p4   = TLorentzVector()
      p4.SetPxPyPzE(Px,Py,Pz,E0)
      wgt0 = 1
      pdgId0 = 11 if(rnd.Uniform()>0.5) else -11
      if(proc=="elaser"): pdgId0 = -11
      
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

h_Px.Write()
h_Py.Write()
h_vx.Write()
h_vy.Write()
h_vz.Write()
tf.Write()
tf.Write()
tf.Close()




