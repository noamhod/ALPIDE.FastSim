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
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-nevents', metavar='events',  required=True,  help='N events')
parser.add_argument('-nbckgrtrk', metavar='background tracks',  required=True,  help='N background tracks')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.p
Nevt  = int(argus.nevents)
### background stuff
NbkgTracks = int(argus.nbckgrtrk) ## tracks to generate
sides = "e+" if(proc=="trident") else "e+e-" ## jsut the default
print("Running with proc=%s and sides=%s" % (proc,sides))

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage = os.path.expandvars("$STORAGEDIR")

#############################################
### read configuration from file
cfgmap = cfg.set(proc,sides,True)

#############################################
NeventsToDraw = 10

### electron mass:
me = 0.51099895/1000. ### GeV
me2 = me*me
cm2m = 1.e-2
cm2um = 1.e4
um2cm = 1.e-4

### magnetic field
B  = 1.0 if(proc=="trident") else 1.7 # Tesla
LB = 1.029   # meters

### possible energies 
Emax = 17.5 if(proc=="trident") else 16 # GeV
Emin = 1.00 if(proc=="trident") else 2 # GeV

### geometry:
zDipoleExit = 202.9
xDipoleExitMinAbs = 1.5 if(proc=="bppp") else 4   ## cm --> TODO: need tuning
xDipoleExitMaxAbs = 25  if(proc=="bppp") else 30  ## cm --> TODO: need tuning
yDipoleExitMin = -0.05 ## cm --> TODO: need tuning
yDipoleExitMax = +0.05 ## cm --> TODO: need tuning
xAbsMargins = 0.025 # cm --> TODO: need tuning
yAbsMargins = 0.025 if(proc=="bppp") else 0.1 # cm --> TODO: need tuning

### stave geometry
Hstave    = 1.5  # cm
Lstave    = 50 #27.12 if(proc=="bppp") else 50 # cm
Rbeampipe = 2.413 #4 # cm
RoffsetBfield = 5.7 if(proc=="bppp") else 14 # cm
xPsideL = -RoffsetBfield-Lstave
xPsideR = -RoffsetBfield       
xEsideL = +RoffsetBfield       
xEsideR = +RoffsetBfield+Lstave
yUp = +Hstave/2.
yDn = -Hstave/2.

#############################################
def getgeometry(dipole=False):
   tfile = TFile(storage+"/data/root/"+proc+"_geometry.root","READ")
   geometry = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
              tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
              tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
              tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2")]
   if(dipole): geometry.append(tfile.Get("TPolyLine3D;1"))
   return geometry
#############################################

tf = TFile( storage+'/data/root/raw_'+proc+'_bkg.root', 'recreate' )
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
   nbkgtrks = 0
   while(nbkgtrks<NbkgTracks):
      ### draw vertex position
      R = cfgmap["Rbeampipe"]-cfgmap["Wbeampipe"] ## for now assume vetex can be only on the beampipe
      SigmaZ = 25 ## cm, around the dipole exit
      phi = rnd.Uniform(0,2*ROOT.TMath.Pi())
      vx0 = R*ROOT.TMath.Cos(phi)
      vy0 = R*ROOT.TMath.Sin(phi)
      vz0 = cfgmap["zDipoleExit"]+rnd.Gaus(0,SigmaZ)
      vz0 = round(vz0,1) ## precision of stepsize in z for propagation in B-filed is 0.1

      
      ### draw the momentum
      # tau = -0.5
      tau = -0.5
      P0 = rnd.Exp(tau)
      rp = TVector3()
      # rp.SetXYZ( rnd.Gaus(0,0.01), rnd.Gaus(0,0.01), rnd.Exp(tau) ) ### for now, cannot go backards in z...
      rp.SetXYZ( rnd.Gaus(0,0.01), rnd.Gaus(0,0.01), rnd.Exp(tau) ) ### for now, cannot go backards in z...
      
      r4 = cfgmap["zLayer1"]*rp.Unit()
      # print("x4=%g, y4=%g" % (r4.X(),r4.Y()))
      # if(abs(r4.X())<=cfgmap["Rbeampipe"]):                       continue ## |x| must be >Rbeampipe
      # if(abs(r4.X())>cfgmap["RoffsetBfield"]+2*cfgmap["Lstave"]): continue ## |x| must be <=Rbeampipe+2*Lstave
      # if(abs(r4.Y())>1.5*cfgmap["Hstave"]/2.):                    continue ## |y| must be within Hstave/2+50%
      
      nbkgtrks+=1
      
      Px = rp.Unit().X()*P0
      Py = rp.Unit().Y()*P0
      Pz = rp.Unit().Z()*P0
      
      E0 = ROOT.TMath.Sqrt(P0*P0+me2)
      p4 = TLorentzVector()
      p4.SetXYZM(Px,Py,Pz,me)
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




