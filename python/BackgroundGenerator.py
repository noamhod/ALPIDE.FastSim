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
parser.add_argument('-nnoisecls', metavar='noise clusters',  required=True,  help='N noise clusters')
parser.add_argument('-nbckgrtrk', metavar='background tracks',  required=True,  help='N background tracks')
parser.add_argument('-s', metavar='sides',   required=False, help='detector side [e+, e-, e+e-]')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.p
Nevt  = int(argus.nevents)
### background stuff
NnoiseClusters = int(argus.nnoisecls) ## 50  ## uniformly distributed clusters in x:y for each layer
NbkgTracks     = int(argus.nbckgrtrk) ## 200 ## uniformly distributed tracks coming from the inner part of the beampipe at the dipole-exit plane

sides = "e+" if(proc=="trident") else "e+e-" ## jsut the default
if(argus.s is not None): sides = argus.s
if(proc=="trident" and "e-" in sides):
   print("ERROR: do not run tracking in the electron side for trident")
   quit()
print("Running with proc=%s and sides=%s" % (proc,sides))

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError

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
B  = 1.4 if(proc=="trident") else 2.0 # Tesla
LB = 1   # meters

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
Lstave    = 27 if(proc=="bppp") else 50 # cm
Rbeampipe = 4 # cm
RoffsetBfield22BPPP = 7.0  # cm for BPPP in B=2.2T
RoffsetBfield20BPPP = 5.7  # cm for BPPP in B=2.0T
RoffsetBfield14BPPP = 4.0  # cm for BPPP in B=1.4T
RoffsetBfield = RoffsetBfield20BPPP if(proc=="bppp") else 14 # cm
xPsideL = -RoffsetBfield-Lstave
xPsideR = -RoffsetBfield       
xEsideL = +RoffsetBfield       
xEsideR = +RoffsetBfield+Lstave
yUp = +Hstave/2.
yDn = -Hstave/2.

### for the histogram
detXmin = xPsideL
detXmax = xEsideR
if(proc=="trident"): detXmax = xPsideR
if(proc=="bppp" and sides=="e+"): detXmax = xPsideR
if(proc=="bppp" and sides=="e-"): detXmin = xEsideL

layers = [1,2,3,4]

### histos
histos = {
   "h_dxrel":TH1D("h_dxrel",";(x_{cls}-x_{tru})/x_{tru};Tracks",200,-0.01,+0.01),
   "h_dyrel":TH1D("h_dyrel",";(y_{cls}-y_{tru})/y_{tru};Tracks",200,-0.25,+0.25),
}


#############################################
def getgeometry(dipole=False):
   tfile = TFile("../data/root/"+proc+"_geometry.root","READ")
   geometry = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
              tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
              tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
              tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2")]
   if(dipole): geometry.append(tfile.Get("TPolyLine3D;1"))
   return geometry

def rUnit2(r1,r2):
   r = (r2-r1).Unit()
   return r ## TVector2

def rUnit3(r1,r2):
   r = (r2-r1).Unit()
   return r ## TVector3

def p3(r1,r2,E):
   r = rUnit3(r1,r2)
   p = TVector3()
   p.SetXYZ(E*r.X(),E*r.Y(),E*r.Z())
   return p ## TVector3

def p4(r1,r2,E):
   p = p3(r1,r2,E)
   tlv = TLorentzVector()
   tlv.SetXYZM(p.Px(),p.Py(),p.Pz(),me)
   return tlv ## TLorentzVector

def xofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   if(dx==0):
      print("ERROR in xofz: dx=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
   a = dx/dz
   b = r1[0]-a*r1[2]
   x = a*z+b
   return x

def yofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dy = r2[1]-r1[1]
   if(dz==0):
      print("ERROR in yofz: dz=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
   a = dy/dz
   b = r1[1]-a*r1[2]
   y = a*z+b
   return y
   
def zofx(r1,r2,x):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   if(dx==0):
      print("ERROR in zofx: dx=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
   a = dz/dx
   b = r1[2]-a*r1[0]
   z = a*x+b
   return z

def countpoints(points3d):
   npoints = 0
   for i in range(points3d.GetN()):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      points3d.GetPoint(i,r[0],r[1],r[2])
      if(abs(r[0])<xDipoleExitMin): continue
      if(r[2]<300 or r[2]>330):     continue
      npoints += 1
   return npoints

def getpoint(pm, i):
   r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   pm.GetPoint(i,r[0],r[1],r[2])
   return r[0],r[1],r[2]

def getparticlename(cluster,i):
   r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   cluster.GetPoint(i,r[0],r[1],r[2])
   if  (r[0]>+1): return "electrons"
   elif(r[0]<-1): return "positrons"
   else:
      print("ERROR: x=",r[0])
      print("  cannot determine particle name")
      quit()
   return "unknown"

def isel(x):
   if(x<0): return False
   return True

### define the data structures to hold clusters and windows
def initpoints():
   points = { 
              "IsSig" : {1:[],2:[],3:[],4:[]},
              "TrkId" : {1:[],2:[],3:[],4:[]},
              "Cls"   : {1:TPolyMarker3D(),2:TPolyMarker3D(),3:TPolyMarker3D(),4:TPolyMarker3D()}
   }   
   return points
   
def getLogicSidesArr():
   sidesarr = ["Eside","Pside"]
   if(proc=="trident"): sidesarr = ["Pside"]
   if(proc=="bppp" and sides=="e+"): sidesarr = ["Pside"]
   if(proc=="bppp" and sides=="e-"): sidesarr = ["Eside"]
   return sidesarr

##########################################################################################
##########################################################################################
##########################################################################################

ROOT.gInterpreter.GenerateDictionary("vector<TPolyMarker3D>", "vector")
clusters_xyz  = ROOT.std.vector( 'TPolyMarker3D*' )()
clusters_type = ROOT.std.vector( int )() ## -1 for random noise or 0 for background track
clusters_id   = ROOT.std.vector( int )() ## the id of the generated background track (-1 for noise culsters)
tF = TFile("../data/root/background_clusters_"+proc+".root","RECREATE")
tF.cd()
tT = TTree("clusters","clusters")

tT.Branch('clusters_xyz',  clusters_xyz)
tT.Branch('clusters_type', clusters_type)
tT.Branch('clusters_id',   clusters_id)

sidesarr = getLogicSidesArr()

n=0 ### init n
for n in range(Nevt):
   ## clear the output vectors
   clusters_xyz.clear()
   clusters_type.clear()
   clusters_id.clear()
   
   ### embed some noise ***clusters***
   ids = [1000,2000,3000,4000] ## assuming there's no chance to have more than 1k noise clusters per layer
   layer2z = {1:300,2:310,3:320,4:330}
   rnd = TRandom()
   rnd.SetSeed()
   for side in sidesarr:
      for kN in range(NnoiseClusters):
         for layer in layers:
            x = rnd.Uniform(xPsideL,xPsideR) if(side=="Pside") else rnd.Uniform(xEsideL,xEsideR)
            y = rnd.Uniform(-0.75,+0.75)
            z = layer2z[layer]
            cluster = ROOT.TPolyMarker3D()
            cluster.SetNextPoint(x,y,z)
            clusters_xyz.push_back(cluster)
            clusters_type.push_back(-1)
            clusters_id.push_back( ids[layer-1] )
            ids[layer-1] += 1
            
   
   ### embed background ***tracks***
   resolution = 0.001 ## cm (10 um)
   rnd = TRandom()
   rnd.SetSeed()
   ids = [10000,20000,30000,40000] ## assuming there's no chance to have more than 100k tracks clusters per layer
   ibkgtracks = 0
   while ibkgtracks<NbkgTracks:
      # production vertex is uniform on the inner perimeter of the beampipe
      R = cfgmap["Rbeampipe"]-cfgmap["Wbeampipe"]
      phi = rnd.Uniform(0,2*ROOT.TMath.Pi())
      x0 = R*ROOT.TMath.Cos(phi)
      y0 = R*ROOT.TMath.Sin(phi)
      z0 = cfgmap["zDipoleExit"]
      # chose a point in the exit of the tracker 
      z4 = cfgmap["zLayer4"]
      x4 = rnd.Uniform(1.1*cfgmap["xPsideL"],   1.1*cfgmap["xEsideR"])
      y4 = rnd.Uniform(-1.1*cfgmap["Hstave"]/2, 1.1*cfgmap["Hstave"]/2)
      # require that the track is not pointing in the reange between the 2 trackr sides
      if(x4>cfgmap["xPsideR"] and x4<cfgmap["xEsideL"]): continue
      # place 3 points on layers 1,2,3 along the line between r4 and r0, with some error (resolution)
      r0 = [x0,y0,z0]
      r4 = [x4,y4,z4]
      z1 = cfgmap["zLayer1"]
      x1 = xofz(r4,r0,z1)+rnd.Gaus(0,resolution)
      y1 = yofz(r4,r0,z1)+rnd.Gaus(0,resolution)
      z2 = cfgmap["zLayer2"]
      x2 = xofz(r4,r0,z2)+rnd.Gaus(0,resolution)
      y2 = yofz(r4,r0,z2)+rnd.Gaus(0,resolution)
      z3 = cfgmap["zLayer3"]
      x3 = xofz(r4,r0,z3)+rnd.Gaus(0,resolution)
      y3 = yofz(r4,r0,z3)+rnd.Gaus(0,resolution)
      r1 = [x1,y1,z1]
      r2 = [x2,y2,z2]
      r3 = [x3,y3,z3]
      cluster = ROOT.TPolyMarker3D(); cluster.SetNextPoint(r1[0],r1[1],r1[2]); clusters_xyz.push_back(cluster)
      clusters_type.push_back(0); clusters_id.push_back(ids[0]); ids[0] += 1 #ibkgtracks)
      cluster = ROOT.TPolyMarker3D(); cluster.SetNextPoint(r2[0],r2[1],r2[2]); clusters_xyz.push_back(cluster)
      clusters_type.push_back(0); clusters_id.push_back(ids[1]); ids[1] += 1 #ibkgtracks)
      cluster = ROOT.TPolyMarker3D(); cluster.SetNextPoint(r3[0],r3[1],r3[2]); clusters_xyz.push_back(cluster)
      clusters_type.push_back(0); clusters_id.push_back(ids[2]); ids[2] += 1 #ibkgtracks)
      cluster = ROOT.TPolyMarker3D(); cluster.SetNextPoint(r4[0],r4[1],r4[2]); clusters_xyz.push_back(cluster)
      clusters_type.push_back(0); clusters_id.push_back(ids[3]); ids[3] += 1 #ibkgtracks)
      ibkgtracks+=1
      
   tT.Fill()
   if(n%100==0 and n>0): print("  processed %d events" % n)
print("Total events processed: ",n)


tF.cd()
tT.Write()
tF.Write()
tF.Close()





