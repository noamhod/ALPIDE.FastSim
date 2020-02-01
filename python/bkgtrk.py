#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TMath, TCanvas, TView, TView3D, TGraph2D, TStyle, TF2, TH1, TPolyLine3D, TRandom
import config as cfg
import geometry as geo
import argparse

parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
argus = parser.parse_args()
proc  = argus.p
sides = "e+e-"

'''
https://root.cern/doc/master/httpgeom_8C.html
https://root.cern.ch/doc/master/classTGeoManager.html
https://root.cern.ch/doc/master/classTGeoTube.html
https://root.cern.ch/doc/master/classTGeoVolume.html
https://root.cern.ch/doc/master/classTGeoTrack.html
https://root-forum.cern.ch/t/geometry-pyroot-tgeovolume-draw-problem-in-c-program-state-has-been-reset/26501/4
https://github.com/root-project/root/blob/master/documentation/users-guide/Geometry.md
'''

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)

cfgmap = cfg.set(proc,sides,True)

def sigtracks(ievt=0):
   tfilename = "../data/root/rec_"+proc+".root"
   tfile = TFile(tfilename,"READ")
   tree = tfile.Get("res")
   tracks = []
   n=0
   for event in tree:
      if(n==ievt):
         for i in range(event.poll.size()): tracks.append( event.poll[i].Clone() )
         break
      n+=1
   return tracks
   
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

def bkgtracks(N):
   tracks = []
   resolution = 0.001 ## cm (10 um)
   rnd = TRandom()
   rnd.SetSeed()
   for i in range(N):
      R = cfgmap["Rbeampipe"]-cfgmap["Wbeampipe"]
      phi = rnd.Uniform(0,2*ROOT.TMath.Pi())
      x0 = R*ROOT.TMath.Cos(phi)
      y0 = R*ROOT.TMath.Sin(phi)
      z0 = cfgmap["zDipoleExit"]
      
      z4 = cfgmap["zLayer4"]
      x4 = rnd.Uniform(1.15*cfgmap["xPsideL"],   1.15*cfgmap["xEsideR"])
      y4 = rnd.Uniform(-1.15*cfgmap["Hstave"]/2, 1.15*cfgmap["Hstave"]/2)
      
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
      
      z5 = 360
      x5 = xofz(r4,r0,z5)
      y5 = yofz(r4,r0,z5)
      
      line = TPolyLine3D()
      line.SetNextPoint(x0,y0,z0)
      line.SetNextPoint(x1,y1,z1)
      line.SetNextPoint(x2,y2,z2)
      line.SetNextPoint(x3,y3,z3)
      line.SetNextPoint(x4,y4,z4)
      line.SetNextPoint(x5,y5,z5)
      line.SetLineColor(ROOT.kRed)
      tracks.append(line)
      
   return tracks
   


stracks = sigtracks()
btracks = bkgtracks(100)

geoluxe  = geo.GeoLUXE(proc,stracks,btracks)
world    = geoluxe.createWorld()
geoluxe.configureGeoManager(world)

tfileout = TFile("../output/root/bkgtrk_"+proc+".root","RECREATE")
tfileout.cd()

cnv1 = TCanvas("","",2000,2000)
view = TView3D.CreateView(1)
view.SetPerspective()
# view.SetParallel()
# view.ShowAxis()
view.SetRange(-80,-50,0, +80,+50,350)
geoluxe.draw(world)
cnv1.Write()
cnv1.SaveAs("../output/pdf/bkgtrk_"+proc+".pdf(")

cnv2 = TCanvas("","",2000,2000)
view = TView3D.CreateView(1)
view.SetPerspective()
# view.SetParallel()
# view.ShowAxis()
view.SetRange(cfgmap["xPsideL"],-10,190, cfgmap["xEsideR"],+10,340)
geoluxe.draw(world)
cnv2.Write()
cnv2.SaveAs("../output/pdf/bkgtrk_"+proc+".pdf)")


tfileout.Write()
tfileout.Close()