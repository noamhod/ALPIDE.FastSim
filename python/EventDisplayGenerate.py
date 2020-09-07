#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TMath, TCanvas, TView, TView3D, TGraph2D, TStyle, TF2, TH1, TPolyLine3D, TRandom
import config as cfg
import EventDisplayBase as geo
import argparse

parser = argparse.ArgumentParser(description='GenerateGeoWithTracks.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-i', metavar='ievent', required=True,  help='event number [integer]')
argus = parser.parse_args()
proc  = argus.p
sides = "e+e-"
ievt  = int(argus.i)

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
storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")

cfgmap = cfg.set(proc,sides,True)

def gettracks():
   tfilename = storage+"/data/root/rec_from_clusters_"+proc+".root"
   tfile = TFile(tfilename,"READ")
   tree = tfile.Get("reco")
   stracks = []
   btracks = []
   n=0
   if(ievt<0 or ievt>tree.GetEntries()):
      print("ERROR: illegal event number:", ievt)
      quit()
   for event in tree:
      if(n==ievt):
         for i in range(event.true_trcklin.size()): stracks.append( event.true_trcklin[i].Clone() )
         for i in range(event.bkgr_trcklin.size()): btracks.append( event.bkgr_trcklin[i].Clone() )
         break
      n+=1
   for trk in stracks: trk.SetLineColor(ROOT.kYellow)
   for trk in btracks: trk.SetLineColor(ROOT.kRed)
   return (stracks,btracks)
   


stracks,btracks = gettracks()

geoluxe  = geo.GeoLUXE(proc,stracks,btracks)
world    = geoluxe.createWorld()
geoluxe.configureGeoManager(world)

tfileout = TFile(storage+"/output/root/bkgtrk_"+proc+".root","RECREATE")
tfileout.cd()

cnv1 = TCanvas("","",2000,2000)
view = TView3D.CreateView(1)
view.SetPerspective()
# view.SetParallel()
# view.ShowAxis()
view.SetRange(-80,-50,0, +80,+50,350)
geoluxe.draw(world)
cnv1.Write()
cnv1.SaveAs(storage+"/output/pdf/bkgtrk_"+proc+".pdf(")

cnv2 = TCanvas("","",2000,2000)
view = TView3D.CreateView(1)
view.SetPerspective()
# view.SetParallel()
# view.ShowAxis()
view.SetRange(cfgmap["xPsideL"],-10,190, cfgmap["xEsideR"],+10,340)
geoluxe.draw(world)
cnv2.Write()
cnv2.SaveAs(storage+"/output/pdf/bkgtrk_"+proc+".pdf)")

tfileout.mkdir("SigTracks/")
tfileout.cd("SigTracks/")
Ns = ROOT.TArrayI(1)
Ns[0] = len(stracks)
tfileout.WriteObject(Ns,"Ns")
for trk in stracks: trk.Write()
tfileout.mkdir("BkgTracks/")
tfileout.cd("BkgTracks/")
Nb = ROOT.TArrayI(1)
Nb[0] = len(btracks)
tfileout.WriteObject(Nb,"Nb")
for trk in btracks: trk.Write()

tfileout.cd()
tfileout.Write()
tfileout.Close()