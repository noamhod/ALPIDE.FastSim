#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TMath, TCanvas, TView, TGraph2D, TStyle, TF2, TH1, TPolyLine3D
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


tfilename = "../data/root/rec_"+proc+".root"
tfile = TFile(tfilename,"READ")
tree = tfile.Get("res")
tracks = []
n=0
for event in tree:
   if(n==0):
      for i in range(event.poll.size()):
         tracks.append( event.poll[i].Clone() )
   else:     break
   n+=1


cfgmap = cfg.set(proc,sides,True)

geoluxe  = geo.GeoLUXE(proc,tracks)
world    = geoluxe.createWorld()
geoluxe.configureGeoManager(world)

cnv = TCanvas("","",2000,2000)
view = TView.CreateView(1)
view.ShowAxis()
view.SetRange(-80,-50,0, +80,+50,350)
geoluxe.draw(world)
cnv.SaveAs("bkgtrk_"+proc+".pdf(")

cnv = TCanvas("","",2000,2000)
view = TView.CreateView(1)
view.ShowAxis()
view.SetRange(cfgmap["xPsideL"],-10,190, cfgmap["xEsideR"],+10,340)
# dipole = world.FindNode("dipole")
# world.RemoveNode(dipole)
geoluxe.draw(world)
cnv.SaveAs("bkgtrk_"+proc+".pdf)")


tfileout = TFile("bkgtrk_"+proc+".root","RECREATE")
tfileout.cd()
world.Write()
tfileout.Write()
tfileout.Close()


