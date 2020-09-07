#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
import argparse

parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-t', metavar='tracks',  required=False, help='show tracks [y/n]')
argus = parser.parse_args()
proc  = argus.p
shwtrks = argus.t

ROOT.gErrorIgnoreLevel = ROOT.kWarning
# ROOT.gErrorIgnoreLevel = ROOT.kError
storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")




ROOT.TEveManager.Create()
ROOT.gEve.RegisterGeometryAlias("Default", "../output/root/GeoLUXE_"+proc+".root")
ROOT.gGeoManager = ROOT.gEve.GetDefaultGeometry()

# node_beampipe = ROOT.gGeoManager.GetTopVolume().FindNode("beampipe_1")
# beampipe = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_beampipe)
# beampipe.UseNodeTrans()
# ROOT.gEve.AddGlobalElement(beampipe)

nodes = ROOT.gGeoManager.GetTopVolume().GetNodes()
for node in nodes:
   element = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node)
   element.UseNodeTrans()
   ROOT.gEve.AddGlobalElement(element)

for i in range(100,17+1):
   node_stave = ROOT.gGeoManager.GetTopVolume().FindNode("stave_"+str(i))
   stave = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_stave)
   stave.UseNodeTrans()
   ROOT.gEve.AddGlobalElement(stave)


if(shwtrks!="n"):
   tracks = ROOT.gGeoManager.GetListOfTracks()
   ntracks = tracks.GetEntries()
   for i in range(ntracks):
      npoints = tracks.At(i).GetNpoints()
      eveline = ROOT.TEveLine()
      for k in range(npoints):
         r = [ROOT.Double(),ROOT.Double(),ROOT.Double(),ROOT.Double()]
         tracks.At(i).GetPoint(k,r[0],r[1],r[2],r[3])
         eveline.SetNextPoint(r[0],r[1],r[2])
      eveline.SetLineColor( tracks.At(i).GetLineColor() )
      # if(npoints<10): eveline.SetLineColor(ROOT.kRed)
      # else:           eveline.SetLineColor(ROOT.kYellow)
      ROOT.gEve.AddElement(eveline)
   

ROOT.gEve.Redraw3D(ROOT.kTRUE)

print("choose `Quit ROOT` from the graphics file menue")
ROOT.gApplication.Run()