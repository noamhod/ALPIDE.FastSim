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
argus = parser.parse_args()
proc  = argus.p

# ROOT.TGeoManager.Import("../output/root/GeoLUXE_"+proc+".root")


ROOT.TEveManager.Create()
ROOT.gEve.RegisterGeometryAlias("Default", "../output/root/GeoLUXE_"+proc+".root")
ROOT.gGeoManager = ROOT.gEve.GetDefaultGeometry()


print("before beampipe")
node_beampipe = ROOT.gGeoManager.GetTopVolume().FindNode("beampipe_1")
beampipe = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_beampipe)
beampipe.UseNodeTrans()
ROOT.gEve.AddGlobalElement(beampipe)
print("after beampipe")
node_dipole_hor1 = ROOT.gGeoManager.GetTopVolume().FindNode("dipole_hor_10")
dipole_hor1 = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_dipole_hor1)
dipole_hor1.UseNodeTrans()
ROOT.gEve.AddGlobalElement(dipole_hor1)
print("after dipole hor 1")
node_dipole_hor2 = ROOT.gGeoManager.GetTopVolume().FindNode("dipole_hor_11")
dipole_hor2 = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_dipole_hor2)
dipole_hor2.UseNodeTrans()
ROOT.gEve.AddGlobalElement(dipole_hor2)
print("after dipole hor 2")
node_dipole_ver1 = ROOT.gGeoManager.GetTopVolume().FindNode("dipole_ver_12")
dipole_ver1 = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_dipole_ver1)
dipole_ver1.UseNodeTrans()
ROOT.gEve.AddGlobalElement(dipole_ver1)
print("after dipole ver 1")
node_dipole_ver2 = ROOT.gGeoManager.GetTopVolume().FindNode("dipole_ver_13")
dipole_ver2 = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_dipole_ver2)
dipole_ver2.UseNodeTrans()
ROOT.gEve.AddGlobalElement(dipole_ver2)
print("after dipole ver 2")

for i in range(100,17+1):
   node_stave = ROOT.gGeoManager.GetTopVolume().FindNode("stave_"+str(i))
   stave = ROOT.TEveGeoTopNode(ROOT.gGeoManager, node_stave)
   stave.UseNodeTrans()
   ROOT.gEve.AddGlobalElement(stave)


tracks = ROOT.gGeoManager.GetListOfTracks()
ntracks = tracks.GetEntries()
for i in range(ntracks):
   npoints = tracks.At(i).GetNpoints()
   eveline = ROOT.TEveLine()
   for k in range(npoints):
      r = [ROOT.Double(),ROOT.Double(),ROOT.Double(),ROOT.Double()]
      tracks.At(i).GetPoint(k,r[0],r[1],r[2],r[3])
      eveline.SetNextPoint(r[0],r[1],r[2])
   if(npoints<10): eveline.SetLineColor(ROOT.kRed)
   else:           eveline.SetLineColor(ROOT.kBlack)
   ROOT.gEve.AddElement(eveline)

# TObjArray* tracksarr = gGeoManager->GetListOfTracks();
# Int_t ntracks = tracksarr->GetEntries();
# cout << "ntracks=" << ntracks << endl;
# for(int i=0 ; i<ntracks ; ++i)
# {
# 	tracksarr->At(i)->Draw("/*");
# 	cout << "drawn i=" << i << endl;
# 	tracksarr->At(i)->Print();
# }

ROOT.gEve.Redraw3D(ROOT.kTRUE)

print("choose `Quit ROOT` from the graphics file menue")
ROOT.gApplication.Run()