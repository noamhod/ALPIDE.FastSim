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

ROOT.TGeoManager.Import("../output/root/GeoLUXE_"+proc+".root")
# ROOT.gGeoManager.DefaultColors()
# ROOT.gGeoManager.SetMaxVisNodes(5000)
# ROOT.gGeoManager->SetVisLevel(1000000)
ROOT.gGeoManager.GetVolume("TOP").Draw("ogl")
# ROOT.gGeoManager.DrawTracks("ogl")

# TObjArray* tracksarr = gGeoManager->GetListOfTracks();
# Int_t ntracks = tracksarr->GetEntries();
# cout << "ntracks=" << ntracks << endl;
# for(int i=0 ; i<ntracks ; ++i)
# {
# 	tracksarr->At(i)->Draw("/*");
# 	cout << "drawn i=" << i << endl;
# 	tracksarr->At(i)->Print();
# }

tfile = ROOT.TFile("../output/root/bkgtrk_"+proc+".root","READ")
Ns = tfile.Get("Ns")[0]
Nb = tfile.Get("Nb")[0]
print("Ns=%g, Nb=%g" % (Ns,Nb))
for i in range(1,Ns+1):
   trk = tfile.Get("SigTracks/TPolyLine3D;"+str(i))
   eveline = TEveLine()
   for i in 
   eveline.SetNextPoint(120*sin(0.2*i), 120*cos(0.2*i), 80-i);
   
for i in range(1,Nb+1):
   trk = tfile.Get("BkgTracks/TPolyLine3D;"+str(i))
   trk.Draw("same")

print("choose `Quit ROOT` from the graphics file menue")
ROOT.gApplication.Run()