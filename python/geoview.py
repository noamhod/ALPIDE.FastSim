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
ROOT.gGeoManager.DrawTracks("ogl")

# TObjArray* tracksarr = gGeoManager->GetListOfTracks();
# Int_t ntracks = tracksarr->GetEntries();
# cout << "ntracks=" << ntracks << endl;
# for(int i=0 ; i<ntracks ; ++i)
# {
# 	tracksarr->At(i)->Draw("/*");
# 	cout << "drawn i=" << i << endl;
# 	tracksarr->At(i)->Print();
# }

print("choose `Quit ROOT` from the graphics file menue")
ROOT.gApplication.Run()

