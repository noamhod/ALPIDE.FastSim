#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ctypes
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TPad, TView, TLatex, TLegend, TGaxis

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

Ly  = 59.92/10 # cm
Lz  = 1440/10 # cm
Lx  = (980-2*(202+70))/10 # cm
zDipoleCenter = 2050/10 # cm
zDipoleExit = zDipoleCenter+Lz/2
stepsize = 0.1
NdipoleExit = int(zDipoleExit/stepsize)

hxE = TH2D("h2_E_vs_x",";x [cm];E [GeV];Tracks", 1000,-20,+20, 200, 0,  17)
hxz = TH2D("h2_z_vs_x",";x [cm];z [cm];Tracks",  1000,-20,+20, 2000,zDipoleExit-2,zDipoleExit+2)
hyz = TH2D("h2_z_vs_y",";y [cm];z [cm];Tracks",  1000,-0.5,+0.5, 2000,zDipoleExit-2,zDipoleExit+2)
hyx = TH2D("h2_y_vs_x",";x [cm];y [cm];Tracks",  1000,-20,+20, 2000,-0.5,+0.5)

#tfile = TFile("../data/root/dig/dig_bppp_0000001.root","READ")
tfile = TFile("../data/root/dig/dig_bppp_.root","READ")
ttree = tfile.Get("dig")
Ntotal = ttree.GetEntries()

Nprint = 2
Nevents = 0
for event in ttree:
   for i in range(event.trkpts_fullrange.size()):
      #print(i)
      #for n in range(event.trkpts_fullrange[i].GetN()):
      for n in range(NdipoleExit-2,NdipoleExit+2):

         xctype = ctypes.c_double() #ROOT.Double()
         yctype = ctypes.c_double() #ROOT.Double()
         zctype = ctypes.c_double() #ROOT.Double()

         event.trkpts_fullrange[i].GetPoint(n,xctype,yctype,zctype)
         x = xctype.value
         y = yctype.value
         z = zctype.value

         E = event.trkp4[i].E()
         if(abs(z-(zDipoleCenter+Lz/2.))<stepsize/2):
            hxz.Fill(x,z)
            hyz.Fill(y,z)
            hyx.Fill(x,y)
            hxE.Fill(x,E)
   if(Nevents%Nprint==0): print("Done %g out of %g" % (Nevents,Ntotal))
   Nevents += 1

cnv = TCanvas("c1","",1000,1000)
hxE.SetMarkerColor(ROOT.kViolet)
hxE.Draw("scat")
cnv.SaveAs("test_bfield.pdf(")

cnv = TCanvas("c1","",1000,1000)
hxz.SetMarkerColor(ROOT.kViolet)
hxz.Draw("scat")
cnv.SaveAs("test_bfield.pdf")

cnv = TCanvas("c1","",1000,1000)
hyz.SetMarkerColor(ROOT.kViolet)
hyz.Draw("scat")
cnv.SaveAs("test_bfield.pdf")

cnv = TCanvas("c1","",1000,1000)
hyx.SetMarkerColor(ROOT.kViolet)
hyx.Draw("scat")
cnv.SaveAs("test_bfield.pdf)")

fout = TFile("test_bfield.root","RECREATE")
fout.cd()
hxE.Write()
hxz.Write()
hyz.Write()
hyx.Write()
fout.Write()
fout.Close()
