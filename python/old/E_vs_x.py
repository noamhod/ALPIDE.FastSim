#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)

#############################################


histos = {}


### get all histos
tfile_bppp    = TFile("../output/root/all_bppp.root","READ")
tfile_trident = TFile("../output/root/all_trident.root","READ")
histos = [ tfile_bppp.Get("h_xE_layer0_gen"),
           tfile_bppp.Get("h_xE_layer1_gen"),
           tfile_bppp.Get("h_xE_layer2_gen"),
           tfile_bppp.Get("h_xE_layer3_gen"),
           tfile_bppp.Get("h_xE_layer4_gen"), 
         ]
histos[0].Add(tfile_trident.Get("h_xE_layer0_gen"))
histos[1].Add(tfile_trident.Get("h_xE_layer1_gen"))
histos[2].Add(tfile_trident.Get("h_xE_layer2_gen"))
histos[3].Add(tfile_trident.Get("h_xE_layer3_gen"))
histos[4].Add(tfile_trident.Get("h_xE_layer4_gen"))


### get the profiles from the 2d histos
profiles = []
for h in histos:
   profiles.append(h.ProfileX())


### define the fit functions
fEle = [ TF1("Ele_EofX0","[0]/x+[1]/x^2+[2]/x^3",1,30), ### good
         TF1("Ele_EofX1","[0]/x+[1]/x^2+[2]/x^3",4,75), ### good
         TF1("Ele_EofX2","[0]/x+[1]/x^2+[2]/x^3",4,75), ### good
         TF1("Ele_EofX3","[0]/x+[1]/x^2+[2]/x^3",4,75), ### good
         TF1("Ele_EofX4","[0]/x+[1]/x^2+[2]/x^3",4,75), ### good
       ]
fPos = [ TF1("Pos_EofX0","[0]/x+[1]/x^2+[2]/x^3",-1,-30), ### good
         TF1("Pos_EofX1","[0]/x+[1]/x^2+[2]/x^3",-4,-75), ### good
         TF1("Pos_EofX2","[0]/x+[1]/x^2+[2]/x^3",-4,-75), ### good
         TF1("Pos_EofX3","[0]/x+[1]/x^2+[2]/x^3",-4,-75), ### good
         TF1("Pos_EofX4","[0]/x+[1]/x^2+[2]/x^3",-4,-75), ### good
       ]


### make the fits 
l=0
for p in profiles:
   p.SetMaximum(20.)
   p.SetMarkerStyle(24)
   p.SetMarkerSize(1)
   
   print("\nElectron side, layer:",l)
   r = p.Fit("Ele_EofX"+str(l),"WEMRS")
   chi2dof = fEle[l].GetChisquare()/fEle[l].GetNDF()
   print("chi2/dof=",chi2dof)
   
   print("\nPositron side, layer:",l)
   r = p.Fit("Pos_EofX"+str(l),"WEMRS")
   chi2dof = fPos[l].GetChisquare()/fPos[l].GetNDF()
   print("chi2/dof=",chi2dof)

   l+=1


### plot the  profiles with the fits
cnv = TCanvas("cnv_xE","",1000,1000)
cnv.Divide(1,5)
p1 = cnv.cd(1)
p2 = cnv.cd(2)
p3 = cnv.cd(3)
p4 = cnv.cd(4)
p5 = cnv.cd(5)
p1.cd()
profiles[4].Draw()
fEle[4].Draw("same")
fPos[4].Draw("same")
p2.cd()
profiles[3].Draw()
fEle[3].Draw("same")
fPos[3].Draw("same")
p3.cd()
profiles[2].Draw()
fEle[2].Draw("same")
fPos[2].Draw("same")
p4.cd()
profiles[1].Draw()
fEle[1].Draw("same")
fPos[1].Draw("same")
p5.cd()
profiles[0].Draw()
fEle[0].Draw("same")
fPos[0].Draw("same")
cnv.SaveAs("../output/pdf/profile_E_vs_x.pdf")


### write fits to root file
foutname = "../output/root/fits_E_vs_x.root"
tfile = TFile(foutname,"RECREATE")
tfile.cd()
for f in fEle: f.Write()
for f in fPos: f.Write()
tfile.Write()
tfile.Close()
print("Fits written to: "+foutname)