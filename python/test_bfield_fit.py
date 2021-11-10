#!/usr/bin/python
import os
import math
import subprocess
import array
from array import array
import numpy as np
import ctypes
import ROOT
from ROOT import TFile, TH1D, TH2D, TF1, TCanvas, TPad, TGraphAsymmErrors

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


f = TFile("test_bfield.root","READ")

hxE = f.Get("h2_E_vs_x")

hxE_q50 = hxE.QuantilesX(0.5)
hxE_q50.SetLineColor(ROOT.kBlack)
hxE_q50.SetMarkerColor(ROOT.kBlack)
hxE_q50.SetMarkerStyle(24)
hxE_q50.SetMarkerSize(0.8)

hxE_q99 = hxE.QuantilesX(0.99)
hxE_q99.SetLineColor(ROOT.kBlack)
hxE_q01 = hxE.QuantilesX(0.01)
hxE_q01.SetLineColor(ROOT.kBlack)


x,y = array('d'), array('d')
eyh,eyl = array('d'), array('d')
exh,exl = array('d'), array('d')
for b in range(1,hxE_q50.GetNbinsX()+1):
   xc  = hxE_q50.GetBinCenter(b)
   q50 = hxE_q50.GetBinContent(b) 
   q99 = hxE_q99.GetBinContent(b)
   q01 = hxE_q01.GetBinContent(b)
   if(q50<1e-3): continue
   x.append( xc )
   y.append( q50 )
   eyh.append( q99-q50 )
   eyl.append( q50-q01 )
   exh.append( 0 )
   exl.append( 0 )

gr = TGraphAsymmErrors(len(x),x,y,exl,exh,eyl,eyh)
gr.SetName("quantiles")
# gr.SetTitle(titles)
gr.SetLineColor( ROOT.kBlack )
gr.SetLineWidth( 1 )
gr.SetMarkerColor( ROOT.kBlack )
gr.SetMarkerStyle( 24 )
gr.SetMarkerSize( 0.9 )

print("\nPositron side:")
fBvsx_pos = TF1("fBvsx_pos", "[0]+[1]/([2]+x)", 1,20)
fBvsx_pos.SetLineColor(ROOT.kRed)
res = gr.Fit(fBvsx_pos,"MERS")
chi2dof = fBvsx_pos.GetChisquare()/fBvsx_pos.GetNDF() if(fBvsx_pos.GetNDF()>0) else -1
print("fBvsx_pos: chi2/Ndof=",chi2dof)

print("\nElectron side:")
fBvsx_ele = TF1("fBvsx_ele", "[0]+[1]/([2]+x)", -20,1)
fBvsx_ele.SetLineColor(ROOT.kRed)
res = gr.Fit(fBvsx_ele,"MERS")
chi2dof = fBvsx_ele.GetChisquare()/fBvsx_ele.GetNDF() if(fBvsx_ele.GetNDF()>0) else -1
print("fBvsx_ele: chi2/Ndof=",chi2dof)

cnv = TCanvas("c1","",1000,1000)
hxE.Draw("scat")
gr.Draw("p")
fBvsx_pos.Draw("same")
fBvsx_ele.Draw("same")
cnv.SaveAs("test_bfield_fit.pdf")

fout = TFile("test_bfield_fit.root","RECREATE")
fout.cd()
hxE.Write()
hxE_q50.Write()
hxE_q99.Write()
hxE_q01.Write()
gr.Write()
fBvsx_pos.Write()
fBvsx_ele.Write()
fout.Write()
fout.Close()
