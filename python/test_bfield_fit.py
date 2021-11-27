#!/usr/bin/python
import os
import math
import subprocess
import array
from array import array
import numpy as np
import ctypes
import ROOT
from ROOT import TFile, TH1D, TH2D, TF1, TCanvas, TPad, TGraphAsymmErrors, TLegend, TLatex

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

def label(text,x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1)
   s.SetTextAlign(13)
   s.SetTextColor(col)
   s.SetTextSize(0.044)
   if(boldit): text = "#font[72]{"+text+"}"
   s.DrawLatex(x,y,text)

Ly  = 59.92/10 # cm
Lz  = 1440/10 # cm
Lx  = (980-2*(202+70))/10 # cm
zDipoleCenter = 2050/10 # cm
zDipoleExit = zDipoleCenter+Lz/2
stepsize = 0.1
NdipoleExit = int(zDipoleExit/stepsize)


f = TFile("test_bfield.root","READ")

hxE = f.Get("h2_E_vs_x")
hxE.GetXaxis().SetTitle("x_{at dipole exit} [cm]")
hxE.GetYaxis().SetTitle("E_{truth} [GeV]")

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
gr.SetMarkerSize( 0.1 )

print("\nPositron side:")
fEvsX_pos = TF1("fEvsX_pos", "[0]+[1]/([2]+x)", 1,20)
fEvsX_pos.SetLineColor(ROOT.kRed)
fEvsX_pos.SetLineWidth(1)
res_pos = gr.Fit(fEvsX_pos,"MERS")
chi2dof = fEvsX_pos.GetChisquare()/fEvsX_pos.GetNDF() if(fEvsX_pos.GetNDF()>0) else -1
print("fEvsX_pos: chi2/Ndof=",chi2dof)

print("\nElectron side:")
fEvsX_ele = TF1("fEvsX_ele", "[0]+[1]/([2]+x)", -20,1)
fEvsX_ele.SetLineColor(ROOT.kBlue)
fEvsX_ele.SetLineWidth(1)
res_ele = gr.Fit(fEvsX_ele,"MERS")
chi2dof = fEvsX_ele.GetChisquare()/fEvsX_ele.GetNDF() if(fEvsX_ele.GetNDF()>0) else -1
print("fEvsX_ele: chi2/Ndof=",chi2dof)


leg1 = TLegend(0.55,0.8,0.87,0.87)
# leg1.SetFillStyle(4000) # will be transparent
leg1.SetFillColor(0)
leg1.SetTextFont(42)
leg1.SetBorderSize(0)
hxE.SetFillColor(hxE.GetMarkerColor())
hxE.SetLineColor(hxE.GetMarkerColor())
leg1.AddEntry(hxE,"Simulation","F")



leg2 = TLegend(0.55,0.6,0.87,0.87)
# leg2.SetFillStyle(4000) # will be transparent
leg2.SetFillColor(0)
leg2.SetTextFont(42)
leg2.SetBorderSize(0)
leg2.AddEntry(gr,"Median^{+Q99%}_{-Q1%}","f")
leg2.AddEntry(fEvsX_pos,"Fit (e^{+} side)","f")
leg2.AddEntry(fEvsX_ele,"Fit (e^{-} side)","f")

cnv = TCanvas("c1","",1200,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetGrid()
hxE.Draw("scat")
leg1.Draw("same")
label("B_{y}^{0} = -0.95 T",0.2,0.83)
label("(non-uniform)",0.2,0.75)
cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetGrid()
hxE.Draw("scat")
gr.Draw("ps")
fEvsX_pos.Draw("same")
fEvsX_ele.Draw("same")
leg2.Draw("same")
label("B_{y}^{0} = -0.95 T",0.2,0.83)
label("(non-uniform)",0.2,0.75)
cnv.SaveAs("test_bfield_fit.pdf")

fout = TFile("test_bfield_fit.root","RECREATE")
fout.cd()
hxE.Write()
hxE_q50.Write()
hxE_q99.Write()
hxE_q01.Write()
gr.Write()
fEvsX_pos.Write()
fEvsX_ele.Write()
fout.Write()
fout.Close()
