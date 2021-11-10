#!/usr/bin/python
import os
import math
import subprocess
import array
from array import array
import numpy as np
from operator import truediv
import ROOT
from ROOT import TFile, TH1D, TH2D, TCanvas, TLegend, TGraphErrors, TF1, TF2, TF3, TLatex, TMultiGraph, TArrow


ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kRainBow) # see https://root.cern.ch/doc/master/classTColor.html#C06
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.18)
ROOT.gStyle.SetPadRightMargin(0.2)


Bofz = "1/((1+exp(([0]-y)/[1]))*(1+exp((y-[2])/[3])))"
Bofx = "1/((1+exp(([4]-x)/[5]))*(1+exp((x-[6])/[7])))"

pz = [-616.0, 28.66,  622.0, 28.91]
px = [-165.0, 7.7, 165.0, 7.7]

Ly  = 59.92
Lz  = 1440
Lx  = 980-2*(202+70)
B0y = -0.95
LxEff = #Lx #2*170
LzEff = Lz

ymin = -Ly/2
ymax = +Ly/2
zmin = -Lz/2
zmax = +Lz/2
xmin = -Lx/2
xmax = +Lx/2



B = TF2("B",str(B0y)+"*"+Bofz+"*"+Bofx, xmin,xmax, zmin,zmax)

B.SetNpx(1000)
B.SetNpy(1000)

B.SetParameter(0,pz[0])
B.SetParameter(1,pz[1])
B.SetParameter(2,pz[2])
B.SetParameter(3,pz[3])

B.SetParameter(4,px[0])
B.SetParameter(5,px[1])
B.SetParameter(6,px[2])
B.SetParameter(7,px[3])

B.GetXaxis().SetTitle("x [mm]")
B.GetYaxis().SetTitle("z [mm]")
B.GetZaxis().SetTitle("B_{y} [T]")
B.SetTitle("")

c = TCanvas("c","",700,1500)
B.Draw("colz")
c.SetTicks(1,1)
c.SetGrid()
c.RedrawAxis()
c.SaveAs("Bfield.pdf")

Beff = B.Integral(-LxEff/2,+LxEff/2,-LzEff/2,+LzEff/2)/(LxEff*LzEff)
print("Beff=",Beff,"[T]")

