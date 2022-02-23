#!/usr/bin/python
import os, sys
import math
import subprocess
import array
from array import array
import numpy as np
import ctypes
import ROOT
from ROOT import *
from copy import copy, deepcopy
import argparse
import pprint

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.11)
ROOT.gStyle.SetPadLeftMargin(0.12)
ROOT.gStyle.SetPadRightMargin(0.12)


def LegendMaker():
   legend1 = TLegend(0.5,0.75,0.80,0.92); ### with the mean
   legend1.SetFillColor(kWhite);
   legend1.SetTextFont(42);
   legend1.SetTextSize(0.046);
   legend1.SetBorderSize(0);
   legend1.SetShadowColor(kWhite);
   legend1.SetFillStyle(0);
   return legend1;


def label(text,x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1)
   s.SetTextFont(132)
   s.SetTextAlign(13)
   s.SetTextColor(col)
   s.SetTextSize(0.05)
   if(boldit): text = "#font[72]{"+text+"}"
   s.DrawLatex(x,y,text)


def graph(name,arrx,arry,col,mrk,xtitle,ytitle,witherr=False):
   htmp = TH1D("htmp","",1,0,1)
   htmp.SetBinErrorOption(ROOT.TH1.kPoisson)
   
   g = TGraphAsymmErrors()
   point = 0
   for point in range(len(arrx)):
      x = arrx[point]
      y = arry[point]
      htmp.SetBinContent(1,y)
      # yerr = math.sqrt(y) if(y>10) else htmp.GetBinError(1);
      yerr = math.sqrt(y)
      xerrdn = 0
      xerrup = 0
      # yerrdn = yerr if(y>10) else y-htmp.GetBinErrorLow(1)
      # yerrup = yerr if(y>10) else htmp.GetBinErrorUp(1)-y
      # print(y, math.sqrt(y), yerrdn, yerrup)
      yerrdn = yerr if(witherr) else 1e-6
      yerrup = yerr if(witherr) else 1e-6
      g.SetPoint(point,x,y)
      g.SetPointError(point, xerrdn,xerrup, yerrdn,yerrup)
   g.SetName(name)
   g.GetXaxis().SetTitle( xtitle )
   g.GetYaxis().SetTitle( ytitle )
   g.SetLineColor( col )
   g.SetLineWidth( 2 )
   g.SetMarkerColor( col )
   g.SetMarkerStyle( mrk )
   g.SetMarkerSize( 1 )
   return g
   

######################################################
######################################################
######################################################

def main():   
   gROOT.LoadMacro("LuxeStyle.C")
   gROOT.LoadMacro("LuxeLabels.C")
   gROOT.SetBatch()
   SetLuxeStyle()
   
   linarity = { "Truth":[], "Seeds":[], "Recos":[], "Selected":[], "Matched":[] }
   linarity["Truth"]    = [ 1.0,  5.0,    10.0, 50.0,   126.75, 250.0,  499.8, 999.4,  4989.5,  9951.2,  40453.5  ]
   linarity["Seeds"]    = [ 6675, 6676.9, 6687, 6709.2, 7140.5, 6920.3, 6863,  7670.3, 11506.4, 16405.7, 38415.25 ]
   linarity["Recos"]    = [ 1.3,  4.2,    8.0,  40.7,   162.5,  293.3,  456.5, 975.3,  4820.5,  9478.9,  33033.5  ]
   linarity["Selected"] = [ 0.9,  4.0,    7.7,  40.5,   111.5,  249.2,  439.8, 935.4,  4484.4,  9317.6,  30951.3  ]
   linarity["Matched"]  = [ 0.9,  4.0,    7.7,  40.5,   111.25, 211.1,  418.8, 828.2,  3778,    7499.9,  19824.25 ]
   
   ytitle = "Reconstructed track multiplicity per BX"
   xtitle = "Truth signal particle multiplicity per BX"
   
   grec = graph("grec", linarity["Truth"], linarity["Recos"],    ROOT.kRed,     24, xtitle,ytitle)
   gsel = graph("gsel", linarity["Truth"], linarity["Selected"], ROOT.kBlack,   20, xtitle,ytitle)
   gmat = graph("gmat", linarity["Truth"], linarity["Matched"],  ROOT.kGreen+2, 24, xtitle,ytitle)
   
   legend = LegendMaker()
   legend.AddEntry(grec,"Reconstructed","p")
   legend.AddEntry(gsel,"Selected","p")
   legend.AddEntry(gmat,"Matched","p")
   
   xmin = 1e-1
   xmax = 1e+5
   ymin = 1e-1
   ymax = 1e+5
   
   mg = TMultiGraph()
   mg.Add( grec,"EP" )
   mg.Add( gsel,"EP" )
   mg.Add( gmat,"EP" )
   mg.GetXaxis().SetTitle( xtitle )
   mg.GetYaxis().SetTitle( ytitle )
   mg.SetMinimum(ymin)
   mg.SetMaximum(ymax)
   mg.GetXaxis().SetLimits(xmin,xmax)
   
   f1 = TF1("f1","x",xmin,xmax)
   f1.SetLineColor(ROOT.kGray)
   f1.SetLineStyle(2)
   f1.SetLineWidth(2)
   
   
   cnv = TCanvas("cnv","",600,600)
   cnv.cd()
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetGrid()
   ROOT.gPad.SetLogx()
   ROOT.gPad.SetLogy()
   mg.Draw("ap")
   f1.Draw("same")
   mg.Draw("p")
   LUXELabel(0.2,0.85,"TDR")
   legend.Draw("sames")
   cnv.RedrawAxis()
   cnv.Update()
   cnv.SaveAs("linearity.pdf")
   


if __name__=="__main__":
    main()
		
		
