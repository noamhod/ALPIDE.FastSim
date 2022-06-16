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


def graph(name,arrx,arry,col,mrk,xtitle,ytitle,errx,erry,witherr=False):
   htmp = TH1D("htmp","",1,0,1)
   htmp.SetBinErrorOption(ROOT.TH1.kPoisson)
   
   g = TGraphAsymmErrors()
   point = 0
   for point in range(len(arrx)):
      x = arrx[point]
      y = arry[point]
      htmp.SetBinContent(1,y)
      ### here forced symmetric error
      yerr = erry[point]
      xerrdn = errx[point]
      xerrup = errx[point]
      yerrdn = yerr
      yerrup = yerr
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
   
   linarity = { "Truth":[], "TruthError":[], "Npix":[], "NpixError":[] }
   
   linarity["Truth"] = [ 499.8, 999.4, 4989.5, 9951.2, 40453.5 ] #126.75,  250,  
   linarity["TruthError"] = [ 0,0,0,0,0 ] #0,0
   linarity["Npix"]  = [ 0.06,  0.93,  1.64,   1.72,   1.75 ] #2.11,  -1.68,  
   linarity["NpixError"] = [ 0.3463450308, 0.2159309193, 0.06685101419, 0.0422737233, 0.01786155703] #1.620009288, 0.605509768, 

   
   ytitle = "(N_{pix}^{s+b}-N_{pix}^{b})/N_{tru} per BX"
   xtitle = "Truth signal particle multiplicity per BX"
   
   gnpix = graph("gnpix", linarity["Truth"], linarity["Npix"], ROOT.kBlack,  20, xtitle,ytitle,linarity["TruthError"],linarity["NpixError"])
   
   xmin = 3e+2
   xmax = 1e+5
   ymin = -1
   ymax = +2.5
   
   mg = TMultiGraph()
   mg.Add( gnpix,"EP" )
   mg.GetXaxis().SetTitle( xtitle )
   mg.GetYaxis().SetTitle( ytitle )
   mg.SetMinimum(ymin)
   mg.SetMaximum(ymax)
   mg.GetXaxis().SetLimits(xmin,xmax)
   
   # f1 = TF1("f1","x",xmin,xmax)
   # f1.SetLineColor(ROOT.kGray)
   # f1.SetLineStyle(2)
   # f1.SetLineWidth(2)
   
   
   cnv = TCanvas("cnv","",600,600)
   cnv.cd()
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetGrid()
   ROOT.gPad.SetLogx()
   mg.Draw("ap")
   LUXELabel(0.5,0.85,"TDR")
   # legend.Draw("sames")
   cnv.RedrawAxis()
   cnv.Update()
   cnv.SaveAs("counter.pdf")
   


if __name__=="__main__":
    main()
		
		
