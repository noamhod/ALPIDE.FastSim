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


def LegendMaker():
  legend1 = TLegend(0.80,0.70,0.92,0.88); ### with the mean
  legend1.SetFillColor(kWhite);
  legend1.SetTextFont(42);
  legend1.SetTextSize(0.046);
  legend1.SetBorderSize(0);
  legend1.SetShadowColor(kWhite);
  legend1.SetFillStyle(0);
  return legend1;

def main():
    
    gROOT.LoadMacro("LuxeStyle.C")
    gROOT.LoadMacro("LuxeLabels.C")
    gROOT.SetBatch()
    SetLuxeStyle()
    
    parser = argparse.ArgumentParser(description='Code to make track reco plots')
    parser.add_argument('-p', action="store", dest="proc", type=str, default="elaser") ### default is phase 2 signal
    parser.add_argument('-s', action="store", dest="sideNeeded", type=str, default="Pside")
    
    args    = parser.parse_args()    
    process = args.proc             ### "elaser/glaser"
    side    = args.sideNeeded        ### "Pside/Eside"

    tfnameSig = "../data/root/rec/"+process+"/phase0/sig/ppw/3.0/rec_"+process+"__"+side+".root"
    tfileSig = TFile(tfnameSig,"READ")
    
    tfnameSigBkg = "../data/root/rec/"+process+"/phase0/sigbkg/ppw/3.0/rec_"+process+"__"+side+".root"
    tfileSigBkg = TFile(tfnameSigBkg,"READ")
    
    tfnameBkg = "../data/root/rec/"+process+"/bkg/rec_"+process+"__"+side+".root"
    tfileBkg = TFile(tfnameBkg,"READ")

    hnames = ["h_px_rec_"+side,
              "h_px_rec_zoom_"+side,
              "h_py_rec_"+side,
              "h_py_rec_zoom_"+side,
              "h_pz_rec_"+side,
              "h_E_rec_all_"+side+"_log0",
              "h_E_rec_all_"+side+"_log1",
              "h_E_rec_all_"+side+"_log2",
              "h_E_rec_all_"+side+"_log3",
              "h_chi2_"+side,
              "h_chi2_mat_"+side,
              "h_chi2_non_"+side,
              "h_chi2dof_"+side,
              "h_SnpSig_"+side,
              "h_TglSig_"+side,
              "h_xVtxSig_"+side,
              "h_yVtxSig_"+side]
    
    
    outname = "signal_and_background_tracking_variable_"+process
    ### plot
    cnv = TCanvas("c","",1200,500)
    cnv.SaveAs(outname+".pdf(")
    for name in hnames:
        sigHist    = tfileSig.Get(name)
        sigBkgHist = tfileSigBkg.Get(name)
        bkgHist    = tfileBkg.Get(name)
        
        sigHist.SetLineColor(2)
        bkgHist.SetLineColor(4)
        sigBkgHist.SetLineColor(kGreen+3)

        
        sigHist.SetLineWidth(2)
        sigBkgHist.SetLineWidth(2)
        bkgHist.SetLineWidth(2)
        
        ## set the range
        nbinsX = sigHist.GetNbinsX()
        sigHist.GetXaxis().SetRangeUser(sigHist.GetXaxis().GetBinLowEdge(1), sigHist.GetBinContent(nbinsX))
        
        legend1 = LegendMaker()
        legend1.AddEntry(sigHist,"sig", "l")
        legend1.AddEntry(bkgHist,"bkg", "l")
        legend1.AddEntry(sigBkgHist,"sig+bkg", "l")
        
        cnv = TCanvas("c_"+name,"",1200,500)
        cnv.cd()
        ROOT.gPad.SetTicks(1,1)
        ROOT.gPad.SetGrid()
        
        LUXELabel(0.2,0.85,"CDR")
        
        sigBkgHist.Draw("hist")
        bkgHist.Draw("hist same")
        sigHist.Draw("hist same")
        
        legend1.Draw("sames")
        
        cnv.SaveAs(outname+".pdf")
        
    cnv = TCanvas("c","",1200,500)
    cnv.SaveAs(outname+".pdf)")
	    
	    
if __name__=="__main__":
    main()
		
		
