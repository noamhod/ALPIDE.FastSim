#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, TCanvas, TLegend, TLatex, TGraph, TMultiGraph

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
storage = os.path.expandvars("$STORAGEDIR")

###############################################################

def minmax(h1,h2,f=1.1):
   hmin = h1.GetMinimum() if(h1.GetMinimum()<h2.GetMinimum()) else h2.GetMinimum()
   hmax = h1.GetMaximum() if(h1.GetMaximum()>h2.GetMaximum()) else h2.GetMaximum()
   h1.SetMaximum(hmax*f)
   h2.SetMaximum(hmax*f)
   return hmin,hmax*f

###############################################################

files = { "jeti40_trident":TFile("/storage/agrp/nhod/output/root/CDRresults_JETI40_e_laser_16.5GeV_.root","READ"),
          "jeti40_bppp"   :TFile("/storage/agrp/nhod/output/root/CDRresults_JETI40_g_laser_16.5GeV_.root","READ"),
          "phase2_trident":TFile("/storage/agrp/nhod/output/root/CDRresults_phaseII_e_laser_16.5GeV_.root","READ"),
          "phase2_bppp"   :TFile("/storage/agrp/nhod/output/root/CDRresults_phaseII_g_laser_16.5GeV_.root","READ"),
}
histos = {}

def geth(fname,hname):
   name = hname+"_"+fname
   if(fname not in files): print(fname," not found in files")
   histos.update( {name:files[fname].Get(hname).Clone(name)} )
   histos[name].SetDirectory(0)
   print("WARNING: setting unoccupoed bins to 1e6 instead of 9999")
   for b in range(histos[name].GetNbinsX()+1):
      if(histos[name].GetBinContent(b)<1e4): histos[name].SetBinError(b,0)
      if(histos[name].GetBinContent(b)>1e4): histos[name].SetBinContent(b,1e-6)

def seth(name,mrk,col):
   histos[name].SetLineColor(col)
   histos[name].SetMarkerColor(col)
   histos[name].SetMarkerStyle(mrk)
   histos[name].GetXaxis().SetTitle("Maximum #xi")
   histos[name].GetYaxis().SetTitle("Positrons/BX")

def LUXE(x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1)
   s.SetTextAlign(42)
   s.SetTextFont(52)
   s.SetTextColor(col)
   s.SetTextSize(0.044)
   s.DrawLatex(x,y,"LUXE #it{CDR}")

def getg(h):
   xlist = []
   ylist = []
   for b in range(h.GetNbinsX()+1):
      y = h.GetBinContent(b)
      x = h.GetBinCenter(b)
      if(y<1e4):
         xlist.append(x)
         ylist.append(y)
   x  = array.array('d', xlist)
   y  = array.array('d', ylist)
   g = TGraph(len(x),x,y)
   g.SetLineColor(h.GetLineColor())
   g.SetMarkerColor(h.GetMarkerColor())
   g.SetMarkerStyle(h.GetMarkerStyle())
   g.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
   g.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
   return g


geth("jeti40_trident","h_xi_acc")
geth("jeti40_bppp","h_xi_acc")
geth("phase2_trident","h_xi_acc")
geth("phase2_bppp","h_xi_acc")

seth("h_xi_acc_jeti40_trident",20,ROOT.kBlack)
seth("h_xi_acc_jeti40_bppp",24,ROOT.kBlack)
seth("h_xi_acc_phase2_trident",21,ROOT.kGray+1)
seth("h_xi_acc_phase2_bppp",25,ROOT.kGray+1)

g_jeti40_trident = getg(histos["h_xi_acc_jeti40_trident"])
g_jeti40_bppp    = getg(histos["h_xi_acc_jeti40_bppp"])
g_phase2_trident = getg(histos["h_xi_acc_phase2_trident"])
g_phase2_bppp    = getg(histos["h_xi_acc_phase2_bppp"])

###############################################################

allpdf  = storage+"/output/pdf/CDRresultsCombined.pdf"

leg_xi_trident = TLegend(0.15,0.63,0.53,0.80)
leg_xi_trident.SetFillStyle(4000) # will be transparent
leg_xi_trident.SetFillColor(0)
leg_xi_trident.SetTextFont(42)
leg_xi_trident.SetBorderSize(0)
leg_xi_trident.AddEntry(histos["h_xi_acc_jeti40_trident"],"e+laser JETI40","p")
leg_xi_trident.AddEntry(histos["h_xi_acc_phase2_trident"],"e+laser PhaseII","p")
leg_xi_trident.AddEntry(histos["h_xi_acc_jeti40_bppp"],"#gamma+laser JETI40","p")
leg_xi_trident.AddEntry(histos["h_xi_acc_phase2_bppp"],"#gamma+laser PhaseII","p")

cnv = TCanvas("cnv","",600,500)
cnv.cd()
cnv.SetTicks(1,1)
cnv.SetLogy()

#histos["h_xi_acc_phase2_trident"].SetMaximum(5e+4)
#histos["h_xi_acc_phase2_trident"].SetMinimum(5e-4)

#histos["h_xi_acc_phase2_trident"].Draw("e1p")
#histos["h_xi_acc_jeti40_trident"].Draw("e1p same")
#histos["h_xi_acc_phase2_bppp"].Draw("e1p same")
#histos["h_xi_acc_jeti40_bppp"].Draw("e1p same")

g_jeti40_trident.SetMaximum(5e+4)
g_jeti40_trident.SetMinimum(5e-4)

mgr = TMultiGraph()
mgr.Add(g_jeti40_trident,"p0")
mgr.Add(g_phase2_trident,"p0")
mgr.Add(g_jeti40_bppp,"p0")
mgr.Add(g_phase2_bppp,"p0")
mgr.Draw("ap0")
mgr.SetMaximum(5e+4)
mgr.SetMinimum(5e-4)
mgr.GetXaxis().SetTitle(g_jeti40_trident.GetXaxis().GetTitle())
mgr.GetYaxis().SetTitle(g_jeti40_trident.GetYaxis().GetTitle())
mgr.GetXaxis().SetRangeUser(0.,6.1)

#g_jeti40_trident.Draw("AP0")
#g_phase2_trident.Draw("AC")
#g_jeti40_bppp.Draw("AC")
#g_phase2_bppp.Draw("AC")

leg_xi_trident.Draw("same")

LUXE(0.18,0.85,ROOT.kBlack)


#s = TLatex()
#s.SetNDC(1);
#s.SetTextAlign(13);
#s.SetTextColor(ROOT.kBlack)
#s.SetTextSize(0.035)
#s.DrawLatex(0.17,0.27,"Stat error only")
#s.DrawLatex(0.17,0.22,"set to Poisson")

cnv.SaveAs(allpdf)
