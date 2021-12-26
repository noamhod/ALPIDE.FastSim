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
   legend1 = TLegend(0.70,0.70,0.88,0.88); ### with the mean
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
   
   fnames = {
      "ss":"../data/root/rec/"+process+"/phase0/allpix/ppw/3.0s/rec_"+process+"__"+side+".root",
      "sb":"../data/root/rec/"+process+"/phase0/allpix/ppw/3.0sb/rec_"+process+"__"+side+".root",
      "bb":"../data/root/rec/"+process+"/bkg/rec_"+process+"__"+side+".root"
   }
   files = {}
   for ftype,fname in fnames.items(): files.update({ftype:TFile(fname,"READ")})
   
   hnames = ["h_all_csize_L1I_"+side,
             "h_all_csize_L2I_"+side,
             "h_all_csize_L3I_"+side,
             "h_all_csize_L4I_"+side,
             "h_all_csize_L1O_"+side,
             "h_all_csize_L2O_"+side,
             "h_all_csize_L3O_"+side,
             "h_all_csize_L4O_"+side,
             
             "h_all_csizex_L1I_"+side,
             "h_all_csizex_L2I_"+side,
             "h_all_csizex_L3I_"+side,
             "h_all_csizex_L4I_"+side,
             "h_all_csizex_L1O_"+side,
             "h_all_csizex_L2O_"+side,
             "h_all_csizex_L3O_"+side,
             "h_all_csizex_L4O_"+side,
             
             "h_all_csizey_L1I_"+side,
             "h_all_csizey_L2I_"+side,
             "h_all_csizey_L3I_"+side,
             "h_all_csizey_L4I_"+side,
             "h_all_csizey_L1O_"+side,
             "h_all_csizey_L2O_"+side,
             "h_all_csizey_L3O_"+side,
             "h_all_csizey_L4O_"+side,
             
             "h_rec_csize_L1I_"+side,
             "h_rec_csize_L2I_"+side,
             "h_rec_csize_L3I_"+side,
             "h_rec_csize_L4I_"+side,
             "h_rec_csize_L1O_"+side,
             "h_rec_csize_L2O_"+side,
             "h_rec_csize_L3O_"+side,
             "h_rec_csize_L4O_"+side,
             
             "h_rec_csizex_L1I_"+side,
             "h_rec_csizex_L2I_"+side,
             "h_rec_csizex_L3I_"+side,
             "h_rec_csizex_L4I_"+side,
             "h_rec_csizex_L1O_"+side,
             "h_rec_csizex_L2O_"+side,
             "h_rec_csizex_L3O_"+side,
             "h_rec_csizex_L4O_"+side,
             
             "h_rec_csizey_L1I_"+side,
             "h_rec_csizey_L2I_"+side,
             "h_rec_csizey_L3I_"+side,
             "h_rec_csizey_L4I_"+side,
             "h_rec_csizey_L1O_"+side,
             "h_rec_csizey_L2O_"+side,
             "h_rec_csizey_L3O_"+side,
             "h_rec_csizey_L4O_"+side,
             
             "h_rec_Nhits_"+side,
             "h_rec_px_"+side,
             "h_rec_px_zoom_"+side,
             "h_rec_py_"+side,
             "h_rec_py_zoom_"+side,
             "h_rec_E_"+side,
             "h_rec_E_"+side+"_log0",
             # "h_rec_E_"+side+"_log1",
             # "h_rec_E_"+side+"_log2",
             # "h_rec_E_"+side+"_log3",
             "h_rec_chi2dof_"+side,
             "h_rec_SnpSig_"+side,
             "h_rec_TglSig_"+side,
             "h_rec_xVtxSig_"+side,
             "h_rec_yVtxSig_"+side,
             "h_tru_rec_dErel_"+side,
             "h_tru_rec_E_ratio_"+side+"_log0",
             
             "h_mat_Nhits_"+side,
             "h_mat_px_"+side,
             "h_mat_px_zoom_"+side,
             "h_mat_py_"+side,
             "h_mat_py_zoom_"+side,
             "h_mat_E_"+side,
             "h_mat_E_"+side+"_log0",
             # "h_mat_E_"+side+"_log1",
             # "h_mat_E_"+side+"_log2",
             # "h_mat_E_"+side+"_log3",
             "h_mat_chi2dof_"+side,
             "h_mat_SnpSig_"+side,
             "h_mat_TglSig_"+side,
             "h_mat_xVtxSig_"+side,
             "h_mat_yVtxSig_"+side,
             "h_tru_mat_dErel_"+side,
             "h_tru_mat_E_ratio_"+side+"_log0",
             
             "h_non_Nhits_"+side,
             "h_non_px_"+side,
             "h_non_px_zoom_"+side,
             "h_non_py_"+side,
             "h_non_py_zoom_"+side,
             "h_non_E_"+side,
             "h_non_E_"+side+"_log0",
             # "h_non_E_"+side+"_log1",
             # "h_non_E_"+side+"_log2",
             # "h_non_E_"+side+"_log3",
             "h_non_chi2dof_"+side,
             "h_non_SnpSig_"+side,
             "h_non_TglSig_"+side,
             "h_non_xVtxSig_"+side,
             "h_non_yVtxSig_"+side,
             "h_tru_non_dErel_"+side,
             "h_tru_non_E_ratio_"+side+"_log0",
             
             "h_sel_Nhits_"+side,
             "h_sel_px_"+side,
             "h_sel_px_zoom_"+side,
             "h_sel_py_"+side,
             "h_sel_py_zoom_"+side,
             "h_sel_E_"+side,
             "h_sel_E_"+side+"_log0",
             # "h_sel_E_"+side+"_log1",
             # "h_sel_E_"+side+"_log2",
             # "h_sel_E_"+side+"_log3",
             "h_sel_chi2dof_"+side,
             "h_sel_SnpSig_"+side,
             "h_sel_TglSig_"+side,
             "h_sel_xVtxSig_"+side,
             "h_sel_yVtxSig_"+side,
             "h_tru_sel_dErel_"+side,
             "h_tru_sel_E_ratio_"+side+"_log0",
             
             "h_cutflow_"+side,
          ]
   
   
   outname = "signal_and_background_tracking_variable_"+process
   ### plot
   cnv = TCanvas("c","",700,500)
   cnv.SaveAs(outname+".pdf(")
   
   for name in hnames:
      htypes = {
         "ss":{"col":ROOT.kRed,     "leg":"Sig"},
         "sb":{"col":ROOT.kGreen+3, "leg":"Sig+Bkg"},
         "bb":{"col":ROOT.kBlue,    "leg":"Bkg"}
      }
      
      histos = {}
      
      cnv = TCanvas("c_"+name,"",700,500)
      cnv.cd()
      ROOT.gPad.SetTicks(1,1)
      ROOT.gPad.SetGrid()
      if("h_tru_" not in name): ROOT.gPad.SetLogy()
      
      LUXELabel(0.2,0.85,"TDR")
      legend = LegendMaker()
      
      hmax = -1e20
      for htyp,attr in htypes.items():
         hname = name+"_"+htyp
         
         ## draw truth just on energy plots
         if("_E_" in name and "_tru_" not in name and htyp=="ss"):
            name1 = name
            name1 = name1.replace("rec","tru").replace("sel","tru").replace("mat","tru").replace("non","tru")
            hname1 = name1+"_"+htyp
            hist = files[htyp].Get(name1).Clone(hname1)
            if(not hist): print(name1,"is null")
            histos.update( {hname1:hist} )
            histos[hname1].SetDirectory(0)
            hmax = histos[hname1].GetMaximum() if(histos[hname1].GetMaximum()>hmax) else hmax
            histos[hname1].SetLineColor(ROOT.kBlack)
            histos[hname1].SetLineWidth(2)
            histos[hname1].SetLineStyle(2)
            legend.AddEntry(histos[hname1],"Tru","l")
         
         hist = files[htyp].Get(name).Clone(hname)
         if(not hist): print(name,"is null")
         histos.update( {hname:hist} )
         histos[hname].SetDirectory(0)
         hmax = histos[hname].GetMaximum() if(histos[hname].GetMaximum()>hmax) else hmax
         histos[hname].SetLineColor(attr["col"])
         histos[hname].SetFillColorAlpha(attr["col"],0.35)
         histos[hname].SetLineWidth(1)
         legend.AddEntry(histos[hname],attr["leg"],"f")
         
      ## set max
      for hname,h in histos.items():
         h.SetMaximum(hmax*1.1)
         if("Nhits" in name or "cutflow" in name):
            h.SetMinimum(0.25)
            h.SetMaximum(hmax*2)
      
      ## draw truth just on energy plots
      if("_E_" in name and "_tru_" not in name):
         name1 = name
         name1 = name1.replace("rec","tru").replace("sel","tru").replace("mat","tru").replace("non","tru")
         htmp = histos[name1+"_ss"].Clone(name1+"_ss_tmp")
         htmp.GetXaxis().SetTitle( histos[name+"_sb"].GetXaxis().GetTitle() )
         htmp.GetYaxis().SetTitle( histos[name+"_sb"].GetYaxis().GetTitle() )
         htmp.Draw("hist")
         histos[name+"_sb"].Draw("hist same")
      else: histos[name+"_sb"].Draw("hist")
      if("mat" not in name and "non" not in name): histos[name+"_bb"].Draw("hist same")
      histos[name+"_ss"].Draw("hist same")
      legend.Draw("sames")
      
      cnv.RedrawAxis()
      cnv.Update()
      cnv.SaveAs(outname+".pdf")
      
   cnv = TCanvas("c","",700,500)
   cnv.SaveAs(outname+".pdf)")

if __name__=="__main__":
    main()
		
		
