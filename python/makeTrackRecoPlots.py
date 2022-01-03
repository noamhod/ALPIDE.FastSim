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
   legend1 = TLegend(0.70,0.70,0.88,0.88); ### with the mean
   legend1.SetFillColor(kWhite);
   legend1.SetTextFont(42);
   legend1.SetTextSize(0.046);
   legend1.SetBorderSize(0);
   legend1.SetShadowColor(kWhite);
   legend1.SetFillStyle(0);
   return legend1;


def hChopperUp(h, binstochop):
   nbinsorig = h.GetNbinsX()
   nbins = nbinsorig-binstochop
   print("nbins=",nbins)
   xbins = []
   xaxis = h.GetXaxis()
   yaxis = h.GetYaxis()
   for b in range(1,nbins+1): xbins.append( xaxis.GetBinLowEdge(b) )
   xbins.append( xaxis.GetBinUpEdge(nbins) )
   arrxbins = array("d", xbins)
   print(arrxbins)
   name   = h.GetName()
   title  = h.GetTitle()
   xtitle = xaxis.GetTitle()
   ytitle = yaxis.GetTitle()
   hChopped = TH1D(name+"_chop",";"+xtitle+";"+ytitle,nbins,arrxbins)
   for b in range(1,nbins+1):
      hChopped.SetBinContent(b, h.GetBinContent(b))
      hChopped.SetBinError(b, h.GetBinError(b))
      hChopped.GetXaxis().SetBinLabel(b, h.GetXaxis().GetBinLabel(b))
   hChopped.SetLineColor(h.GetLineColor())
   hChopped.SetFillColor(h.GetFillColor())
   hChopped.SetLineStyle(h.GetLineStyle())
   hChopped.SetLineWidth(h.GetLineWidth())
   return hChopped

# def hChopperUp(h,binstokeep=-1):
#    hChopped = TH1D()
#    hChopped.SetName(h.GetName()+"_chop")
#    for b in range(1,h.GetNbinsX()):
#       label = h.GetXaxis().GetBinLabel(b)
#       print(b,label)
#       if(label==""): break
#       hChopped.Fill(label,h.GetBinContent(b))
#    hChopped.SetTitle(";"+h.GetXaxis().GetTitle()+";"+h.GetYaxis().GetTitle())
#    hChopped.SetLineColor(h.GetLineColor())
#    hChopped.SetFillColor(h.GetFillColor())
#    hChopped.SetLineStyle(h.GetLineStyle())
#    hChopped.SetLineWidth(h.GetLineWidth())
#    return hChopped


def main():   
   gROOT.LoadMacro("LuxeStyle.C")
   gROOT.LoadMacro("LuxeLabels.C")
   gROOT.SetBatch()
   SetLuxeStyle()
   
   parser = argparse.ArgumentParser(description='makeTrackRecoPlots.py...')
   parser.add_argument('-proc', metavar='process', required=True,  help='process [elaser,glaser]')
   parser.add_argument('-smpl', metavar='sample',  required=True,  help='sample [phase0/allpix/ppw/3.0]')
   parser.add_argument('-side', metavar='side',    required=False, help='side [,Pside,Eside]')
   argus  = parser.parse_args()
   process = argus.proc
   sample  = argus.smpl
   sidepe  = argus.side
   
   storage = os.getenv("$STORAGEDIR","../")
   print("storage=",storage)
   pdfsdir = storage+"/output/pdf/rec/"+process+"/"+sample+"/"
   os.makedirs(pdfsdir, exist_ok=True)
   
   outname = "signal_and_background_tracking_variable_"+process
   
   sides = ["Eside","Pside"] if(sidepe==None) else [sidepe]
   
   for side in sides:
      
      fnames = {
         "ss":"../data/root/rec/"+process+"/"+sample+"s/rec_"+process+"__"+side+".root",
         "sb":"../data/root/rec/"+process+"/"+sample+"sb/rec_"+process+"__"+side+".root",
         "bb":"../data/root/rec/"+process+"/bkg/rec_"+process+"__"+side+".root"
      }
      
      files = {}
      nBXs  = {}
      for ftype,fname in fnames.items():
         files.update({ftype:TFile(fname,"READ")})
         nBXs.update({ftype:files[ftype].Get("h_nBX_"+side).GetBinContent(1)})
      
      
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
                "h_rec_pz_"+side,
                "h_rec_pz_"+side+"_log0",
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
                "h_resol_rec_dpzrel_"+side,
                "h_resol_rec_dErel_"+side,
                "h_ratio_rec_pz_"+side+"_log0",
                "h_ratio_rec_pz_"+side+"_log1",
                "h_ratio_rec_E_"+side+"_log0",
                "h_ratio_rec_E_"+side+"_log1",
                
                "h_mat_Nhits_"+side,
                "h_mat_px_"+side,
                "h_mat_px_zoom_"+side,
                "h_mat_py_"+side,
                "h_mat_py_zoom_"+side,
                "h_mat_pz_"+side,
                "h_mat_pz_"+side+"_log0",
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
                "h_resol_mat_dpzrel_"+side,
                "h_resol_mat_dErel_"+side,
                "h_ratio_mat_pz_"+side+"_log0",
                "h_ratio_mat_pz_"+side+"_log1",
                "h_ratio_mat_E_"+side+"_log0",
                "h_ratio_mat_E_"+side+"_log1",
                
                "h_non_Nhits_"+side,
                "h_non_px_"+side,
                "h_non_px_zoom_"+side,
                "h_non_py_"+side,
                "h_non_py_zoom_"+side,
                "h_non_pz_"+side,
                "h_non_pz_"+side+"_log0",
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
                # "h_resol_non_dpzrel_"+side,
                # "h_resol_non_dErel_"+side,
                "h_ratio_non_pz_"+side+"_log0",
                "h_ratio_non_pz_"+side+"_log1",
                "h_ratio_non_E_"+side+"_log0",
                "h_ratio_non_E_"+side+"_log1",
                
                "h_sel_Nhits_"+side,
                "h_sel_px_"+side,
                "h_sel_px_zoom_"+side,
                "h_sel_py_"+side,
                "h_sel_py_zoom_"+side,
                "h_sel_pz_"+side,
                "h_sel_pz_"+side+"_log0",
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
                "h_resol_sel_dpzrel_"+side,
                "h_resol_sel_dErel_"+side,
                "h_ratio_sel_pz_"+side+"_log0",
                "h_ratio_sel_pz_"+side+"_log1",
                "h_ratio_sel_E_"+side+"_log0",
                "h_ratio_sel_E_"+side+"_log1",
                
                "h_cutflow_"+side,
             ]
      
      ### plot
      cnv = TCanvas("c","",700,500)
      cnv.SaveAs(outname+".pdf(")
      
      for name in hnames:
         print(name)
         
         htypes = {
            "ss":{"col":ROOT.kRed,     "leg":"Sig"},
            "sb":{"col":ROOT.kGreen+3, "leg":"Sig+Bkg"},
            "bb":{"col":ROOT.kBlue,    "leg":"Bkg"}
         }
         
         histos = {}
         
         cnv = TCanvas("c_"+name,"",700,500)
         if("cutflow" in name): cnv = TCanvas("c_"+name,"",900,500)
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
            if(("_E_" in name or "_pz_" in name) and "_ratio_" not in name and htyp=="ss"):
               name1 = name
               name1 = name1.replace("rec","tru").replace("sel","tru").replace("mat","tru").replace("non","tru")
               hname1 = name1+"_"+htyp
               # print("getting",htyp,name1)
               hist = files[htyp].Get(name1).Clone(hname1)
               if(not hist): print(name1,"is null")
               histos.update( {hname1:hist} )
               histos[hname1].SetDirectory(0)
               hmax = histos[hname1].GetMaximum() if(histos[hname1].GetMaximum()>hmax) else hmax
               histos[hname1].SetLineColor(ROOT.kBlack)
               histos[hname1].SetLineWidth(2)
               histos[hname1].SetLineStyle(2)
               legend.AddEntry(histos[hname1],"Tru","l")
            
            # print("getting",htyp,name)
            hist = files[htyp].Get(name).Clone(hname)
            if(not hist): print(name,"is null")
            histos.update( {hname:hist} )
            histos[hname].SetDirectory(0)
            hmax = histos[hname].GetMaximum() if(histos[hname].GetMaximum()>hmax) else hmax
            histos[hname].SetLineColor(attr["col"])
            histos[hname].SetFillColorAlpha(attr["col"],0.35)
            histos[hname].SetLineWidth(1)
            legend.AddEntry(histos[hname],attr["leg"],"f")
            if("cutflow" in name):
               hchop = hChopperUp(histos[hname],14)
               # hchop = hChopperUp(histos[hname])
               histos.update( {hname+"_chop":hchop} )
               
      
         ## normalise to nBX:
         for hname,h in histos.items():
            if("_ratio_" in name): continue
            nBX = -1
            if("ss" in hname): nBX = nBXs["ss"]
            if("sb" in hname): nBX = nBXs["sb"]
            if("bb" in hname): nBX = nBXs["bb"]
            h.Scale(1./nBX) 
      
         ## set max
         for hname,h in histos.items():
            h.SetMinimum(0.1)
            h.SetMaximum(hmax*1.1)
            if("Nhits" in name or "cutflow" in name): h.SetMaximum(hmax*2)
         
         ## draw truth just on energy plots
         if(("_E_" in name or "_pz_" in name) and "_ratio_" not in name):
            name1 = name
            name1 = name1.replace("rec","tru").replace("sel","tru").replace("mat","tru").replace("non","tru")
            htmp = histos[name1+"_ss"].Clone(name1+"_ss_tmp")
            htmp.GetXaxis().SetTitle( histos[name+"_sb"].GetXaxis().GetTitle() )
            htmp.GetYaxis().SetTitle( histos[name+"_sb"].GetYaxis().GetTitle() )
            htmp.Draw("hist")
            histos[name+"_sb"].Draw("hist same")
         else:
            if("cutflow" in name): histos[name+"_sb_chop"].Draw("hist")
            else:                  histos[name+"_sb"].Draw("hist")
         if("mat" not in name and "non" not in name):
            if("cutflow" in name): histos[name+"_bb_chop"].Draw("hist same")
            else:                  histos[name+"_bb"].Draw("hist same")
         if("cutflow" in name): histos[name+"_ss_chop"].Draw("hist same")
         else:                  histos[name+"_ss"].Draw("hist same")
         legend.Draw("sames")
         
         cnv.RedrawAxis()
         cnv.Update()
         cnv.SaveAs(outname+".pdf")
         cnv.SaveAs(pdfsdir+name.replace("h_","")+".pdf")
         
      cnv = TCanvas("c","",700,500)
      cnv.SaveAs(outname+".pdf)")

if __name__=="__main__":
    main()
		
		
