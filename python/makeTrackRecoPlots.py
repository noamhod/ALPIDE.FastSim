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
   legend1 = TLegend(0.77,0.80,0.92,0.92); ### with the mean
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


def isTruth(name,truthvars):
   for var in truthvars:
      if(var in name): return True
   return False


def GetMeanSigma(fitf,name):
   xmin = fitf.GetXmin()
   xmax = fitf.GetXmax()
   mean = fitf.Mean(xmin,xmax)
   sigma = ROOT.TMath.Sqrt(fitf.Variance(xmin,xmax))
   xmin = mean-5*sigma
   xmax = mean+5*sigma
   mean = fitf.Mean(xmin,xmax)
   sigma = ROOT.TMath.Sqrt(fitf.Variance(xmin,xmax))
   print("%s: mean=%.5f, sigma=%.5f" % (name,mean,sigma))
   return mean,sigma


def TrippleGausFitResE(hResE, cname, fsinglename, fallname):
   cnv = TCanvas(cname,"",700,500)
   cnv.SetTicks(1,1)
   g1 = TF1("g1", "gaus", -0.03,+0.03); g1.SetLineColor(ROOT.kViolet)
   g2 = TF1("g2", "gaus", -0.03,+0.03); g2.SetLineColor(ROOT.kGreen+2)
   g3 = TF1("g3", "gaus", -0.03,+0.03); g3.SetLineColor(ROOT.kYellow+2)
   hResE.Fit(g1,"EMRS")
   hResE.Fit(g2,"EMRS")
   hResE.Fit(g3,"EMRS")
   fite = TF1("fite", "gaus(0)+gaus(3)+gaus(6)", -0.03,+0.03)
   fite.SetLineColor(ROOT.kBlack)
   fite.SetLineWidth(2)
   fite.SetParameter(0,g1.GetParameter(0))
   fite.SetParameter(1,g1.GetParameter(1))
   fite.SetParameter(2,g1.GetParameter(2))
   fite.SetParameter(3,g2.GetParameter(0))
   fite.SetParameter(4,g2.GetParameter(1))
   fite.SetParameter(5,g2.GetParameter(2))
   fite.SetParameter(6,g3.GetParameter(0))
   fite.SetParameter(7,g3.GetParameter(1))
   fite.SetParameter(8,g3.GetParameter(2))
   res = hResE.Fit(fite,"EMRS")
   chi2dof = fite.GetChisquare()/fite.GetNDF() if(fite.GetNDF()>0) else -1
   print("Res(E rec:tru) chi2/Ndof=",chi2dof)
   hResE.Draw("hist")
   fite.Draw("same")
   mean,sigma = GetMeanSigma(fite,"E")
   s = ROOT.TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.04)
   s.DrawLatex(0.7,0.85,ROOT.Form("Mean=%.6f" % (mean)))
   s.DrawLatex(0.7,0.78,ROOT.Form("#sigma=%.4f" % (sigma)))
   s.DrawLatex(0.7,0.71,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
   # LUXE(0.18,0.85,ROOT.kBlack)
   cnv.SaveAs(fsinglename)
   cnv.SaveAs(fallname)

def DoubleGausFitResE(hResE, cname, fsinglename, fallname):
   cnv = TCanvas(cname,"",700,500)
   cnv.SetTicks(1,1)
   g1 = TF1("g1", "gaus", -0.03,+0.03); g1.SetLineColor(ROOT.kViolet)
   g2 = TF1("g2", "gaus", -0.03,+0.03); g2.SetLineColor(ROOT.kGreen+2)
   hResE.Fit(g1,"EMRS")
   hResE.Fit(g2,"EMRS")
   fite = TF1("fite", "gaus(0)+gaus(3)", -0.03,+0.03)
   fite.SetLineColor(ROOT.kBlack)
   fite.SetLineWidth(2)
   fite.SetParameter(0,g1.GetParameter(0))
   fite.SetParameter(1,g1.GetParameter(1))
   fite.SetParameter(2,g1.GetParameter(2))
   fite.SetParameter(3,g2.GetParameter(0))
   fite.SetParameter(4,g2.GetParameter(1))
   fite.SetParameter(5,g2.GetParameter(2))
   res = hResE.Fit(fite,"EMRS")
   chi2dof = fite.GetChisquare()/fite.GetNDF() if(fite.GetNDF()>0) else -1
   print("Res(E rec:tru) chi2/Ndof=",chi2dof)
   hResE.Draw("hist")
   fite.Draw("same")
   mean,sigma = GetMeanSigma(fite,"E")
   s = ROOT.TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.04)
   s.DrawLatex(0.7,0.85,ROOT.Form("Mean=%.6f" % (mean)))
   s.DrawLatex(0.7,0.78,ROOT.Form("#sigma=%.4f" % (sigma)))
   s.DrawLatex(0.7,0.71,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
   # LUXE(0.18,0.85,ROOT.kBlack)
   cnv.SaveAs(fsinglename)
   cnv.SaveAs(fallname)
   
def SingleGausFitResE(hResE, cname, fsinglename, fallname):
   cnv = TCanvas(cname,"",700,500)
   cnv.SetTicks(1,1)
   fite = TF1("g1", "gaus", -0.03,+0.03)
   fite.SetLineColor(ROOT.kBlack)
   fite.SetLineWidth(2)
   hResE.Fit(fite,"EMRS")
   res = hResE.Fit(fite,"EMRS")
   chi2dof = fite.GetChisquare()/fite.GetNDF() if(fite.GetNDF()>0) else -1
   print("Res(E rec:tru) chi2/Ndof=",chi2dof)
   hResE.Draw("hist")
   fite.Draw("same")
   mean  = fite.GetParameter("Mean")
   sigma = fite.GetParameter("Sigma")
   s = ROOT.TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.04)
   s.DrawLatex(0.7,0.85,ROOT.Form("Mean=%.6f" % (mean)))
   s.DrawLatex(0.7,0.78,ROOT.Form("#sigma=%.4f" % (sigma)))
   s.DrawLatex(0.7,0.71,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
   # LUXE(0.18,0.85,ROOT.kBlack)
   cnv.SaveAs(fsinglename)
   cnv.SaveAs(fallname)


def GetHminmax(h):
   hmin = +1e20
   hmax = -1e20
   for b in range(1,h.GetNbinsX()):
      y = h.GetBinContent(b)
      if(y<hmin and y>0): hmin = y
      if(y>hmax):         hmax = y
   return hmin,hmax


######################################################
######################################################
######################################################

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
   
   outname = "signal_and_background_tracking_variable_"+process+"_"+sample.replace("/","_")
   
   sides = ["Eside","Pside"] if(sidepe==None) else [sidepe]
   
   for side in sides:
      
      fnames = {
         "ss":"../data/root/rec/"+process+"/"+sample+"s/rec_"+process+"__"+side+".root",
         "sb":"../data/root/rec/"+process+"/"+sample+"sb/rec_"+process+"__"+side+".root",
         "bb":"../data/root/rec/"+process+"/"+sample+"b/rec_"+process+"__"+side+".root"
      }
      if("linearity" not in sample):
         fnames["bb"] = "../data/root/rec/"+process+"/bkg/rec_"+process+"__"+side+".root"

      print(fnames)
      
      files = {}
      nBXs  = {}
      for ftype,fname in fnames.items():
         files.update({ftype:TFile(fname,"READ")})
         nBXs.update({ftype:files[ftype].Get("h_nBX_"+side).GetBinContent(1)})
      
      plottruth = ["_E_", "_pz_", "_px_", "_py_"]
      
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
                "h_resol_rec_dx1_"+side,
                "h_resol_rec_dx4_"+side,
                "h_resol_rec_dy1_"+side,
                "h_resol_rec_dy4_"+side,
                "h_ratio_rec_pz_"+side+"_log0",
                "h_ratio_rec_pz_"+side+"_log1",
                "h_ratio_rec_pz_"+side+"_log2",
                "h_ratio_rec_pz_"+side+"_log3",
                "h_ratio_rec_E_"+side+"_log0",
                "h_ratio_rec_E_"+side+"_log1",
                "h_ratio_rec_E_"+side+"_log2",
                "h_ratio_rec_E_"+side+"_log3",
                "h_eff_rec_pz_"+side+"_log0",
                "h_eff_rec_pz_"+side+"_log1",
                "h_eff_rec_pz_"+side+"_log2",
                "h_eff_rec_pz_"+side+"_log3",
                "h_eff_rec_E_"+side+"_log0",
                "h_eff_rec_E_"+side+"_log1",
                "h_eff_rec_E_"+side+"_log2",
                "h_eff_rec_E_"+side+"_log3",
                
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
                "h_resol_mat_dx1_"+side,
                "h_resol_mat_dx4_"+side,
                "h_resol_mat_dy1_"+side,
                "h_resol_mat_dy4_"+side,
                "h_ratio_mat_pz_"+side+"_log0",
                "h_ratio_mat_pz_"+side+"_log1",
                "h_ratio_mat_pz_"+side+"_log2",
                "h_ratio_mat_pz_"+side+"_log3",
                "h_ratio_mat_E_"+side+"_log0",
                "h_ratio_mat_E_"+side+"_log1",
                "h_ratio_mat_E_"+side+"_log2",
                "h_ratio_mat_E_"+side+"_log3",
                "h_eff_mat_pz_"+side+"_log0",
                "h_eff_mat_pz_"+side+"_log1",
                "h_eff_mat_pz_"+side+"_log2",
                "h_eff_mat_pz_"+side+"_log3",
                "h_eff_mat_E_"+side+"_log0",
                "h_eff_mat_E_"+side+"_log1",
                "h_eff_mat_E_"+side+"_log2",
                "h_eff_mat_E_"+side+"_log3",
                
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
                "h_ratio_non_pz_"+side+"_log2",
                "h_ratio_non_pz_"+side+"_log3",
                "h_ratio_non_E_"+side+"_log0",
                "h_ratio_non_E_"+side+"_log1",
                "h_ratio_non_E_"+side+"_log2",
                "h_ratio_non_E_"+side+"_log3",
                "h_eff_non_pz_"+side+"_log0",
                "h_eff_non_pz_"+side+"_log1",
                "h_eff_non_pz_"+side+"_log2",
                "h_eff_non_pz_"+side+"_log3",
                "h_eff_non_E_"+side+"_log0",
                "h_eff_non_E_"+side+"_log1",
                "h_eff_non_E_"+side+"_log2",
                "h_eff_non_E_"+side+"_log3",
                
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
                "h_resol_sel_dx1_"+side,
                "h_resol_sel_dx4_"+side,
                "h_resol_sel_dy1_"+side,
                "h_resol_sel_dy4_"+side,
                "h_ratio_sel_pz_"+side+"_log0",
                "h_ratio_sel_pz_"+side+"_log1",
                "h_ratio_sel_pz_"+side+"_log2",
                "h_ratio_sel_pz_"+side+"_log3",
                "h_ratio_sel_E_"+side+"_log0",
                "h_ratio_sel_E_"+side+"_log1",
                "h_ratio_sel_E_"+side+"_log2",
                "h_ratio_sel_E_"+side+"_log3",
                "h_eff_sel_pz_"+side+"_log0",
                "h_eff_sel_pz_"+side+"_log1",
                "h_eff_sel_pz_"+side+"_log2",
                "h_eff_sel_pz_"+side+"_log3",
                "h_eff_sel_E_"+side+"_log0",
                "h_eff_sel_E_"+side+"_log1",
                "h_eff_sel_E_"+side+"_log2",
                "h_eff_sel_E_"+side+"_log3",
                
                "h_cutflow_"+side,
             ]
      
      ### plot
      cnv = TCanvas("c","",700,500)
      cnv.SaveAs(outname+".pdf(")
      
      for name in hnames:
         # print(name)
         
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
         if("ratio" not in name and "eff" not in name and "resol" not in name): ROOT.gPad.SetLogy()
         
         LUXELabel(0.2,0.85,"TDR")
         legend = LegendMaker()
         
         for htyp,attr in htypes.items():
            hname = name+"_"+htyp
            
            ## draw truth just on energy plots
            istruth = isTruth(name,plottruth)
            if(istruth and "_ratio_" not in name and "_eff_" not in name and htyp=="ss"):
               name1 = name
               name1 = name1.replace("rec","tru").replace("sel","tru").replace("mat","tru").replace("non","tru")
               hname1 = name1+"_"+htyp
               # print("getting",htyp,name1)
               hist = files[htyp].Get(name1).Clone(hname1)
               if(not hist): print(name1,"is null")
               histos.update( {hname1:hist} )
               histos[hname1].SetDirectory(0)
               histos[hname1].SetLineColor(ROOT.kBlack)
               histos[hname1].SetLineWidth(2)
               histos[hname1].SetLineStyle(2)
               legend.AddEntry(histos[hname1],"Tru","l")
            
            # print("getting",htyp,name)
            hist = files[htyp].Get(name).Clone(hname)
            if(not hist): print(name,"is null")
            histos.update( {hname:hist} )
            histos[hname].SetDirectory(0)
            histos[hname].SetLineColor(attr["col"])
            histos[hname].SetFillColorAlpha(attr["col"],0.35)
            histos[hname].SetLineWidth(1)
            legend.AddEntry(histos[hname],attr["leg"],"f")
            if("cutflow" in name):
               hchop = hChopperUp(histos[hname],15)
               histos.update( {hname+"_chop":hchop} )
         
      
         ## normalise to nBX:
         hmax = -1e20
         for hname,h in histos.items():
            if("_ratio_" not in hname and "_eff_" not in hname):
               nBX = -1
               if("ss" in hname): nBX = nBXs["ss"]
               if("sb" in hname): nBX = nBXs["sb"]
               if("bb" in hname): nBX = nBXs["bb"]
               h.Scale(1./nBX)
               if("_resol_" in hname and "bb" not in hname):
                  h.Scale(1./h.Integral())
                  h.GetYaxis().SetTitle( h.GetYaxis().GetTitle()+" [Normalised]" )
            hmax = h.GetMaximum() if(h.GetMaximum()>hmax) else hmax

      
         ## set max
         for hname,h in histos.items():
            if("_ratio_" in hname or "_eff_" in hname or "_resol_" in hname):
               h.SetMinimum(0.0)
               h.SetMaximum(hmax*1.25)
            else:
               h.SetMinimum(0.1)
               if("Nhits" in name or "cutflow" in name): h.SetMaximum(hmax*2)
               else:                                     h.SetMaximum(hmax*1.25)

         
         ## draw truth just on energy plots
         if(istruth and "_ratio_" not in name and "_eff_" not in name):
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
      
      ### resolution fits:
      fitnames = [
         "h_resol_rec_dErel_"+side,
         "h_resol_mat_dErel_"+side,
         "h_resol_sel_dErel_"+side,
      ]
      for name in fitnames:
         histos = {}
         for htyp,attr in htypes.items():
            if(htyp=="bb"): continue
            hname = name+"_"+htyp
            hist = files[htyp].Get(name).Clone(hname)
            if(not hist): print(name,"is null")
            histos.update( {hname:hist} )
            histos[hname].SetDirectory(0)
            # hmax = histos[hname].GetMaximum() if(histos[hname].GetMaximum()>hmax) else hmax
            histos[hname].SetLineColor(attr["col"])
            histos[hname].SetFillColorAlpha(attr["col"],0.35)
            histos[hname].SetLineWidth(1)
            
            ## normalise to nBX:
            h.Scale(1./nBXs[htyp])
            
            hmin,hmax = GetHminmax(histos[hname])
            # print(name,histos[hname].GetMaximum(),hmax)
            histos[hname].SetMinimum(0)
            histos[hname].SetMaximum(1.1*hmax)
            # TrippleGausFitResE(histos[hname],"cnv_resE",outname+".pdf",pdfsdir+name.replace("h_","fit_")+".pdf")
            DoubleGausFitResE(histos[hname],"cnv_resE",outname+".pdf",pdfsdir+name.replace("h_","fit_")+".pdf")
            # SingleGausFitResE(histos[hname],"cnv_resE",outname+".pdf",pdfsdir+name.replace("h_","fit_")+".pdf")
      
      cnv = TCanvas("c","",700,500)
      cnv.SaveAs(outname+".pdf)")

if __name__=="__main__":
    main()
		
		
