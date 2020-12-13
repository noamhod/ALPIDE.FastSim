#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
argus = parser.parse_args()
proc  = argus.p

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)
# storage = ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage = os.path.expandvars("$STORAGEDIR")

#############################################
def minmax(h1,h2,islogy=False,f=1.1):
   if(islogy): f=2.
   h1min = +1e10
   h2min = +1e10
   h1max = -1e10
   h2max = -1e10
   for b in range(h1.GetNbinsX()+1):
      h1min = h1.GetBinContent(b) if(h1.GetBinContent(b)<h1min) else h1min
      h1max = h1.GetBinContent(b) if(h1.GetBinContent(b)>h1max) else h1max
   for b in range(h2.GetNbinsX()+1):
      h2min = h2.GetBinContent(b) if(h2.GetBinContent(b)<h2min) else h2min
      h2max = h2.GetBinContent(b) if(h2.GetBinContent(b)>h2max) else h2max
   hmin = h1min if(h1min<h2min) else h2min
   hmax = h1max if(h1max>h2max) else h2max
   h1.SetMaximum(hmax*f)
   h2.SetMaximum(hmax*f)
   return hmin*f,hmax*f
   
def chopped(h,xmin,xmax):
   nbins0 = h.GetNbinsX()
   bmin = h.FindBin(xmin) if(xmin!=-9999) else 1
   bmax = h.FindBin(xmax) if(xmax!=-9999) else nbins0
   hxmin = h.GetXaxis().GetBinLowEdge(bmin)
   hxmax = h.GetXaxis().GetBinLowEdge(bmax)
   nbins = h.GetNbinsX() - bmin - (h.GetNbinsX()-bmax)
   print("nbins0=%g, nbins=%g" % (nbins0,nbins))
   hnew = TH1D(h.GetName()+"_copped",";"+h.GetXaxis().GetTitle()+";"+h.GetYaxis().GetTitle(),nbins,hxmin,hxmax)
   for b in range(bmin,bmax):
      y = h.GetBinContent(b)
      x = h.GetBinCenter(b)
      bnew = hnew.FindBin(x)
      hnew.SetBinContent(bnew,y)
   return hnew
   
#############################################

process = proc
print("process=",process)
print("storage=",storage)
tfile = TFile(storage+"/output/root/analysis_reco_"+process+".root","READ")
fn = storage+"/output/pdf/analysis_reco_"+process+"_"


leg_kin = TLegend(0.63,0.74,0.90,0.88)
leg_kin.SetFillStyle(4000) # will be transparent
leg_kin.SetFillColor(0)
leg_kin.SetTextFont(42)
leg_kin.SetBorderSize(0)

leg_eff = TLegend(0.30,0.43,0.77,0.6)
leg_eff.SetFillStyle(4000) # will be transparent
leg_eff.SetFillColor(0)
leg_eff.SetTextFont(42)
leg_eff.SetBorderSize(0)

leg_ntrks = TLegend(0.60,0.73,0.87,0.87)
leg_ntrks.SetFillStyle(4000) # will be transparent
leg_ntrks.SetFillColor(0)
leg_ntrks.SetTextFont(42)
leg_ntrks.SetBorderSize(0)

cnv = TCanvas("cnv_E","",500,500)
cnv.cd(1)
ROOT.gPad.SetLogy()
ROOT.gPad.SetTicks(1,1)
hmin,hmax = minmax(tfile.Get("h_E_sel"),tfile.Get("h_E_rec"),True)
hmin,hmax = minmax(tfile.Get("h_E_rec"),tfile.Get("h_E_sig"),True)
tfile.Get("h_E_sig").SetLineColor(ROOT.kGray+1)
tfile.Get("h_E_sig").SetFillColor(ROOT.kGray+1)
tfile.Get("h_E_sig").Draw("hist")
tfile.Get("h_E_rec").SetLineWidth(2)
tfile.Get("h_E_rec").SetLineColor(ROOT.kBlack)
tfile.Get("h_E_rec").Draw("hist same")
tfile.Get("h_E_sel").SetLineWidth(2)
tfile.Get("h_E_sel").SetLineColor(ROOT.kRed)
tfile.Get("h_E_sel").Draw("hist same")
leg_kin.AddEntry(tfile.Get("h_E_sig"),"Signal","f")
leg_kin.AddEntry(tfile.Get("h_E_rec"),"Reconstructed","l")
leg_kin.AddEntry(tfile.Get("h_E_sel"),"Selected","l")
leg_kin.Draw("same")
cnv.SaveAs(fn+"selection_energy.pdf")


cnv = TCanvas("cnv_kin","",1000,1000)
cnv.Divide(2,2)
cnv.cd(1)
ROOT.gPad.SetLogy()
ROOT.gPad.SetTicks(1,1)
# hmin,hmax = minmax(tfile.Get("h_E_sel"),tfile.Get("h_E_rec"),True)
# hmin,hmax = minmax(tfile.Get("h_E_rec"),tfile.Get("h_E_sig"),True)
tfile.Get("h_E_sig").Draw("hist")
tfile.Get("h_E_rec").Draw("hist same")
tfile.Get("h_E_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_E_sel").Draw("hist same")
# leg_kin.AddEntry(tfile.Get("h_E_sig"),"Signal","l")
# leg_kin.AddEntry(tfile.Get("h_E_rec"),"Reconstructed","l")
# leg_kin.AddEntry(tfile.Get("h_E_sel"),"Selected","l")
leg_kin.Draw("same")
cnv.cd(2)
ROOT.gPad.SetLogy()
ROOT.gPad.SetTicks(1,1)
hmin,hmax = minmax(tfile.Get("h_pz_sel"),tfile.Get("h_pz_rec"),True)
hmin,hmax = minmax(tfile.Get("h_pz_rec"),tfile.Get("h_pz_sig"),True)
tfile.Get("h_pz_sig").Draw("hist")
tfile.Get("h_pz_rec").Draw("hist same")
tfile.Get("h_pz_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_pz_sel").Draw("hist same")
leg_kin.Draw("same")
cnv.cd(3)
ROOT.gPad.SetLogy()
ROOT.gPad.SetTicks(1,1)
hmin,hmax = minmax(tfile.Get("h_px_zoom_sel"),tfile.Get("h_px_zoom_rec"),True)
hmin,hmax = minmax(tfile.Get("h_px_zoom_rec"),tfile.Get("h_px_zoom_sig"),True)
tfile.Get("h_px_zoom_sig").Draw("hist")
tfile.Get("h_px_zoom_rec").Draw("hist same")
tfile.Get("h_px_zoom_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_px_zoom_sel").Draw("hist same")
leg_kin.Draw("same")
cnv.cd(4)
ROOT.gPad.SetLogy()
ROOT.gPad.SetTicks(1,1)
hmin,hmax = minmax(tfile.Get("h_py_zoom_sel"),tfile.Get("h_py_zoom_rec"),True)
hmin,hmax = minmax(tfile.Get("h_py_zoom_rec"),tfile.Get("h_py_zoom_sig"),True)
tfile.Get("h_py_zoom_sig").Draw("hist")
tfile.Get("h_py_zoom_rec").Draw("hist same")
tfile.Get("h_py_zoom_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_py_zoom_sel").Draw("hist same")
leg_kin.Draw("same")
cnv.SaveAs(fn+"selection_kinematics.pdf")


# cnv = TCanvas("cnv","",1000,500)
# cnv.Divide(2,1)
# cnv.cd(1)
# ROOT.gPad.SetTicks(1,1)
# tfile.Get("h_E_tru_eff").SetMinimum(0)
# tfile.Get("h_E_tru_eff").SetMaximum(1.05)
# tfile.Get("h_E_tru_eff").SetLineColor(ROOT.kBlack)
# tfile.Get("h_E_tru_eff").SetMarkerColor(ROOT.kBlack)
# tfile.Get("h_E_tru_eff").SetMarkerStyle(20)
# tfile.Get("h_E_tru_eff").Draw("ep")
# tfile.Get("h_E_tru_eff_sel").SetMinimum(0)
# tfile.Get("h_E_tru_eff_sel").SetMaximum(1.05)
# tfile.Get("h_E_tru_eff_sel").SetLineColor(ROOT.kGray+1)
# tfile.Get("h_E_tru_eff_sel").SetMarkerColor(ROOT.kGray+1)
# tfile.Get("h_E_tru_eff_sel").SetMarkerStyle(24)
# tfile.Get("h_E_tru_eff_sel").Draw("ep same")
# leg_eff.AddEntry(tfile.Get("h_E_tru_eff"),"Reconstruction","epl")
# leg_eff.AddEntry(tfile.Get("h_E_tru_eff_sel"),"Selected","epl")
# leg_eff.Draw("same")
# cnv.cd(2)
# ROOT.gPad.SetTicks(1,1)
# tfile.Get("h_E_tru_eff_acc").SetMinimum(0)
# tfile.Get("h_E_tru_eff_acc").SetMaximum(1.05)
# tfile.Get("h_E_tru_eff_acc").SetLineColor(ROOT.kBlack)
# tfile.Get("h_E_tru_eff_acc").SetMarkerColor(ROOT.kBlack)
# tfile.Get("h_E_tru_eff_acc").SetMarkerStyle(20)
# tfile.Get("h_E_tru_eff_acc").Draw("ep")
# tfile.Get("h_E_tru_eff_acc_sel").SetMinimum(0)
# tfile.Get("h_E_tru_eff_acc_sel").SetMaximum(1.05)
# tfile.Get("h_E_tru_eff_acc_sel").SetLineColor(ROOT.kGray+1)
# tfile.Get("h_E_tru_eff_acc_sel").SetMarkerColor(ROOT.kGray+1)
# tfile.Get("h_E_tru_eff_acc_sel").SetMarkerStyle(24)
# tfile.Get("h_E_tru_eff_acc_sel").Draw("ep same")
# leg_eff.Draw("same")
# cnv.SaveAs(fn+"selection_eff.pdf")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetTicks(1,1)
tfile.Get("h_E_tru_eff").SetMinimum(0)
tfile.Get("h_E_tru_eff").SetMaximum(1.05)
tfile.Get("h_E_tru_eff").SetLineColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerStyle(20)
# tfile.Get("h_E_tru_eff").Draw("ep")
tfile.Get("h_E_tru_eff_sel").SetMinimum(0)
tfile.Get("h_E_tru_eff_sel").SetMaximum(1.05)
tfile.Get("h_E_tru_eff_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_sel").SetMarkerColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_sel").SetMarkerStyle(24)
# tfile.Get("h_E_tru_eff_sel").Draw("ep same")
rp = ROOT.TRatioPlot(tfile.Get("h_E_tru_eff_sel"),tfile.Get("h_E_tru_eff"))
rp.SetH1DrawOpt("ep")
rp.SetH2DrawOpt("ep")
rp.SetGraphDrawOpt("ALX")
ratioMax = 1.008 #maximum value on ratio plot Y axis
ratioMin = 0.965 #minimum value on ratio plot Y axis
rp.Draw("ep0 nohide")
rp.GetLowerRefGraph().SetMaximum(ratioMax)
rp.GetLowerRefGraph().SetMinimum(ratioMin)
# rp.SetLeftMargin(0.13)
rp.GetLowerRefGraph().GetYaxis().SetTitleOffset(1.5)
# rp.GetLowerRefGraph().GetYaxis().SetLabelSize(0.025)
rp.GetLowYaxis().SetNdivisions(506)
rp.SetSeparationMargin(0.0)
rp.GetLowerRefYaxis().SetTitle("Sel/Rec")
leg_eff.AddEntry(tfile.Get("h_E_tru_eff"),"Reconstruction","epl")
leg_eff.AddEntry(tfile.Get("h_E_tru_eff_sel"),"Selected","epl")
leg_eff.Draw("same")
leg_eff.Draw("same")
cnv.SaveAs(fn+"selection_eff_only.pdf")


cnv = TCanvas("cnv","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetTicks(1,1)
tfile.Get("h_E_tru_eff").SetMinimum(0)
tfile.Get("h_E_tru_eff").SetMaximum(1.05)
tfile.Get("h_E_tru_eff").SetLineColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerStyle(20)
# tfile.Get("h_E_tru_eff").Draw("ep")
tfile.Get("h_E_tru_eff_sel").SetMinimum(0)
tfile.Get("h_E_tru_eff_sel").SetMaximum(1.05)
tfile.Get("h_E_tru_eff_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_sel").SetMarkerColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_sel").SetMarkerStyle(24)
# tfile.Get("h_E_tru_eff_sel").Draw("ep same")
rp = ROOT.TRatioPlot(tfile.Get("h_E_tru_eff_sel"),tfile.Get("h_E_tru_eff"))
rp.SetH1DrawOpt("ep")
rp.SetH2DrawOpt("ep")
rp.SetGraphDrawOpt("ALX")
ratioMax = 1.008 #maximum value on ratio plot Y axis
ratioMin = 0.965 #minimum value on ratio plot Y axis
rp.Draw("ep0 nohide")
rp.GetLowerRefGraph().SetMaximum(ratioMax)
rp.GetLowerRefGraph().SetMinimum(ratioMin)
# rp.SetLeftMargin(0.13)
rp.GetLowerRefGraph().GetYaxis().SetTitleOffset(1.5)
# rp.GetLowerRefGraph().GetYaxis().SetLabelSize(0.025)
rp.GetLowYaxis().SetNdivisions(506)
rp.SetSeparationMargin(0.0)
rp.GetLowerRefYaxis().SetTitle("Sel/Rec")
# leg_eff.AddEntry(tfile.Get("h_E_tru_eff"),"Reconstruction","epl")
# leg_eff.AddEntry(tfile.Get("h_E_tru_eff_sel"),"Selected","epl")
leg_eff.Draw("same")

cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
tfile.Get("h_E_tru_eff_acc").SetMinimum(0)
tfile.Get("h_E_tru_eff_acc").SetMaximum(1.05)
tfile.Get("h_E_tru_eff_acc").SetLineColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff_acc").SetMarkerColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff_acc").SetMarkerStyle(20)
# tfile.Get("h_E_tru_eff_acc").Draw("ep")
tfile.Get("h_E_tru_eff_acc_sel").SetMinimum(0)
tfile.Get("h_E_tru_eff_acc_sel").SetMaximum(1.05)
tfile.Get("h_E_tru_eff_acc_sel").SetLineColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_acc_sel").SetMarkerColor(ROOT.kGray+1)
tfile.Get("h_E_tru_eff_acc_sel").SetMarkerStyle(24)
# tfile.Get("h_E_tru_eff_acc_sel").Draw("ep same")
rp_acc = ROOT.TRatioPlot(tfile.Get("h_E_tru_eff_acc_sel"),tfile.Get("h_E_tru_eff_acc"))
rp_acc.SetH1DrawOpt("ep")
rp_acc.SetH2DrawOpt("ep")
rp_acc.SetGraphDrawOpt("ALX")
ratioMax = 1.008 #maximum value on ratio plot Y axis
ratioMin = 0.965 #minimum value on ratio plot Y axis
rp_acc.Draw("ep0 nohide")
rp_acc.GetLowerRefGraph().SetMaximum(ratioMax)
rp_acc.GetLowerRefGraph().SetMinimum(ratioMin)
rp_acc.GetLowerRefGraph().GetYaxis().SetTitleOffset(1.55)
# rp_acc.GetLowerRefGraph().GetYaxis().SetLabelSize(0.025)
rp_acc.GetLowYaxis().SetNdivisions(506)
rp_acc.SetSeparationMargin(0.0)
rp_acc.GetLowerRefYaxis().SetTitle("Sel/Rec")
leg_eff.Draw("same")
cnv.SaveAs(fn+"selection_eff.pdf")




cnv = TCanvas("cnv_ntrks","",500,500)
cnv.cd()
ROOT.gPad.SetTicks(1,1)
hsel = chopped(tfile.Get("h_ntrks_sel"),20,200)
hrec = chopped(tfile.Get("h_ntrks_rec"),20,200)
hsig = chopped(tfile.Get("h_ntrks_sig"),20,200)
hmin,hmax = minmax(hsel,hrec,False)
hmin,hmax = minmax(hrec,hsig,False)
hsel.SetLineColor(ROOT.kGray+1)
hsel.SetLineWidth(2)
hrec.SetLineColor(ROOT.kBlack)
hrec.SetLineWidth(2)
hsig.SetLineColor(ROOT.kRed)
hsig.SetFillColor(ROOT.kRed)
hsig.Draw("hist")
hrec.Draw("hist same")
hsel.Draw("hist same")
leg_ntrks.AddEntry(hsig,"Signal","f")
leg_ntrks.AddEntry(hrec,"Reconstructed","l")
leg_ntrks.AddEntry(hsel,"Selected","l")
leg_ntrks.Draw("same")
cnv.SaveAs(fn+"selection_ntrks.pdf")



xmineff = 1.5
xmaxeff = 15.5

cnv = TCanvas("cnv","",500,500)
cnv.cd()
ROOT.gPad.SetTicks(1,1)
tfile.Get("h_E_tru_eff").SetMinimum(0)
tfile.Get("h_E_tru_eff").SetMaximum(1.05)
tfile.Get("h_E_tru_eff").SetLineColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerColor(ROOT.kBlack)
tfile.Get("h_E_tru_eff").SetMarkerStyle(20)

# sfturnon = '[0]*ROOT::Math::erfc([1]*x-[2]) + [3]/([4]+ROOT::Math::exp(-[4]*(x-1)))'
sfturnon = '[0]*(1+ROOT::Math::erfc([1]*(x-1))) + [2]/([2]+ROOT::Math::exp(-[2]*(x-1)))'
positive = []

fturnon_guess = TF1("turnon_guess", sfturnon,xmineff,xmaxeff)
fturnon_guess.SetLineColor(ROOT.kBlue)
fturnon_guess.SetLineWidth(2)
fturnon_guess.SetParameter(0, 1)
fturnon_guess.SetParameter(1, 1)
fturnon_guess.SetParameter(2, 1)
# fturnon_guess.SetParameter(3, 1)
# fturnon_guess.SetParameter(4, 1)
# fturnon_guess.SetParameter(5, 1)
# fturnon_guess.SetParameter(6, 1)


fturnon = TF1("turnon", sfturnon, xmineff,xmaxeff)
for pos in positive: fturnon.SetParLimits(pos,0,1000)
fturnon.SetLineWidth(2)
fturnon.SetLineColor(ROOT.kRed)

for n in range(fturnon.GetNpar()): fturnon.SetParameter(n,fturnon_guess.GetParameter(n))

res = tfile.Get("h_E_tru_eff").Fit(fturnon,"MERS")
chi2dof = fturnon.GetChisquare()/fturnon.GetNDF() if(fturnon.GetNDF()>0) else -1
print("Turnon: chi2/Ndof=",chi2dof)
tfile.Get("h_E_tru_eff").Draw("ep")
fturnon.Draw("same")
fturnon_guess.Draw("same")
cnv.SaveAs(fn+"E_tru_eff_turnonfit.pdf")

tfileout = TFile(fn.replace("pdf","root")+"accxeff.root","RECREATE")
tfileout.cd()
fturnon.Write()
tfileout.Write()
tfileout.Close()




