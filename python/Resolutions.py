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
ROOT.gStyle.SetPadLeftMargin(0.13)

#############################################

def GetLogBinning(nbins,xmin,xmax):
   logmin  = math.log10(xmin)
   logmax  = math.log10(xmax)
   logbinwidth = (logmax-logmin)/nbins
   # Bin edges in GeV
   xbins = [xmin,] #the lowest edge first
   for i in range(1,nbins+1):
      xbins.append( ROOT.TMath.Power( 10,(logmin + i*logbinwidth) ) )
   arrxbins = array.array("d", xbins)
   return nbins, arrxbins
   
def GetMeany(h,bx):
   av = 0
   norm = 0
   for by in range(h.GetNbinsY()+1):
      av += h.GetYaxis().GetBinCenter(by)*h.GetBinContent(bx,by)
      norm += h.GetBinContent(bx,by)
   if(norm==0): return -99999,-99999
   return av/norm,norm
   
def GetRMSy(h,bx):
   av,norm = GetMeany(h,bx)
   if(av==-99999 and norm==-99999): return -99999,-99999
   rms = 0
   for by in range(h.GetNbinsY()+1):
      n = h.GetBinContent(bx,by)
      y = h.GetYaxis().GetBinCenter(by)
      rms += n*(y-av)*(y-av)
   rms = math.sqrt(rms/norm)
   return rms,av
   
def GetSigma(fitf,name):
   xmin = fitf.GetXmin()
   xmax = fitf.GetXmax()
   mean = fitf.Mean(xmin,xmax)
   sigma = ROOT.TMath.Sqrt(fitf.Variance(xmin,xmax))
   xmin = mean-5*sigma
   xmax = mean+5*sigma
   mean = fitf.Mean(xmin,xmax)
   sigma = ROOT.TMath.Sqrt(fitf.Variance(xmin,xmax))
   print("%s: mean=%.5f, sigma=%.5f" % (name,mean,sigma))
   return sigma

#############################################

process = proc
tfile = TFile("../output/root/analysis_reco_"+process+".root","READ")
fn = "../output/pdf/analysis_reco_"+process+"_"

hx = tfile.Get("h_dx_vs_xtru_recvstru_logbins")
hy = tfile.Get("h_dy_vs_ytru_recvstru_logbins")

cnv = TCanvas("cnv_x","",500,500)
cnv.SaveAs(fn+"x_tests.pdf(")
nxbinstmp,xbinstmp = GetLogBinning(50,1,70)
xbins = []
for n in reversed(range(len(xbinstmp))): xbins.append(-xbinstmp[n])
for n in range(len(xbinstmp)):           xbins.append(+xbinstmp[n])
hxRMS = TH1D("hxRMS","",len(xbins)-1,array.array("d",xbins))
for bx in range(hx.GetNbinsX()+1):
   cnv = TCanvas("cnv_"+str(bx),"",500,500)
   htmp = TH1D("htmp_"+str(bx),"",hx.GetNbinsY(),hx.GetYaxis().GetXmin(),hx.GetYaxis().GetXmax())
   for by in range(hx.GetNbinsY()+1): htmp.SetBinContent(by,hx.GetBinContent(bx,by))
   x    = hx.GetXaxis().GetBinCenter(bx)
   mean = htmp.GetMean()
   rms  = htmp.GetRMS()
   rmserr = htmp.GetRMSError()
   htmp.Draw()
   rmstmp,meantmp = GetRMSy(hx,bx)
   s = ROOT.TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.025)
   s.DrawLatex(0.15,0.85,ROOT.Form("x=%.3f" % (x)))
   s.DrawLatex(0.15,0.80,ROOT.Form("TH1: mean=%.6f, rms=%.4f" % (mean,rms)))
   s.DrawLatex(0.15,0.75,ROOT.Form("TMP: mean=%.6f, rms=%.4f" % (meantmp,rmstmp)))
   cnv.SaveAs(fn+"x_tests.pdf")
   hxRMS.SetBinContent(bx,rms)
   hxRMS.SetBinError(bx,rmserr)
   # print("x=%g, mean=%g, rms=%g" % (x,mean,rms))
cnv = TCanvas("cnv_b","",500,500)
cnv.SaveAs(fn+"x_tests.pdf)")


cnv = TCanvas("cnv_y","",500,500)
cnv.SaveAs(fn+"y_tests.pdf(")

nybinstmp,ybinstmp = GetLogBinning(50,0.001,0.4)
ybins = []
for n in reversed(range(len(ybinstmp))): ybins.append(-ybinstmp[n])
for n in range(len(ybinstmp)):           ybins.append(+ybinstmp[n])  
hyRMS = TH1D("hyRMS","",len(xbins)-1,array.array("d",ybins))
for bx in range(hy.GetNbinsX()+1):
   cnv = TCanvas("cnv_"+str(bx),"",500,500)
   htmp = TH1D("htmp_"+str(bx),"",hy.GetNbinsY(),hy.GetYaxis().GetXmin(),hy.GetYaxis().GetXmax())
   for by in range(hy.GetNbinsY()+1): htmp.SetBinContent(by,hy.GetBinContent(bx,by))
   y    = hy.GetXaxis().GetBinCenter(bx)
   mean = htmp.GetMean()
   rms  = htmp.GetRMS()
   rmserr = htmp.GetRMSError()
   hyRMS.SetBinContent(bx,rms)
   hyRMS.SetBinError(bx,rmserr)
   htmp.Draw()
   rmstmp,meantmp = GetRMSy(hy,bx)
   s = ROOT.TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.025)
   s.DrawLatex(0.15,0.85,ROOT.Form("y=%.3f" % (y)))
   s.DrawLatex(0.15,0.80,ROOT.Form("TH1: mean=%.6f, rms=%.4f" % (mean,rms)))
   s.DrawLatex(0.15,0.75,ROOT.Form("TMP: mean=%.6f, rms=%.4f" % (meantmp,rmstmp)))
   cnv.SaveAs(fn+"y_tests.pdf")
   # print("y=%g, mean=%g, rms=%g" % (y,mean,rms))
cnv = TCanvas("cnv_b","",500,500)
cnv.SaveAs(fn+"y_tests.pdf)")


cnv = TCanvas("cnv_res_xy_vs_xy_log_log","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetRightMargin(0.05)
hx.Draw("col")
hxRMS.SetLineColor(ROOT.kRed)
hxRMS.SetMarkerColor(ROOT.kRed)
hxRMS.SetMarkerStyle(24)
hxRMS.SetMarkerSize(0.6)
hxRMS.Draw("hist ep same")
cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetRightMargin(0.05)
hy.Draw("col")
hyRMS.SetLineColor(ROOT.kRed)
hyRMS.SetMarkerColor(ROOT.kRed)
hyRMS.SetMarkerStyle(24)
hyRMS.SetMarkerSize(0.6)
hyRMS.Draw("hist ep same")
cnv.SaveAs(fn+"res_xy_vs_xy_log_profile.pdf")


peak0 = "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-[2])*(x[0]-[2])+ .25*[1]*[1])"
peak3 = "(0.5*[3]*[4]/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-[5])*(x[0]-[5])+ .25*[4]*[4])"
peak6 = "(0.5*[6]*[7]/TMath::Pi()) / TMath::Max(1.e-10, (x[0]-[8])*(x[0]-[8])+ .25*[7]*[7])"

peakfunc0 = "gaus"
peakfunc6 = "gaus(6)"

cnv = TCanvas("cnv_res_E","",500,500)
hrese = tfile.Get("h_res_E")
g1 = TF1("g1", "gaus", -0.01,+0.01);    g1.SetLineColor(ROOT.kViolet)
g2 = TF1("g2", "gaus", -0.01,+0.01);    g2.SetLineColor(ROOT.kGreen+2)
l0 = TF1("l0", peakfunc0, -0.01,+0.01); l0.SetLineColor(ROOT.kYellow+2)
hrese.Fit(g1,"EMRS")
hrese.Fit(g2,"EMRS")
hrese.Fit(l0,"EMRS")
fite = TF1("fite", "gaus(0)+gaus(3)+"+peakfunc6, -0.01,+0.01)
fite.SetLineWidth(1)
fite.SetParameter(0,g1.GetParameter(0))
fite.SetParameter(1,g1.GetParameter(1))
fite.SetParameter(2,g1.GetParameter(2))
fite.SetParameter(3,g2.GetParameter(0))
fite.SetParameter(4,g2.GetParameter(1))
fite.SetParameter(5,g2.GetParameter(2))
fite.SetParameter(6,l0.GetParameter(0))
fite.SetParameter(7,l0.GetParameter(1))
fite.SetParameter(8,l0.GetParameter(2))
res = hrese.Fit(fite,"EMRS")
chi2dof = fite.GetChisquare()/fite.GetNDF() if(fite.GetNDF()>0) else -1
print("Res(E rec:tru) chi2/Ndof=",chi2dof)
hrese.Draw("hist")
fite.Draw("same")
sigma = GetSigma(fite,"E")
s = ROOT.TLatex()
s.SetNDC(1);
s.SetTextAlign(13);
s.SetTextColor(ROOT.kBlack)
s.SetTextSize(0.03)
s.DrawLatex(0.2,0.85,ROOT.Form("#sigma(E)=%.4f" % (sigma)))
s.DrawLatex(0.2,0.78,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
cnv.SaveAs(fn+"res_E.pdf")

cnv = TCanvas("cnv_res_pi","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
hrespx = tfile.Get("h_res_px")
g1 = TF1("g1", "gaus", -6,+6);    g1.SetLineColor(ROOT.kViolet)
g2 = TF1("g2", "gaus", -6,+6);    g2.SetLineColor(ROOT.kGreen+2)
l0 = TF1("l0", peakfunc0, -6,+6); l0.SetLineColor(ROOT.kYellow+2)
hrespx.Fit(g1,"EMRS")
hrespx.Fit(g2,"EMRS")
hrespx.Fit(l0,"EMRS")
fitpx = TF1("fitpx", "gaus(0)+gaus(3)+"+peakfunc6, -6,+6)
fitpx.SetLineWidth(1)
fitpx.SetParameter(0,g1.GetParameter(0))
fitpx.SetParameter(1,g1.GetParameter(1))
fitpx.SetParameter(2,g1.GetParameter(2))
fitpx.SetParameter(3,g2.GetParameter(0))
fitpx.SetParameter(4,g2.GetParameter(1))
fitpx.SetParameter(5,g2.GetParameter(2))
fitpx.SetParameter(6,l0.GetParameter(0))
fitpx.SetParameter(7,l0.GetParameter(1))
fitpx.SetParameter(8,l0.GetParameter(2))
res = hrespx.Fit(fitpx,"EMRS")
chi2dof = fitpx.GetChisquare()/fitpx.GetNDF() if(fitpx.GetNDF()>0) else -1
print("Res(Px rec:tru) chi2/Ndof=",chi2dof)
hrespx.Draw("hist")
sigma = GetSigma(fitpx,"Px")
s = ROOT.TLatex()
s.SetNDC(1);
s.SetTextAlign(13);
s.SetTextColor(ROOT.kBlack)
s.SetTextSize(0.04)
s.DrawLatex(0.2,0.85,ROOT.Form("#sigma(p_{x})=%.4f" % (sigma)))
s.DrawLatex(0.2,0.78,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
fitpx.Draw("same")

cnv.cd(2)
hrespy = tfile.Get("h_res_py")
g1 = TF1("g1", "gaus", -0.1,+0.1);    g1.SetLineColor(ROOT.kViolet)
g2 = TF1("g2", "gaus", -0.1,+0.1);    g2.SetLineColor(ROOT.kGreen+2)
l0 = TF1("l0", peakfunc0, -0.1,+0.1); l0.SetLineColor(ROOT.kYellow+2)
hrespy.Fit(g1,"EMRS")
hrespy.Fit(g2,"EMRS")
hrespy.Fit(l0,"EMRS")
fitpy = TF1("fitpy", "gaus(0)+gaus(3)+"+peakfunc6, -0.1,+0.1)
fitpy.SetLineWidth(1)
fitpy.SetParameter(0,g1.GetParameter(0))
fitpy.SetParameter(1,g1.GetParameter(1))
fitpy.SetParameter(2,g1.GetParameter(2))
fitpy.SetParameter(3,g2.GetParameter(0))
fitpy.SetParameter(4,g2.GetParameter(1))
fitpy.SetParameter(5,g2.GetParameter(2))
fitpy.SetParameter(6,l0.GetParameter(0))
fitpy.SetParameter(7,l0.GetParameter(1))
fitpy.SetParameter(8,l0.GetParameter(2))
res = hrespy.Fit(fitpy,"EMRS")
chi2dof = fitpy.GetChisquare()/fitpy.GetNDF() if(fitpy.GetNDF()>0) else -1
print("Res(Py rec:tru) chi2/Ndof=",chi2dof)
hrespy.Draw("hist")
sigma = GetSigma(fitpy,"Py")
s = ROOT.TLatex()
s.SetNDC(1);
s.SetTextAlign(13);
s.SetTextColor(ROOT.kBlack)
s.SetTextSize(0.04)
s.DrawLatex(0.2,0.85,ROOT.Form("#sigma(p_{y})=%.4f" % (sigma)))
s.DrawLatex(0.2,0.78,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
fitpy.Draw("same")

cnv.cd(3)
hrespz = tfile.Get("h_res_pz")
g1 = TF1("g1", "gaus", -0.01,+0.01);    g1.SetLineColor(ROOT.kViolet)
g2 = TF1("g2", "gaus", -0.01,+0.01);    g2.SetLineColor(ROOT.kGreen+2)
l0 = TF1("l0", peakfunc0, -0.01,+0.01); l0.SetLineColor(ROOT.kYellow+2)
hrespz.Fit(g1,"EMRS")
hrespz.Fit(g2,"EMRS")
hrespz.Fit(l0,"EMRS")
fitpz = TF1("fitpz", "gaus(0)+gaus(3)+"+peakfunc6, -0.01,+0.01)
fitpz.SetLineWidth(1)
fitpz.SetParameter(0,g1.GetParameter(0))
fitpz.SetParameter(1,g1.GetParameter(1))
fitpz.SetParameter(2,g1.GetParameter(2))
fitpz.SetParameter(3,g2.GetParameter(0))
fitpz.SetParameter(4,g2.GetParameter(1))
fitpz.SetParameter(5,g2.GetParameter(2))
fitpz.SetParameter(6,l0.GetParameter(0))
fitpz.SetParameter(7,l0.GetParameter(1))
fitpz.SetParameter(8,l0.GetParameter(2))
res = hrespz.Fit(fitpz,"EMRS")
chi2dof = fitpz.GetChisquare()/fitpz.GetNDF() if(fitpz.GetNDF()>0) else -1
print("Res(Pz rec:tru) chi2/Ndof=",chi2dof)
hrespz.Draw("hist")
fitpz.Draw("same")
sigma = GetSigma(fitpz,"Pz")
s = ROOT.TLatex()
s.SetNDC(1);
s.SetTextAlign(13);
s.SetTextColor(ROOT.kBlack)
s.SetTextSize(0.04)
s.DrawLatex(0.2,0.85,ROOT.Form("#sigma(p_{z})=%.4f" % (sigma)))
s.DrawLatex(0.2,0.78,ROOT.Form("#chi^{2}/N_{DOF}=%.1f" % (chi2dof)))
cnv.SaveAs(fn+"res_pi.pdf")


# cnv.cd(3)
# hresx = tfile.Get("h_res_x_recvstru")
# g1x = TF1("g1x", "gaus", -0.0004,+0.0004);      g1x.SetLineColor(ROOT.kViolet)
# g2x = TF1("g2x", "gaus", -0.0004,+0.0004);      g2x.SetLineColor(ROOT.kGreen+2)
# l0x = TF1("l0x", peakfunc0, -0.0004,+0.0004); l0x.SetLineColor(ROOT.kYellow+2)
# hresx.Fit(g1x,"EMRS")
# hresx.Fit(g2x,"EMRS")
# hresx.Fit(l0x,"EMRS")
# fitx = TF1("fitx", "gaus(0)+gaus(3)+"+peakfunc6, -0.0004,+0.0004)
# fitx.SetParameter(0,g1x.GetParameter(0))
# fitx.SetParameter(1,g1x.GetParameter(1))
# fitx.SetParameter(2,g1x.GetParameter(2))
# fitx.SetParameter(3,g2x.GetParameter(0))
# fitx.SetParameter(4,g2x.GetParameter(1))
# fitx.SetParameter(5,g2x.GetParameter(2))
# fitx.SetParameter(6,l0x.GetParameter(0))
# fitx.SetParameter(7,l0x.GetParameter(1))
# fitx.SetParameter(8,l0x.GetParameter(2))
# # fitx = TF1("fitx", "gaus", -0.00005,+0.00005)
# res = hresx.Fit(fitx,"EMRS")
# chi2dof = fitx.GetChisquare()/fitx.GetNDF() if(fitx.GetNDF()>0) else -1
# print("Res(x rec:tru) chi2/Ndof=",chi2dof)
# hresx.Draw("hist")
# g1x.Draw("same")
# g2x.Draw("same")
# l0x.Draw("same")
# fitx.Draw("same")
#
# cnv.cd(4)
# hresy = tfile.Get("h_res_y_recvstru")
# g1y = TF1("g1y", "gaus", -0.03,+0.03);          g1y.SetLineColor(ROOT.kViolet)
# g2y = TF1("g2y", "gaus", -0.03,+0.03);          g2y.SetLineColor(ROOT.kGreen+2)
# l0y = TF1("l0y", peakfunc0, -0.0004,+0.0004); l0y.SetLineColor(ROOT.kYellow+2)
# hresy.Fit(g1y,"EMRS")
# hresy.Fit(g2y,"EMRS")
# hresy.Fit(l0y,"EMRS")
# fity = TF1("fity", "gaus(0)+gaus(3)+"+peakfunc6, -0.03,+0.03)
# fity.SetParameter(0,g1y.GetParameter(0))
# fity.SetParameter(1,g1y.GetParameter(1))
# fity.SetParameter(2,g1y.GetParameter(2))
# fity.SetParameter(3,g2y.GetParameter(0))
# fity.SetParameter(4,g2y.GetParameter(1))
# fity.SetParameter(5,g2y.GetParameter(2))
# fity.SetParameter(6,l0y.GetParameter(0))
# fity.SetParameter(7,l0y.GetParameter(1))
# fity.SetParameter(8,l0y.GetParameter(2))
# # fity = TF1("fity", "gaus", -0.01,+0.01)
# res = hresy.Fit(fity,"EMRS")
# chi2dof = fity.GetChisquare()/fity.GetNDF() if(fity.GetNDF()>0) else -1
# print("Res(y rec:tru) chi2/Ndof=",chi2dof)
# hresy.Draw("hist")
# g1y.Draw("same")
# g2y.Draw("same")
# l0y.Draw("same")
# fity.Draw("same")
# cnv.SaveAs(fn+"res_xy_vs_xy_log_profile.pdf")
