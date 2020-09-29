#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, TCanvas, TLegend, TLatex
import argparse
parser = argparse.ArgumentParser(description='truthtests.py...')
parser.add_argument('-p', metavar='process',   required=True,  help='physics process [trident or bppp]')
parser.add_argument('-d', metavar='prefix', required=True,  help='e.g. trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/')
argus = parser.parse_args()
proc = argus.p
prefix = argus.d

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")

###############################################################

signals = {}
signals.update( {"trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/" :{}} )
signals.update( {"bppp/IPstrong_V1.1.00/JETI40/g_laser/16.5GeV/"    :{}} )
signals.update( {"trident/IPstrong_V1.1.00/phaseII/e_laser/16.5GeV/":{}} )
signals.update( {"bppp/IPstrong_V1.1.00/phaseII/g_laser/16.5GeV/"   :{}} )

signals["trident/IPstrong_V1.1.00/JETI40/e_laser/16.5GeV/"] = {
   "w0_3000nm":   {"xi":5.12, "chi":0.900, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_3500nm":   {"xi":4.44, "chi":0.790, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_4000nm":   {"xi":3.88, "chi":0.690, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_4500nm":   {"xi":3.45, "chi":0.610, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_5000nm":   {"xi":3.10, "chi":0.550, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_8000nm":   {"xi":1.94, "chi":0.340, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_20000nm":  {"xi":0.78, "chi":0.138, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_50000nm":  {"xi":0.31, "chi":0.055, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_100000nm": {"xi":0.15, "chi":0.028, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
}
signals["bppp/IPstrong_V1.1.00/JETI40/g_laser/16.5GeV/"] = {
   "w0_3000nm":   {"xi":5.12, "chi":0.900, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_3500nm":   {"xi":4.44, "chi":0.790, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_4000nm":   {"xi":3.88, "chi":0.690, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_4500nm":   {"xi":3.45, "chi":0.610, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_5000nm":   {"xi":3.10, "chi":0.550, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_8000nm":   {"xi":1.94, "chi":0.340, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
}
signals["trident/IPstrong_V1.1.00/phaseII/e_laser/16.5GeV/"] = {
   "w0_8000nm":   {"xi":6.27, "chi":1.11, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_9000nm":   {"xi":5.57, "chi":0.99, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_10000nm":  {"xi":5.01, "chi":0.89, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_11000nm":  {"xi":4.56, "chi":0.81, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_12000nm":  {"xi":4.18, "chi":0.74, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
}
signals["bppp/IPstrong_V1.1.00/phaseII/g_laser/16.5GeV/"] = {
   "w0_8000nm":   {"xi":6.27, "chi":1.11, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_9000nm":   {"xi":5.57, "chi":0.99, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_10000nm":  {"xi":5.01, "chi":0.89, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_11000nm":  {"xi":4.56, "chi":0.81, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
   "w0_12000nm":  {"xi":4.18, "chi":0.74, "Ntru":0, "Nrec":0, "dNtru":0, "dNrec":0},
}

process   = proc
spotsizes = signals[prefix]
histos    = {}
foutname  = prefix.replace(proc+"/IPstrong_V1.1.00/","").replace("/","_")

###############################################################

def minmax(h1,h2,f=1.1):
   hmin = h1.GetMinimum() if(h1.GetMinimum()<h2.GetMinimum()) else h2.GetMinimum()
   hmax = h1.GetMaximum() if(h1.GetMaximum()>h2.GetMaximum()) else h2.GetMaximum()
   h1.SetMaximum(hmax*f)
   h2.SetMaximum(hmax*f)
   return hmin,hmax*f
   
def sumentries(h):
   N = 0
   for b in range(0,h.GetNbinsX()+1): N += h.GetBinContent(b)
   return N
   
def sumerrors(h):
   D = 0
   for b in range(0,h.GetNbinsX()+1): D += h.GetBinError(b)
   return D
   
def getN(h):
   N = 0
   for b in range(1,h.GetNbinsX()+1):
      y = h.GetBinContent(b)
      x = h.GetXaxis().GetBinLowEdge(b)
      N += x*y
   return N
      

###############################################################

for spotsize,properites in spotsizes.items():
   xi  = str(properites["xi"])
   chi = str(properites["chi"])
   tfname = storage+"/output/root/raw/"+prefix+"/"+spotsize+"/truthtest_"+proc+".root"
   print("reading histos from "+storage+"/output/root/raw/"+prefix+"/"+spotsize+"/truthtest_trident.root")
   tf = TFile(tfname,"READ")
   
   spotsizeum = str(int(spotsize.replace("w0_","").replace("nm",""))/1000)+" #mum"
   title = "JETI40, "+spotsizeum+" spot #Rightarrow #xi="+xi+", #chi="+chi
   
   h_tru_E = tf.Get("h_E_positrons").Clone(spotsize+"_tru_E")
   h_tru_E.SetTitle(title)
   h_tru_E.GetXaxis().SetTitle("#it{E}_{e^{+}} [GeV]")
   h_tru_E.GetYaxis().SetTitle("N{e^{+}} per BX per Shot")
   # h_tru_E.Scale(1./h_tru_E.GetXaxis().GetBinWidth(1))
   h_tru_E.SetLineColor(ROOT.kRed)
   h_tru_E.SetDirectory(0)
   h_rec_E = tf.Get("h_E_positrons_rec_emul").Clone(spotsize+"_rec_E")
   h_rec_E.SetTitle(title)
   h_rec_E.GetXaxis().SetTitle("#it{E}_{e^{+}} [GeV]")
   h_rec_E.GetYaxis().SetTitle("N_{e^{+}} per BX per Shot")
   # h_rec_E.Scale(1./h_rec_E.GetXaxis().GetBinWidth(1))
   h_rec_E.SetLineColor(ROOT.kBlack)
   h_rec_E.SetMarkerColor(ROOT.kBlack)
   h_rec_E.SetMarkerStyle(20)
   h_rec_E.SetDirectory(0)
   
   h_tru_ntrks = tf.Get("h_ntrks_positrons").Clone(spotsize+"_tru_ntrks")
   h_tru_ntrks.SetTitle(title)
   h_tru_ntrks.GetXaxis().SetTitle("e^{+} multiplicity")
   h_tru_ntrks.GetYaxis().SetTitle("N_{e^{+}} per BX per Shot")
   h_tru_ntrks.SetLineColor(ROOT.kRed)
   h_tru_ntrks.SetDirectory(0)
   h_rec_ntrks = tf.Get("h_ntrks_positrons_rec_emul").Clone(spotsize+"_rec_ntrks")
   h_rec_ntrks.SetTitle(title)
   h_rec_ntrks.GetXaxis().SetTitle("e^{+} multiplicity")
   h_rec_ntrks.GetYaxis().SetTitle("N_{e^{+}} per BX per Shot")
   h_rec_ntrks.SetLineColor(ROOT.kBlack)
   h_rec_ntrks.SetMarkerColor(ROOT.kBlack)
   h_rec_ntrks.SetMarkerStyle(20)
   h_rec_ntrks.SetDirectory(0)
   
   dNtru = ROOT.Double()
   dNrec = ROOT.Double()
   
   # properites["Ntru"] = h_tru_E.IntegralAndError(1,h_tru_E.GetNbinsX(),dNtru)
   # properites["Nrec"] = h_rec_E.IntegralAndError(1,h_rec_E.GetNbinsX(),dNrec)
   properites["Ntru"] = sumentries(h_tru_E)
   properites["Nrec"] = sumentries(h_rec_E)
   properites["dNtru"] = sumerrors(h_tru_E)
   properites["dNrec"] = sumerrors(h_rec_E)
   # print("NtruInt="+str(properites["Ntru"])+", NtruGet="+str(getN(h_tru_ntrks)))
   # print("NrecInt="+str(properites["Nrec"])+", NrecGet="+str(getN(h_rec_ntrks)))
   
   histos.update({ spotsize+"_tru_E"     : h_tru_E })
   histos.update({ spotsize+"_rec_E"     : h_rec_E })
   histos.update({ spotsize+"_tru_ntrks" : h_tru_ntrks })
   histos.update({ spotsize+"_rec_ntrks" : h_rec_ntrks })

print(spotsizes)
###############################################################

leg_kin = TLegend(0.50,0.70,0.87,0.87)
leg_kin.SetFillStyle(4000) # will be transparent
leg_kin.SetFillColor(0)
leg_kin.SetTextFont(42)
leg_kin.SetBorderSize(0)
leg_kin.AddEntry(histos["w0_8000nm_tru_E"],"Truth signal","l")  ### I assume this spot is always used
leg_kin.AddEntry(histos["w0_8000nm_rec_E"],"Emul. Recon.","lp") ### I assume this spot is always used


###############################################################

xis    = []
chis   = []
Ntrus  = []
Nrecs  = []
dNtrus = []
dNrecs = []
for spotsize,properties in spotsizes.items():
   xis.append( properties["xi"] )
   chis.append( properties["chi"] )
   Ntrus.append( properties["Ntru"] )
   Nrecs.append( properties["Nrec"] )
   dNtrus.append( properties["dNtru"] )
   dNrecs.append( properties["dNrec"] )
h_xi_tru   = TH1D("h_xi_tru", ";Peak #xi;N_{e^{+}}/BX/Shot", len(spotsizes)*1000,0,max(xis)*1.2)
h_xi_rec   = TH1D("h_xi_rec", ";Peak #xi;N_{e^{+}}/BX/Shot", len(spotsizes)*1000,0,max(xis)*1.2)
h_chi_tru  = TH1D("h_chi_tru",";Peak #chi;N_{e^{+}}/BX/Shot",len(spotsizes)*1000,0,max(chis)*1.2)
h_chi_rec  = TH1D("h_chi_rec",";Peak #chi;N_{e^{+}}/BX/Shot",len(spotsizes)*1000,0,max(chis)*1.2)
for i in range(len(spotsizes)):
   bx_xi  = h_xi_tru.FindBin(xis[i])
   bx_chi = h_chi_tru.FindBin(chis[i])
   h_xi_tru.SetBinContent(bx_xi,Ntrus[i])
   h_xi_rec.SetBinContent(bx_xi,Nrecs[i])
   h_chi_tru.SetBinContent(bx_chi,Ntrus[i])
   h_chi_rec.SetBinContent(bx_chi,Nrecs[i])
h_xi_tru.SetMinimum(0)
h_xi_rec.SetMinimum(0)
h_chi_tru.SetMinimum(0)
h_chi_rec.SetMinimum(0)
hmin,hmax = minmax(h_xi_tru,h_xi_rec,1.75)
h_xi_tru.SetBinErrorOption(ROOT.TH1.kPoisson)# if(hmax>1) else ROOT.TH1.kNormal)
h_xi_rec.SetBinErrorOption(ROOT.TH1.kPoisson)# if(hmax>1) else ROOT.TH1.kNormal)
hmin,hmax = minmax(h_chi_tru,h_chi_rec,1.75)
h_chi_tru.SetBinErrorOption(ROOT.TH1.kPoisson)# if(hmax>1) else ROOT.TH1.kNormal)
h_chi_rec.SetBinErrorOption(ROOT.TH1.kPoisson)# if(hmax>1) else ROOT.TH1.kNormal)
for b in range(1,h_xi_tru.GetNbinsX()+1):
   if(h_xi_tru.GetBinContent(b)<=0): h_xi_tru.SetBinContent(b,+99999)
   if(h_xi_rec.GetBinContent(b)<=0): h_xi_rec.SetBinContent(b,+99999)
   if(h_chi_tru.GetBinContent(b)<=0): h_chi_tru.SetBinContent(b,+99999)
   if(h_chi_rec.GetBinContent(b)<=0): h_chi_rec.SetBinContent(b,+99999)
h_xi_tru.SetLineColor(ROOT.kRed)
h_xi_tru.SetMarkerColor(ROOT.kRed)
h_xi_tru.SetMarkerStyle(24)
h_xi_rec.SetLineColor(ROOT.kBlack)
h_xi_rec.SetMarkerColor(ROOT.kBlack)
h_xi_rec.SetMarkerStyle(20)
h_chi_tru.SetLineColor(ROOT.kRed)
h_chi_tru.SetMarkerColor(ROOT.kRed)
h_chi_tru.SetMarkerStyle(24)
h_chi_rec.SetLineColor(ROOT.kBlack)
h_chi_rec.SetMarkerColor(ROOT.kBlack)
h_chi_rec.SetMarkerStyle(20)



###############################################################

##### summarize pdf plots
allpdf  = storage+"/output/pdf/CDRresults_"+foutname+".pdf"
fn = allpdf.replace(".pdf","_")

cnv = TCanvas("cnv","",500,500)
cnv.SaveAs(allpdf+"(")

for spotsize,properties in spotsizes.items():
   cnv = TCanvas("cnv","",500,500)
   # cnv = TCanvas("cnv","",1000,500)
   # cnv.Divide(2,1)
   # cnv.cd(1)
   # ROOT.gPad.SetTicks(1,1)
   # hmin,hmax = minmax(histos[spotsize+"_tru_ntrks"],histos[spotsize+"_rec_ntrks"])
   # histos[spotsize+"_tru_ntrks"].Draw("hist")
   # histos[spotsize+"_rec_ntrks"].Draw("e same")
   # leg_kin.Draw("same")
   
   # cnv.cd(2)
   # ROOT.gPad.SetTicks(1,1)
   cnv.SetTicks(1,1)
   hmin,hmax = minmax(histos[spotsize+"_tru_E"],histos[spotsize+"_rec_E"])
   histos[spotsize+"_tru_E"].Draw("hist")
   histos[spotsize+"_rec_E"].Draw("e same")
   leg_kin.Draw("same")
   
   s = TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(ROOT.kBlack)
   s.SetTextSize(0.035)
   s.DrawLatex(0.47,0.65,ROOT.Form("Tru: #SigmaN_{e^{+}}/BX/Shot=%.2f" % (properties["Ntru"])))
   s.DrawLatex(0.47,0.60,ROOT.Form("Rec: #SigmaN_{e^{+}}/BX/Shot=%.2f" % (properties["Nrec"])))
   
   cnv.SaveAs(fn+"_"+spotsize+"_CDR.pdf")
   cnv.SaveAs(allpdf)


leg_xichi = TLegend(0.15,0.70,0.50,0.87)
leg_xichi.SetFillStyle(4000) # will be transparent
leg_xichi.SetFillColor(0)
leg_xichi.SetTextFont(42)
leg_xichi.SetBorderSize(0)
leg_xichi.AddEntry(h_xi_tru,"Truth signal","lp")
leg_xichi.AddEntry(h_xi_rec,"Emul. Recon.","lp")

s = TLatex()
s.SetNDC(1);
s.SetTextAlign(13);
s.SetTextColor(ROOT.kBlack)
s.SetTextSize(0.035)

cnv = TCanvas("cnv","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetTicks(1,1)
h_xi_tru.Draw("e1p")
h_xi_rec.Draw("e1p same")
leg_xichi.Draw("same")
s.DrawLatex(0.20,0.50,"Stat error only")
s.DrawLatex(0.20,0.45,"set to Poisson")
cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
h_chi_tru.Draw("e1p")
h_chi_rec.Draw("e1p same")
leg_xichi.Draw("same")
s.DrawLatex(0.20,0.50,"Stat error only")
s.DrawLatex(0.20,0.45,"set to Poisson")
cnv.SaveAs(fn+"_xi_chi_CDR.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.SaveAs(allpdf+")")