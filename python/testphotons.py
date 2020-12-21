#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, TCanvas
import argparse
parser = argparse.ArgumentParser(description='truthtests.py...')
parser.add_argument('-p', metavar='process',   required=True,  help='physics process [trident or bppp]')
parser.add_argument('-d', metavar='directory', required=True,  help='the relative path to the raw root files')
argus = parser.parse_args()
proc = argus.p
path = argus.d

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage = os.path.expandvars("$STORAGEDIR")

meGeV = 0.5109989461/1000.
um2cm = 1.e-4

histos = {}
process = proc
basepath = storage+"/data/root/raw/"+path
targetdir_pdf  = basepath.replace("data/root","output/pdf")
targetdir_root = basepath.replace("data/root","output/root")
p = subprocess.Popen("mkdir -p "+targetdir_pdf, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p "+targetdir_root, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
allpdf  = targetdir_pdf+"/truthphotons_"+process+".pdf"
allroot = targetdir_root+"/truthtest_"+process+".root"



avgE = {"ele":0, "pos":0, "gam":0, "gam05":0, "gam1":0}
ntot = {"ele":0, "pos":0, "gam":0, "gam05":0, "gam1":0}

## book
def Book(process):
   histos.update( {"h_ntot":TH1D("h_ntot","",5,0,5)} )
   histos["h_ntot"].GetXaxis().SetBinLabel(1,"e-")
   histos["h_ntot"].GetXaxis().SetBinLabel(2,"#gamma")
   histos["h_ntot"].GetXaxis().SetBinLabel(3,"#gamma(>0.5 GeV)")
   histos["h_ntot"].GetXaxis().SetBinLabel(4,"#gamma(>1.0 GeV)")
   histos["h_ntot"].GetXaxis().SetBinLabel(5,"e+")

   histos.update( {"h_avgE":TH1D("h_avgE","",5,0,5)} )
   histos["h_avgE"].GetXaxis().SetBinLabel(1,"e-")
   histos["h_avgE"].GetXaxis().SetBinLabel(2,"#gamma")
   histos["h_avgE"].GetXaxis().SetBinLabel(3,"#gamma(>0.5 GeV)")
   histos["h_avgE"].GetXaxis().SetBinLabel(4,"#gamma(>1.0 GeV)")
   histos["h_avgE"].GetXaxis().SetBinLabel(5,"e+")

   Emin = 0
   Emax = 20
   Ebins = 80
   histos.update( {"h_E_photons":TH1D("h_E_photons",";#it{E} [GeV];N_{#gamma}/BX",Ebins,Emin,Emax)} )
   histos.update( {"h_E_electrons":TH1D("h_E_electrons",";#it{E} [GeV];N_{e-}/BX",Ebins,Emin,Emax)} )
   histos.update( {"h_E_positrons":TH1D("h_E_positrons",";#it{E} [GeV];N_{e+}/BX",Ebins,Emin,Emax)} )
   histos.update( {"h_Efine_photons":TH1D("h_Efine_photons",";#it{E} [GeV];N_{#gamma}/BX",4*Ebins,Emin,Emax)} )
   histos.update( {"h_Efine_electrons":TH1D("h_Efine_electrons",";#it{E} [GeV];N_{e-}/BX",4*Ebins,Emin,Emax)} )
   histos.update( {"h_Efine_positrons":TH1D("h_Efine_positrons",";#it{E} [GeV];N_{e+}/BX",4*Ebins,Emin,Emax)} )
   histos.update( {"h_Evfine_photons":TH1D("h_Evfine_photons",";#it{E} [GeV];N_{#gamma}/BX",10*Ebins,Emin,Emax)} )
   histos.update( {"h_Evfine_electrons":TH1D("h_Evfine_electrons",";#it{E} [GeV];N_{e-}/BX",10*Ebins,Emin,Emax)} )
   histos.update( {"h_Evfine_positrons":TH1D("h_Evfine_positrons",";#it{E} [GeV];N_{e+}/BX",10*Ebins,Emin,Emax)} )
 
   ## sumw2
   for hname,hist in histos.items(): hist.Sumw2()


## analysis
def Analyze(n,event):
   nele = 0
   npos = 0
   ngam = 0
   ngam05 = 0
   ngam1 = 0
   ### loop on all tracks
   for j in range(event.E.size()):
      wgt = event.wgt[j]
      p = TLorentzVector()
      p.SetXYZM(event.px[j],event.py[j],event.pz[j],meGeV)
      Etru = event.E[j]
      ### only photons
      if(event.pdgId[j]==22):
         ngam += wgt
         ntot["gam"] += wgt
         avgE["gam"] += wgt*Etru if(not np.isnan(Etru)) else 0
         if(Etru>0.5):
            ngam05 += wgt
            ntot["gam05"] += wgt
            avgE["gam05"] += wgt*Etru if(not np.isnan(Etru)) else 0
         if(Etru>1.0):
            ngam1 += wgt
            ntot["gam1"] += wgt
            avgE["gam1"] += wgt*Etru if(not np.isnan(Etru)) else 0
         histos["h_E_photons"].Fill(Etru,wgt)
         histos["h_Efine_photons"].Fill(Etru,wgt)
         histos["h_Evfine_photons"].Fill(Etru,wgt)
      ### only electrons
      if(event.pdgId[j]==11):
         nele += wgt
         ntot["ele"] += wgt
         avgE["ele"] += wgt*Etru if(not np.isnan(Etru)) else 0
         histos["h_E_electrons"].Fill(Etru,wgt)
         histos["h_Efine_electrons"].Fill(Etru,wgt)
         histos["h_Evfine_electrons"].Fill(Etru,wgt)
      ### only positrons
      if(event.pdgId[j]==-11):
         npos += wgt
         ntot["pos"] += wgt
         avgE["pos"] += wgt*Etru if(not np.isnan(Etru)) else 0
         histos["h_E_positrons"].Fill(Etru,wgt)
         histos["h_Efine_positrons"].Fill(Etru,wgt)
         histos["h_Evfine_positrons"].Fill(Etru,wgt)

   histos["h_ntot"].AddBinContent(1,nele)
   histos["h_ntot"].AddBinContent(2,ngam)
   histos["h_ntot"].AddBinContent(3,ngam05)
   histos["h_ntot"].AddBinContent(4,ngam1)
   histos["h_ntot"].AddBinContent(5,npos)


## event loop
def EventLoop(tree):
   nevents = tree.GetEntries()
   print("with %d events" % nevents)
   n=0 ### init n
   for event in tree:
      if(n%10==0 and n>0):
         print("  processed %d events" % n)
      Analyze(n,event)
      n+=1
   print("Total events processed: ",n)

   ### before normalising
   histos["h_avgE"].SetBinContent(1,avgE["ele"]/ntot["ele"]     if(ntot["ele"]>0)   else 0)
   histos["h_avgE"].SetBinContent(2,avgE["gam"]/ntot["gam"]     if(ntot["gam"]>0)   else 0)
   histos["h_avgE"].SetBinContent(3,avgE["gam05"]/ntot["gam05"] if(ntot["gam05"]>0) else 0)
   histos["h_avgE"].SetBinContent(4,avgE["gam1"]/ntot["gam1"]   if(ntot["gam1"]>0)  else 0)
   histos["h_avgE"].SetBinContent(5,avgE["pos"]/ntot["pos"]     if(ntot["pos"]>0)   else 0)

   ### scale to the BX and to the max bunch charge
   for hname,hist in histos.items():
      if(hname=="h_avgE"): continue
      hist.Scale(4./nevents)

   return nevents

###############################################################

## file and tree
def Run(tfilename,ttreename):
   print("getting tree from ",tfilename)
   tfile = TFile(tfilename,"READ")
   tree = tfile.Get(ttreename)
   nevents = EventLoop(tree)
   return nevents

###############################################################

## actually run
tfilename = basepath+"/raw_"+process+".root"
Book(process)
nevents = Run(tfilename,"tt")
print("nevents=",nevents)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
cnv.SetTicks(1,1)
histos["h_ntot"].SetMinimum(0.9)
histos["h_ntot"].Draw("hist text")
cnv.SaveAs(allpdf+"(")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
cnv.SetTicks(1,1)
histos["h_E_photons"].SetMinimum(0.9)
histos["h_E_photons"].Draw("hist")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
cnv.SetTicks(1,1)
histos["h_E_electrons"].SetMinimum(0.9)
histos["h_E_electrons"].Draw("hist")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
cnv.SetTicks(1,1)
histos["h_E_positrons"].SetMinimum(0.9)
histos["h_E_positrons"].Draw("hist")
cnv.SaveAs(allpdf+")")

print(avgE)
print(ntot)

######################################
tfileout = TFile(allroot,"RECREATE")
tfileout.cd()
for hname,hist in histos.items(): hist.Write()
tfileout.Write()
tfileout.Close()
