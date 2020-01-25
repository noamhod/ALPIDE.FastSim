#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-s', metavar='sides',   required=False, help='detector side [e+, e-, e+e-]')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.p
sides = "e+" if(proc=="trident") else "e+e-" ## jsut the default
if(argus.s is not None): sides = argus.s
if(proc=="trident" and "e-" in sides):
   print("ERROR: do not run tracking in the electron side for trident")
   quit()
print("Running with proc=%s and sides=%s" % (proc,sides))

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError

#############################################

### electron mass:
me = 0.51099895/1000. ### GeV
me2 = me*me
cm2m = 1.e-2
cm2um = 1.e4
um2cm = 1.e-4

### magnetic field
B  = 1.4 if(proc=="trident") else 2.0 # Tesla
LB = 1   # meters

### possible energies 
Emax = 17.5 if(proc=="trident") else 16 # GeV
Emin = 1.00 if(proc=="trident") else 2 # GeV

### geometry:
zDipoleExit = 202.9
xDipoleExitMinAbs = 1.5 if(proc=="bppp") else 4   ## cm --> TODO: need tuning
xDipoleExitMaxAbs = 25  if(proc=="bppp") else 30  ## cm --> TODO: need tuning
xDipoleExitAbs = (xDipoleExitMaxAbs-xDipoleExitMinAbs)/2.
yDipoleExitMin = -0.05 ## cm --> TODO: need tuning
yDipoleExitMax = +0.05 ## cm --> TODO: need tuning
xAbsMargins = 0.025 # cm --> TODO: need tuning
yAbsMargins = 0.025 if(proc=="bppp") else 0.1 # cm --> TODO: need tuning

### stave geometry
Hstave    = 1.5  # cm
Lstave    = 27 if(proc=="bppp") else 50 # cm
Rbeampipe = 4 # cm
RoffsetBfield22BPPP = 7.0  # cm for BPPP in B=2.2T
RoffsetBfield20BPPP = 5.7  # cm for BPPP in B=2.0T
RoffsetBfield14BPPP = 4.0  # cm for BPPP in B=1.4T
RoffsetBfield = RoffsetBfield20BPPP if(proc=="bppp") else 14 # cm
xPsideL = -RoffsetBfield-Lstave
xPsideR = -RoffsetBfield       
xEsideL = +RoffsetBfield       
xEsideR = +RoffsetBfield+Lstave
yUp = +Hstave/2.
yDn = -Hstave/2.

### for the histogram
detXmin = xPsideL
detXmax = xEsideR
if(proc=="trident"): detXmax = xPsideR
if(proc=="bppp" and sides=="e+"): detXmax = xPsideR
if(proc=="bppp" and sides=="e-"): detXmin = xEsideL

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

nbins_E, bins_E = GetLogBinning(40,1,15) if(proc=="bppp") else GetLogBinning(20,1,4.5)
#############################################

histos = { "h_residuals_xz_sig": TH1D("residuals_xz_sig",";residuals_{xz};Tracks", 500,0,0.5),
           "h_residuals_yz_sig": TH1D("residuals_yz_sig",";residuals_{yz};Tracks", 500,0,500), 
           "h_residuals_xz_bkg": TH1D("residuals_xz_bkg",";residuals_{xz};Tracks", 500,0,0.5),
           "h_residuals_yz_bkg": TH1D("residuals_yz_bkg",";residuals_{yz};Tracks", 500,0,500),
           
           "h_svd_dd0_sig": TH1D("svd_dd0_sig",";svd_{dd0};Tracks", 500,21,24),
           "h_svd_dd0_bkg": TH1D("svd_dd0_bkg",";svd_{dd0};Tracks", 500,21,24), 
           "h_svd_dd1_sig": TH1D("svd_dd1_sig",";svd_{dd1};Tracks", 500,0,0.1),
           "h_svd_dd1_bkg": TH1D("svd_dd1_bkg",";svd_{dd1};Tracks", 500,0,0.1),
           "h_svd_dd2_sig": TH1D("svd_dd2_sig",";svd_{dd2};Tracks", 500,0,0.05),
           "h_svd_dd2_bkg": TH1D("svd_dd2_bkg",";svd_{dd2};Tracks", 500,0,0.05),
           
           "h_prob_xz_sig": TH1D("prob_xz_sig",";prob_{xz};Tracks", 500,0,1.0),
           "h_prob_yz_sig": TH1D("prob_yz_sig",";prob_{yz};Tracks", 500,0,1.0), 
           "h_prob_xz_bkg": TH1D("prob_xz_bkg",";prob_{xz};Tracks", 500,0,1.0),
           "h_prob_yz_bkg": TH1D("prob_yz_bkg",";prob_{yz};Tracks", 500,0,1.0),
           
           "h_chi2ndf_xz_sig": TH1D("chi2ndf_xz_sig",";chi2ndf_{xz};Tracks", 500,0,0.001),
           "h_chi2ndf_yz_sig": TH1D("chi2ndf_yz_sig",";chi2ndf_{yz};Tracks", 500,0,0.001), 
           "h_chi2ndf_xz_bkg": TH1D("chi2ndf_xz_bkg",";chi2ndf_{xz};Tracks", 500,0,0.001),
           "h_chi2ndf_yz_bkg": TH1D("chi2ndf_yz_bkg",";chi2ndf_{yz};Tracks", 500,0,0.001),
           
           "h_seed_resE" : TH1D("seed_resE", ";(E_{seed}-E_{gen})/E_{gen};Tracks",    100,-0.05,+0.05),
           "h_seed_resPz": TH1D("seed_resPz",";(Pz_{seed}-Pz_{gen})/Pz_{gen};Tracks", 100,-0.05,+0.05), 
           "h_seed_resPy": TH1D("seed_resPy",";(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,-10,+10),
           
           "h_seed_resE_vs_x"  : TH2D("seed_resE_vs_x",  ";x;(E_{seed}-E_{gen})/E_{gen};Tracks",    100,detXmin,detXmax, 100,-0.05,+0.05),
           "h_seed_resPy_vs_x" : TH2D("seed_resPy_vs_x", ";x;(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           
           "h_N_sigacc":        TH1D("N_sigacc",        ";Track multiplicity;Events", 100,30,330),
           "h_N_all_seeds":     TH1D("N_all_seeds",     ";Track multiplicity;Events", 100,30,330),
           "h_N_matched_seeds": TH1D("N_matched_seeds", ";Track multiplicity;Events", 100,30,330),
           "h_N_good_seeds":    TH1D("N_good_seeds",    ";Track multiplicity;Events", 100,30,330),
           
           "h_seeding_score": TH1D("h_seeding_score", ";N_{seeds}^{matched}/N_{signa}^{in.acc} [%];Events", 20,91,101),
           "h_seeding_pool":  TH1D("h_seeding_pool",  ";N_{seeds}^{all}/N_{signa}^{in.acc} [%];Events", 50,90,590),
           
           "h_gen_E_all"     : TH1D("gen_E_all",     ";E_{gen} [GeV];Tracks",  nbins_E, bins_E),
           "h_gen_E_mat_all" : TH1D("gen_E_mat_all", ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_gen_E_mat_sel" : TH1D("gen_E_mat_sel", ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_eff_E_mat_all" : TH1D("eff_E_mat_all", ";E_{gen} [GeV];Efficiency", nbins_E, bins_E),
           "h_eff_E_mat_sel" : TH1D("eff_E_mat_sel", ";E_{gen} [GeV];Efficiency", nbins_E, bins_E),
}
histos["h_gen_E_all"].Sumw2()
histos["h_gen_E_mat_all"].Sumw2()
histos["h_gen_E_mat_sel"].Sumw2()
histos["h_eff_E_mat_all"].Sumw2()
histos["h_eff_E_mat_sel"].Sumw2()


frootin = "../data/root/seeds_"+proc+".root"
print("Opening root file:",frootin)
intfile = TFile(frootin,"READ")
intree = intfile.Get("seeds")
nevents = intree.GetEntries()
print("with %d events" % nevents)
n=0 ### init n
for event in intree:
   # Nsigall = 
   Nsigacc = event.iGen.size()
   Nseeds = event.eSeed.size()
   Nmatched = 0
   Ngood = 0
   
   for j in range(event.eGen.size()):
      histos["h_gen_E_all"].Fill(event.eGen[j])
   
   for i in range(event.eSeed.size()):
      # event.svd0Seed
      # event.svd1Seed
      # event.svd2Seed
      # event.chi2xzSeed
      # event.chi2yzSeed
      # event.residxzSeed
      # event.residyzSeed
      # event.issigSeed
      # event.iGenMatch
      #
      # event.x1Seed
      # event.y1Seed
      # event.z1Seed
      # event.x2Seed
      # event.y2Seed
      # event.z2Seed
      # event.x3Seed
      # event.y3Seed
      # event.z3Seed
      # event.x4Seed
      # event.y4Seed
      # event.z4Seed
      #
      # event.pxSeed
      # event.pySeed
      # event.pzSeed
      # event.eSeed
      #
      # event.pxGen
      # event.pyGen
      # event.pzGen
      # event.eGen
      # event.qGen
      # event.iGen
      
      ### get the matching seed if matched
      igen = event.iGenMatch[i]
      foundmatch = (igen>=0 and event.issigSeed[i])
      pgen = TLorentzVector()
      if(foundmatch): pgen.SetPxPyPzE(event.pxGen[igen],event.pyGen[igen],event.pzGen[igen],event.eGen[igen])
      
      if(event.issigSeed[i]): 
         Nmatched += 1
      
      ### 2d fit
      if(event.issigSeed[i]):
         histos["h_chi2ndf_xz_sig"].Fill(event.chi2xzSeed[i])
         histos["h_chi2ndf_yz_sig"].Fill(event.chi2yzSeed[i])
      else:
         histos["h_chi2ndf_xz_bkg"].Fill(event.chi2xzSeed[i])
         histos["h_chi2ndf_yz_bkg"].Fill(event.chi2yzSeed[i])

      ### a single 3d fit
      if(event.issigSeed[i]):
         histos["h_residuals_xz_sig"].Fill(event.residxzSeed[i])
         histos["h_residuals_yz_sig"].Fill(event.residyzSeed[i])
      else:
         histos["h_residuals_xz_bkg"].Fill(event.residxzSeed[i])
         histos["h_residuals_yz_bkg"].Fill(event.residyzSeed[i])
         
      ### the SVD alg
      if(event.issigSeed[i]):
         histos["h_svd_dd0_sig"].Fill(event.svd0Seed[i])
         histos["h_svd_dd1_sig"].Fill(event.svd1Seed[i])
         histos["h_svd_dd2_sig"].Fill(event.svd2Seed[i])
      else:
         histos["h_svd_dd0_bkg"].Fill(event.svd0Seed[i])
         histos["h_svd_dd1_bkg"].Fill(event.svd1Seed[i])
         histos["h_svd_dd2_bkg"].Fill(event.svd2Seed[i])
      
      ### for efficiencies before selection  
      if(foundmatch): histos["h_gen_E_mat_all"].Fill(event.eGen[igen])
      
      ### cut on some quality
      isgood = (event.svd1Seed[i]<0.005 and event.svd2Seed[i]<0.0025)
      if(not isgood): continue
      Ngood += 1
      
      ### get the seed
      pseed = TLorentzVector()
      pseed.SetPxPyPzE(event.pxSeed[i],event.pySeed[i],event.pzSeed[i],event.eSeed[i])
      
      ### get the matching seed
      if(foundmatch):
         ### for efficiencies after selection
         histos["h_gen_E_mat_sel"].Fill(event.eGen[igen])
         ### check momentum perforrmance of seeding
         resE = (pseed.E()-pgen.E())/pgen.E()
         resPz = (pseed.Pz()-pgen.Pz())/pgen.Pz()
         resPy = (pseed.Py()-pgen.Py())/pgen.Py()
         histos["h_seed_resE"].Fill(resE)
         histos["h_seed_resPz"].Fill(resPz)
         histos["h_seed_resPy"].Fill(resPy)
         ### check spatial perforrmance of seeding    
         histos["h_seed_resE_vs_x"].Fill(event.x4Seed[i],resE)
         histos["h_seed_resPy_vs_x"].Fill(event.x4Seed[i],resPy)
   
   ### seeding stat       
   histos["h_N_sigacc"].Fill(Nsigacc)       
   histos["h_N_all_seeds"].Fill(Nseeds)       
   histos["h_N_matched_seeds"].Fill(Nmatched)       
   histos["h_N_good_seeds"].Fill(Ngood)
   histos["h_seeding_score"].Fill(Nmatched/Nsigacc*100)
   histos["h_seeding_pool"].Fill(Nseeds/Nsigacc*100)
   
   # print("Event: %g --> Nsigacc=%g, Nseeds=%g, Nmatched=%g, Ngood=%g --> Seeds matching performance: Nmatched/Nsigacc=%5.1f%%" % (n,Nsigacc,Nseeds,Nmatched,Ngood,Nmatched/Nsigacc*100))
   if(n%100==0 and n>0): print("  processed %d events" % n)
   n+=1
print("Total events processed: ",n-1)





fpdf = "../output/pdf/seedingdemoplots_"+proc+".pdf"

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetLogy()
histos["h_chi2ndf_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_xz_sig"].Draw()
histos["h_chi2ndf_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_xz_bkg"].Draw("same")
cnv.cd(2)
ROOT.gPad.SetLogy()
histos["h_chi2ndf_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_yz_sig"].Draw()
histos["h_chi2ndf_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_yz_bkg"].Draw("same")
cnv.SaveAs(fpdf+"(")

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetLogy()
histos["h_residuals_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_xz_sig"].Draw()
histos["h_residuals_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_xz_bkg"].Draw("same")
cnv.cd(2)
ROOT.gPad.SetLogy()
histos["h_residuals_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_yz_sig"].Draw()
histos["h_residuals_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_yz_bkg"].Draw("same")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
ROOT.gPad.SetLogy()
histos["h_svd_dd0_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd0_sig"].Draw()
histos["h_svd_dd0_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd0_bkg"].Draw("same")
cnv.cd(2)
ROOT.gPad.SetLogy()
histos["h_svd_dd1_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd1_sig"].Draw()
histos["h_svd_dd1_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd1_bkg"].Draw("same")
cnv.cd(3)
ROOT.gPad.SetLogy()
histos["h_svd_dd2_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd2_sig"].Draw()
histos["h_svd_dd2_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd2_bkg"].Draw("same")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1000,1000)
cnv.Divide(2,2)
cnv.cd(1); histos["h_seed_resE"].Draw("hist")
cnv.cd(2); histos["h_seed_resPy"].Draw("hist")
cnv.cd(3); histos["h_seed_resE_vs_x"].Draw("col")
cnv.cd(4); histos["h_seed_resPy_vs_x"].Draw("col")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
histos["h_N_sigacc"].SetLineColor(ROOT.kBlack); histos["h_N_sigacc"].Draw()
histos["h_N_all_seeds"].SetLineColor(ROOT.kBlue); histos["h_N_all_seeds"].Draw("same")
histos["h_N_matched_seeds"].SetLineColor(ROOT.kGreen); histos["h_N_matched_seeds"].Draw("same")
histos["h_N_good_seeds"].SetLineColor(ROOT.kRed); histos["h_N_good_seeds"].Draw("same")
cnv.cd(2)
histos["h_seeding_pool"].Draw("hist")
cnv.cd(3)
histos["h_seeding_score"].Draw("hist")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_gen_E_all"].SetLineColor(ROOT.kBlack);     histos["h_gen_E_all"].Draw("hist")
histos["h_gen_E_mat_all"].SetLineColor(ROOT.kGreen); histos["h_gen_E_mat_all"].Draw("hist same")
histos["h_gen_E_mat_sel"].SetLineColor(ROOT.kRed);   histos["h_gen_E_mat_sel"].Draw("hist same")
cnv.cd(2)
histos["h_eff_E_mat_all"].Divide(histos["h_gen_E_mat_all"],histos["h_gen_E_all"])
histos["h_eff_E_mat_sel"].Divide(histos["h_gen_E_mat_sel"],histos["h_gen_E_all"])
for b in range(histos["h_eff_E_mat_all"].GetNbinsX()+1):
   if(histos["h_eff_E_mat_all"].GetBinContent(b)>1): histos["h_eff_E_mat_all"].SetBinContent(b,1)
   if(histos["h_eff_E_mat_sel"].GetBinContent(b)>1): histos["h_eff_E_mat_sel"].SetBinContent(b,1)
histos["h_eff_E_mat_all"].SetLineColor(ROOT.kGreen); histos["h_eff_E_mat_all"].SetMinimum(0.0); histos["h_eff_E_mat_all"].SetMaximum(1.1)
histos["h_eff_E_mat_sel"].SetLineColor(ROOT.kRed);   histos["h_eff_E_mat_sel"].SetMinimum(0.0); histos["h_eff_E_mat_sel"].SetMaximum(1.1)
histos["h_eff_E_mat_all"].Draw("ep")
histos["h_eff_E_mat_sel"].Draw("ep same")
cnv.SaveAs(fpdf+")")





