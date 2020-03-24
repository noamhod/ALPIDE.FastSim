#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse

parser = argparse.ArgumentParser(description='RecoDemoPlot.py...')
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
           
           "h_seed_vs_gen_resE" : TH1D("seed_vs_gen_resE", ";(E_{seed}-E_{gen})/E_{gen};Tracks",    100,-0.05,+0.05),
           "h_seed_vs_gen_resPz": TH1D("seed_vs_gen_resPz",";(Pz_{seed}-Pz_{gen})/Pz_{gen};Tracks", 100,-0.05,+0.05), 
           "h_seed_vs_gen_resPy": TH1D("seed_vs_gen_resPy",";(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,-10,+10),
           "h_seed_vs_gen_resPx": TH1D("seed_vs_gen_resPx",";(Px_{seed}-Px_{gen})/Px_{gen};Tracks", 100,-10,+10),
           
           "h_rec_vs_gen_resE" : TH1D("rec_vs_gen_resE", ";(E_{rec}-E_{gen})/E_{gen};Tracks",    100,-0.05,+0.05),
           "h_rec_vs_gen_resPz": TH1D("rec_vs_gen_resPz",";(Pz_{rec}-Pz_{gen})/Pz_{gen};Tracks", 100,-0.05,+0.05), 
           "h_rec_vs_gen_resPy": TH1D("rec_vs_gen_resPy",";(Py_{rec}-Py_{gen})/Py_{gen};Tracks", 100,-10,+10),
           "h_rec_vs_gen_resPx": TH1D("rec_vs_gen_resPx",";(Px_{rec}-Px_{gen})/Px_{gen};Tracks", 100,-10,+10),
           
           "h_rec_vs_seed_resE" : TH1D("rec_vs_seed_resE", ";(E_{rec}-E_{seed})/E_{seed};Tracks",    100,-0.05,+0.05),
           "h_rec_vs_seed_resPz": TH1D("rec_vs_seed_resPz",";(Pz_{rec}-Pz_{seed})/Pz_{seed};Tracks", 100,-0.05,+0.05), 
           "h_rec_vs_seed_resPy": TH1D("rec_vs_seed_resPy",";(Py_{rec}-Py_{seed})/Py_{seed};Tracks", 100,-10,+10),
           "h_rec_vs_seed_resPx": TH1D("rec_vs_seed_resPx",";(Px_{rec}-Px_{seed})/Px_{seed};Tracks", 100,-10,+10),
           
           "h_seed_vs_gen_resE_vs_x"  : TH2D("seed_vs_gen_resE_vs_x",  ";x;(E_{seed}-E_{gen})/E_{gen};Tracks",    100,detXmin,detXmax, 100,-0.05,+0.05),
           "h_seed_vs_gen_resPy_vs_x" : TH2D("seed_vs_gen_resPy_vs_x", ";x;(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           "h_seed_vs_gen_resPx_vs_x" : TH2D("seed_vs_gen_resPx_vs_x", ";x;(Px_{seed}-Px_{gen})/Px_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),

           "h_rec_vs_gen_resE_vs_x"  : TH2D("rec_vs_gen_resE_vs_x",  ";x;(E_{rec}-E_{gen})/E_{gen};Tracks",    100,detXmin,detXmax, 100,-0.05,+0.05),
           "h_rec_vs_gen_resPy_vs_x" : TH2D("rec_vs_gen_resPy_vs_x", ";x;(Py_{rec}-Py_{gen})/Py_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           "h_rec_vs_gen_resPx_vs_x" : TH2D("rec_vs_gen_resPx_vs_x", ";x;(Px_{rec}-Px_{gen})/Px_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           
           "h_rec_vs_seed_resE_vs_x"  : TH2D("rec_vs_seed_resE_vs_x",  ";x;(E_{rec}-E_{seed})/E_{seed};Tracks",    100,detXmin,detXmax, 100,-0.05,+0.05),
           "h_rec_vs_seed_resPy_vs_x" : TH2D("rec_vs_seed_resPy_vs_x", ";x;(Py_{rec}-Py_{seed})/Py_{seed};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           "h_rec_vs_seed_resPx_vs_x" : TH2D("rec_vs_seed_resPx_vs_x", ";x;(Px_{rec}-Px_{seed})/Px_{seed};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           
           "h_N_sigacc":        TH1D("N_sigacc",        ";Track multiplicity;Events", 100,30,330),
           "h_N_all_seeds":     TH1D("N_all_seeds",     ";Track multiplicity;Events", 100,30,330),
           "h_N_matched_seeds": TH1D("N_matched_seeds", ";Track multiplicity;Events", 100,30,330),
           "h_N_good_seeds":    TH1D("N_good_seeds",    ";Track multiplicity;Events", 100,30,330),
           "h_N_all_rec":       TH1D("N_all_rec",     ";Track multiplicity;Events", 100,30,330),
           "h_N_matched_rec":   TH1D("N_matched_rec", ";Track multiplicity;Events", 100,30,330),
           "h_N_good_rec":      TH1D("N_good_rec",      ";Track multiplicity;Events", 100,30,330),
           
           "h_seeding_score": TH1D("h_seeding_score", ";N_{seeds}^{matched}/N_{signa}^{in.acc} [%];Events", 20,91,101),
           "h_seeding_pool":  TH1D("h_seeding_pool",  ";N_{seeds}^{all}/N_{signa}^{in.acc} [%];Events", 50,90,590),
           
           "h_gen_E_all"                 : TH1D("gen_E_all",                 ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_gen_E_seedMatched2gen_all" : TH1D("gen_E_seedMatched2gen_all", ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_gen_E_seedMatched2gen_sel" : TH1D("gen_E_seedMatched2gen_sel", ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_gen_E_recMatched2gen_all"  : TH1D("gen_E_recMatched2gen_all",  ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           "h_gen_E_recMatched2gen_sel"  : TH1D("gen_E_recMatched2gen_sel",  ";E_{gen} [GeV];Tracks", nbins_E, bins_E),
           
           "h_eff_E_seedMatched2gen_all" : TH1D("eff_E_seedMatched2gen_all", "Seeding efficiency;E_{gen} [GeV];Efficiency", nbins_E, bins_E),
           "h_eff_E_seedMatched2gen_sel" : TH1D("eff_E_seedMatched2gen_sel", "Seeding efficiency;E_{gen} [GeV];Efficiency", nbins_E, bins_E),
           "h_eff_E_recMatched2gen_all"  : TH1D("eff_E_recMatched2gen_all",  "Reconstruction efficiency;E_{gen} [GeV];Efficiency", nbins_E, bins_E),
           "h_eff_E_recMatched2gen_sel"  : TH1D("eff_E_recMatched2gen_sel",  "Reconstruction efficiency;E_{gen} [GeV];Efficiency", nbins_E, bins_E),
           
}
histos["h_gen_E_all"].Sumw2()
histos["h_gen_E_seedMatched2gen_all"].Sumw2()
histos["h_gen_E_seedMatched2gen_sel"].Sumw2()
histos["h_gen_E_recMatched2gen_all"].Sumw2()
histos["h_gen_E_recMatched2gen_sel"].Sumw2()
histos["h_eff_E_seedMatched2gen_all"].Sumw2()
histos["h_eff_E_seedMatched2gen_sel"].Sumw2()
histos["h_eff_E_recMatched2gen_all"].Sumw2()
histos["h_eff_E_recMatched2gen_sel"].Sumw2()


frootin = "../data/root/rec_from_seeds_"+proc+".root"
print("Opening root file:",frootin)
intfile = TFile(frootin,"READ")
intree = intfile.Get("res")

# event.n_gen
# event.j_gen
# event.q_gen
# event.p_gen
#
# event.n_seed
# event.jgen_seed
# event.q_seed
# event.p_seed
# event.acc_seed
# event.polm_seed
# event.poll_seed
# event.svd0_seed
# event.svd1_seed
# event.svd2_seed
# event.chi2xz_seed
# event.chi2yz_seed
# event.residxz_seed
# event.residyz_seed
# event.issig_seed
# event.x1_seed
# event.y1_seed
# event.z1_seed
# event.x2_seed
# event.y2_seed
# event.z2_seed
# event.x3_seed
# event.y3_seed
# event.z3_seed
# event.x4_seed
# event.y4_seed
# event.z4_seed
#
# event.n_rec
# event.jgen_rec
# event.jseed_rec
# event.q_rec
# event.p_rec
# event.acc_rec
# event.polm_rec
# event.poll_rec
# event.issig_rec

nevents = intree.GetEntries()
print("with %d events" % nevents)
n=0 ### init n
for event in intree:
   # Nsigall = 
   Nsigacc = event.n_gen
   Nseeds  = event.n_seed
   Nrec    = event.n_rec
   N_seed_matched = 0
   N_seed_good    = 0
   N_rec_matched  = 0
   N_rec_good     = 0
   
   for j in range(Nsigacc):
      histos["h_gen_E_all"].Fill(event.p_gen[j].E())
   
   for i in range(Nseeds):
      ### get the matching seed if matched
      igen = event.jgen_seed[i]
      foundmatch = (igen>=0 and event.issig_seed[i])
      
      if(event.issig_seed[i]): N_seed_matched += 1
      
      ### 2d fit
      if(event.issig_seed[i]):
         histos["h_chi2ndf_xz_sig"].Fill(event.chi2xz_seed[i])
         histos["h_chi2ndf_yz_sig"].Fill(event.chi2yz_seed[i])
      else:
         histos["h_chi2ndf_xz_bkg"].Fill(event.chi2xz_seed[i])
         histos["h_chi2ndf_yz_bkg"].Fill(event.chi2yz_seed[i])

      ### a single 3d fit
      if(event.issig_seed[i]):
         histos["h_residuals_xz_sig"].Fill(event.residxz_seed[i])
         histos["h_residuals_yz_sig"].Fill(event.residyz_seed[i])
      else:
         histos["h_residuals_xz_bkg"].Fill(event.residxz_seed[i])
         histos["h_residuals_yz_bkg"].Fill(event.residyz_seed[i])
         
      ### the SVD alg
      if(event.issig_seed[i]):
         histos["h_svd_dd0_sig"].Fill(event.svd0_seed[i])
         histos["h_svd_dd1_sig"].Fill(event.svd1_seed[i])
         histos["h_svd_dd2_sig"].Fill(event.svd2_seed[i])
      else:
         histos["h_svd_dd0_bkg"].Fill(event.svd0_seed[i])
         histos["h_svd_dd1_bkg"].Fill(event.svd1_seed[i])
         histos["h_svd_dd2_bkg"].Fill(event.svd2_seed[i])
      
      ### for efficiencies before selection  
      if(foundmatch): histos["h_gen_E_seedMatched2gen_all"].Fill(event.p_gen[igen].E())
      
      ### cut on some quality
      isgood = (event.svd1_seed[i]<0.005 and event.svd2_seed[i]<0.0025)
      if(not isgood): continue
      N_seed_good += 1
      
      ### get the matching seed
      if(foundmatch):
         ### for efficiencies after selection
         histos["h_gen_E_seedMatched2gen_sel"].Fill(event.p_gen[igen].E())
         ### check momentum perforrmance of seeding
         resE  = (event.p_seed[i].E() -event.p_gen[igen].E()) /event.p_gen[igen].E()
         resPz = (event.p_seed[i].Pz()-event.p_gen[igen].Pz())/event.p_gen[igen].Pz()
         resPy = (event.p_seed[i].Py()-event.p_gen[igen].Py())/event.p_gen[igen].Py()
         resPx = (event.p_seed[i].Px()-event.p_gen[igen].Px())/event.p_gen[igen].Px()
         histos["h_seed_vs_gen_resE"].Fill(resE)
         histos["h_seed_vs_gen_resPz"].Fill(resPz)
         histos["h_seed_vs_gen_resPy"].Fill(resPy)
         histos["h_seed_vs_gen_resPx"].Fill(resPx)
         ### check spatial perforrmance of seeding    
         histos["h_seed_vs_gen_resE_vs_x"].Fill(event.x4_seed[i],resE)
         histos["h_seed_vs_gen_resPy_vs_x"].Fill(event.x4_seed[i],resPy)
         histos["h_seed_vs_gen_resPx_vs_x"].Fill(event.x4_seed[i],resPx)
   
   
   for i in range(Nrec):
      ### get the matching seed if matched
      igen  = event.jgen_rec[i]
      iseed = event.jseed_rec[i]
      foundmatch_gen  = (igen>=0  and event.issig_rec[i])
      foundmatch_seed = (iseed>=0 and event.issig_rec[i])
      
      if(event.issig_rec[i]): N_rec_matched += 1
      
      ### for efficiencies before selection  
      if(foundmatch_gen): histos["h_gen_E_recMatched2gen_all"].Fill(event.p_gen[igen].E())
      
      ### cut on some quality
      isgood = (event.svd1_seed[iseed]<0.005 and event.svd2_seed[iseed]<0.0025)
      if(not isgood): continue
      N_rec_good += 1
      
      ### get the matching seed
      if(foundmatch_gen):
         ### for efficiencies after selection
         histos["h_gen_E_recMatched2gen_sel"].Fill(event.p_gen[igen].E())
         
         ### check momentum perforrmance of reconstruction
         resE  = (event.p_rec[i].E() -event.p_gen[igen].E()) /event.p_gen[igen].E()
         resPz = (event.p_rec[i].Pz()-event.p_gen[igen].Pz())/event.p_gen[igen].Pz()
         resPy = (event.p_rec[i].Py()-event.p_gen[igen].Py())/event.p_gen[igen].Py()
         resPx = (event.p_rec[i].Px()-event.p_gen[igen].Px())/event.p_gen[igen].Px()
         histos["h_rec_vs_gen_resE"].Fill(resE)
         histos["h_rec_vs_gen_resPz"].Fill(resPz)
         histos["h_rec_vs_gen_resPy"].Fill(resPy)
         histos["h_rec_vs_gen_resPx"].Fill(resPx)
         ### check spatial perforrmance of reconstruction    
         histos["h_rec_vs_gen_resE_vs_x"].Fill(event.x4_seed[iseed],resE)
         histos["h_rec_vs_gen_resPy_vs_x"].Fill(event.x4_seed[iseed],resPy)
         histos["h_rec_vs_gen_resPx_vs_x"].Fill(event.x4_seed[iseed],resPx)
         
         resE  = (event.p_rec[i].E() -event.p_seed[iseed].E()) /event.p_seed[iseed].E()
         resPz = (event.p_rec[i].Pz()-event.p_seed[iseed].Pz())/event.p_seed[iseed].Pz()
         resPy = (event.p_rec[i].Py()-event.p_seed[iseed].Py())/event.p_seed[iseed].Py()
         resPx = (event.p_rec[i].Px()-event.p_seed[iseed].Px())/event.p_seed[iseed].Px()
         histos["h_rec_vs_seed_resE"].Fill(resE)
         histos["h_rec_vs_seed_resPz"].Fill(resPz)
         histos["h_rec_vs_seed_resPy"].Fill(resPy)
         histos["h_rec_vs_seed_resPx"].Fill(resPx)
         ### check spatial perforrmance of reconstruction
         histos["h_rec_vs_seed_resE_vs_x"].Fill(event.x4_seed[iseed],resE)
         histos["h_rec_vs_seed_resPy_vs_x"].Fill(event.x4_seed[iseed],resPy)
         histos["h_rec_vs_seed_resPx_vs_x"].Fill(event.x4_seed[iseed],resPx)
   
   
   ### seeding stat       
   histos["h_N_sigacc"].Fill(Nsigacc)       
   histos["h_N_all_seeds"].Fill(Nseeds)       
   histos["h_N_matched_seeds"].Fill(N_seed_matched)       
   histos["h_N_good_seeds"].Fill(N_seed_good)
   histos["h_N_all_rec"].Fill(Nrec)
   histos["h_N_matched_rec"].Fill(N_rec_matched)       
   histos["h_N_good_rec"].Fill(N_rec_good)
   histos["h_seeding_score"].Fill(N_seed_matched/Nsigacc*100)
   histos["h_seeding_pool"].Fill(Nseeds/Nsigacc*100)
   
   # print("Event: %g --> Ngen=%g, Nseed=%g, N_seed_mat=%g, N_seed_good=%g, N_rec_good=%g --> Matching: N_seed_mat/Ngenacc=%5.1f%% " % (n,Nsigacc,Nseeds,N_seed_matched,N_seed_good,N_rec_good,N_seed_matched/Nsigacc*100))
   if(n%100==0 and n>0): print("  processed %d events" % n)
   n+=1
print("Total events processed: ",n-1)





fpdf = "../output/pdf/recodemoplots_"+proc+".pdf"

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

cnv = TCanvas("","",1500,1000)
cnv.Divide(3,2)
cnv.cd(1); histos["h_seed_vs_gen_resE"].Draw("hist")
cnv.cd(2); histos["h_seed_vs_gen_resPy"].Draw("hist")
cnv.cd(3); histos["h_seed_vs_gen_resPx"].Draw("hist")
cnv.cd(4); histos["h_seed_vs_gen_resE_vs_x"].Draw("col")
cnv.cd(5); histos["h_seed_vs_gen_resPy_vs_x"].Draw("col")
cnv.cd(6); histos["h_seed_vs_gen_resPx_vs_x"].Draw("col")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1500,1000)
cnv.Divide(3,2)
cnv.cd(1); histos["h_rec_vs_gen_resE"].Draw("hist")
cnv.cd(2); histos["h_rec_vs_gen_resPy"].Draw("hist")
cnv.cd(3); histos["h_rec_vs_gen_resPx"].Draw("hist")
cnv.cd(4); histos["h_rec_vs_gen_resE_vs_x"].Draw("col")
cnv.cd(5); histos["h_rec_vs_gen_resPy_vs_x"].Draw("col")
cnv.cd(6); histos["h_rec_vs_gen_resPx_vs_x"].Draw("col")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1500,1000)
cnv.Divide(3,2)
cnv.cd(1); histos["h_rec_vs_seed_resE"].Draw("hist")
cnv.cd(2); histos["h_rec_vs_seed_resPy"].Draw("hist")
cnv.cd(3); histos["h_rec_vs_seed_resPx"].Draw("hist")
cnv.cd(4); histos["h_rec_vs_seed_resE_vs_x"].Draw("col")
cnv.cd(5); histos["h_rec_vs_seed_resPy_vs_x"].Draw("col")
cnv.cd(6); histos["h_rec_vs_seed_resPx_vs_x"].Draw("col")
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

cnv = TCanvas("","",500,500)
histos["h_N_sigacc"].SetLineColor(ROOT.kBlack); histos["h_N_sigacc"].Draw()
histos["h_N_all_rec"].SetLineColor(ROOT.kBlue); histos["h_N_all_rec"].Draw("same")
histos["h_N_matched_rec"].SetLineColor(ROOT.kGreen); histos["h_N_matched_rec"].Draw("same")
histos["h_N_good_rec"].SetLineColor(ROOT.kRed); histos["h_N_good_rec"].Draw("same")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_gen_E_all"].SetLineColor(ROOT.kBlack); histos["h_gen_E_all"].Draw("hist")
histos["h_gen_E_seedMatched2gen_all"].SetLineColor(ROOT.kGreen); histos["h_gen_E_seedMatched2gen_all"].Draw("hist same")
histos["h_gen_E_seedMatched2gen_sel"].SetLineColor(ROOT.kRed);   histos["h_gen_E_seedMatched2gen_sel"].Draw("hist same")
cnv.cd(2)
histos["h_eff_E_seedMatched2gen_all"].Divide(histos["h_gen_E_seedMatched2gen_all"],histos["h_gen_E_all"])
histos["h_eff_E_seedMatched2gen_sel"].Divide(histos["h_gen_E_seedMatched2gen_sel"],histos["h_gen_E_all"])
for b in range(histos["h_eff_E_seedMatched2gen_all"].GetNbinsX()+1):
   if(histos["h_eff_E_seedMatched2gen_all"].GetBinContent(b)>1): histos["h_eff_E_seedMatched2gen_all"].SetBinContent(b,1)
   if(histos["h_eff_E_seedMatched2gen_sel"].GetBinContent(b)>1): histos["h_eff_E_seedMatched2gen_sel"].SetBinContent(b,1)
histos["h_eff_E_seedMatched2gen_all"].SetLineColor(ROOT.kGreen); histos["h_eff_E_seedMatched2gen_all"].SetMinimum(0.0); histos["h_eff_E_seedMatched2gen_all"].SetMaximum(1.1)
histos["h_eff_E_seedMatched2gen_sel"].SetLineColor(ROOT.kRed);   histos["h_eff_E_seedMatched2gen_sel"].SetMinimum(0.0); histos["h_eff_E_seedMatched2gen_sel"].SetMaximum(1.1)
histos["h_eff_E_seedMatched2gen_all"].Draw("ep")
histos["h_eff_E_seedMatched2gen_sel"].Draw("ep same")
cnv.SaveAs(fpdf)

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_gen_E_all"].SetLineColor(ROOT.kBlack); histos["h_gen_E_all"].Draw("hist")
histos["h_gen_E_recMatched2gen_all"].SetLineColor(ROOT.kGreen); histos["h_gen_E_recMatched2gen_all"].Draw("hist same")
histos["h_gen_E_recMatched2gen_sel"].SetLineColor(ROOT.kRed);   histos["h_gen_E_recMatched2gen_sel"].Draw("hist same")
cnv.cd(2)
histos["h_eff_E_recMatched2gen_all"].Divide(histos["h_gen_E_recMatched2gen_all"],histos["h_gen_E_all"])
histos["h_eff_E_recMatched2gen_sel"].Divide(histos["h_gen_E_recMatched2gen_sel"],histos["h_gen_E_all"])
for b in range(histos["h_eff_E_recMatched2gen_all"].GetNbinsX()+1):
   if(histos["h_eff_E_recMatched2gen_all"].GetBinContent(b)>1): histos["h_eff_E_recMatched2gen_all"].SetBinContent(b,1)
   if(histos["h_eff_E_recMatched2gen_sel"].GetBinContent(b)>1): histos["h_eff_E_recMatched2gen_sel"].SetBinContent(b,1)
histos["h_eff_E_recMatched2gen_all"].SetLineColor(ROOT.kGreen); histos["h_eff_E_recMatched2gen_all"].SetMinimum(0.0); histos["h_eff_E_recMatched2gen_all"].SetMaximum(1.1)
histos["h_eff_E_recMatched2gen_sel"].SetLineColor(ROOT.kRed);   histos["h_eff_E_recMatched2gen_sel"].SetMinimum(0.0); histos["h_eff_E_recMatched2gen_sel"].SetMaximum(1.1)
histos["h_eff_E_recMatched2gen_all"].Draw("ep")
histos["h_eff_E_recMatched2gen_sel"].Draw("ep same")
cnv.SaveAs(fpdf+")")





