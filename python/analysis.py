#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-e', metavar='energy', required=False,  help='beam energy')
argus = parser.parse_args()
proc  = argus.p
ebeam = argus.e

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)

#############################################

tracks = []
points = []
truthtracks = []
truthpoints = []
histos = {}

def minmax(h1,h2,f=1.1):
   hmin = h1.GetMinimum() if(h1.GetMinimum()<h2.GetMinimum()) else h2.GetMinimum()
   hmax = h1.GetMaximum() if(h1.GetMaximum()>h2.GetMaximum()) else h2.GetMaximum()
   h1.SetMaximum(hmax*f)
   h2.SetMaximum(hmax*f)
   return hmin,hmax*f

def label(text,x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1);
   s.SetTextAlign(13);
   s.SetTextColor(col);
   s.SetTextSize(0.044);
   if(boldit): text = "#font[72]{"+text+"}"
   s.DrawLatex(x,y,text);

def BookHistos(process):
   Eeemax   = 20  if(process=="bppp") else 25
   nEeebins = 200 if(process=="bppp") else 250
   histos.update( {"h_Eee_rec":TH1D("h_Eee_rec",";E_{ee} [GeV];Tracks", nEeebins,0,Eeemax)} )
   histos.update( {"h_Eee_gen":TH1D("h_Eee_gen",";E_{ee} [GeV];Tracks", nEeebins,0,Eeemax)} )
   histos.update( {"h_Mee_rec":TH1D("h_Mee_rec",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   histos.update( {"h_Mee_gen":TH1D("h_Mee_gen",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   dEeemax = 15 if(process=="bppp") else 20
   histos.update( {"h_dEee_rec":TH1D("h_dEee_rec",";|E_{1}-E_{2}| [GeV];Tracks", 100,0,dEeemax)} )
   histos.update( {"h_dEee_gen":TH1D("h_dEee_gen",";|E_{1}-E_{2}| [GeV];Tracks", 100,0,dEeemax)} )
   histos.update( {"h_E_rec":TH1D("h_E_rec",";E [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_E_gen":TH1D("h_E_gen",";E [GeV];Tracks", 200,0,20)} )
   ntrkmax = 200  if(process=="bppp") else 1250
   ntrkmin = 0    if(process=="bppp") else 1000
   ntrkbins = 100 if(process=="bppp") else 50
   histos.update( {"h_ntrks_rec":TH1D("h_ntrks_rec",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_gen":TH1D("h_ntrks_gen",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_acc":TH1D("h_ntrks_acc",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )

   histos.update( {"h_res_M":TH1D("h_res_M",";(m_{ee}^{rec}-m_{ee}^{gen})/m_{ee}^{gen};Tracks", 100,-1.5,+2.)} )
   histos.update( {"h_res_E":TH1D("h_res_E",";(E_{rec}-E_{gen})/E_{gen};Tracks", 200,-1,+1)} )
   histos.update( {"h_res_E_fine":TH1D("h_res_E_fine",";(E_{rec}-E_{gen})/E_{gen};Tracks", 200,-0.05,+0.05)} )
   histos.update( {"h_res_x":TH1D("h_res_x",";(x_{rec}-x_{gen})/x_{gen};Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_res_y":TH1D("h_res_y",";(y_{rec}-y_{gen})/y_{gen};Tracks", 200,-0.25,+0.25)} )
   histos.update( {"h_res_px":TH1D("h_res_px",";(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};Tracks", 200,-5,+5)} )
   histos.update( {"h_res_py":TH1D("h_res_py",";(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};Tracks", 200,-0.5,+0.5)} )
   histos.update( {"h_res_pz":TH1D("h_res_pz",";(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};Tracks", 200,-1,+1)} )
   histos.update( {"h_res_pz_fine":TH1D("h_res_pz_fine",";(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};Tracks", 200,-0.05,+0.05)} )
   histos.update( {"h_res_pT":TH1D("h_res_pT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_theta":TH1D("h_res_theta",";(#theta_{rec}-#theta_{gen})/#theta_{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_phi":TH1D("h_res_phi",";(#phi_{rec}-#phi_{gen})/#phi_{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_E_vs_E":TH2D("h_res_E_vs_E",";E_{gen} [GeV];(E_{rec}-E_{gen})/E_{gen};Tracks", 34, 0, 17, 100,-1,+1)} )
   histos.update( {"h_res_theta_vs_E":TH2D("h_res_theta_vs_E",";E_{gen} [GeV];(#theta_{rec}-#theta_{gen})/#theta_{gen};Tracks", 34, 0, 17, 100,-1.5,+2)} )
   histos.update( {"h_res_phi_vs_E":TH2D("h_res_phi_vs_E",";E_{gen} [GeV];(#phi_{rec}-#phi_{gen})/#phi_{gen};Tracks", 34, 0, 17, 100,-1.5,+2)} )
   histos.update( {"h_res_x_vs_x":TH2D("h_res_x_vs_x",";x_{gen} [cm];(x_{rec}-x_{gen})/x_{gen};Tracks", 180, -40, +40, 100,0,+0.01)} )
   histos.update( {"h_res_px_vs_px":TH2D("h_res_px_vs_px",";p_x^{gen} [cm];(p_x^{rec}-p_x^{gen})/p_x^{gen};Tracks", 100,-0.003,+0.003, 100,-6,+6)} ) ### new

   histos.update( {"h_res_y_vs_x_1":TH2D("h_res_y_vs_x_1","Layer1;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_2":TH2D("h_res_y_vs_x_2","Layer2;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_3":TH2D("h_res_y_vs_x_3","Layer3;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_4":TH2D("h_res_y_vs_x_4","Layer4;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )

   histos.update( {"h_xy_layer1_rec":TH2D("h_xy_layer1_rec","Layer 1: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_rec":TH2D("h_xy_layer2_rec","Layer 2: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_rec":TH2D("h_xy_layer3_rec","Layer 3: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_rec":TH2D("h_xy_layer4_rec","Layer 4: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_gen":TH2D("h_xy_layer1_gen","Layer 1: generated tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_gen":TH2D("h_xy_layer2_gen","Layer 2: generated tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_gen":TH2D("h_xy_layer3_gen","Layer 3: generated tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_gen":TH2D("h_xy_layer4_gen","Layer 4: generated tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xE_layer0_gen":TH2D("h_xE_layer0_gen","Dipole exit;x [cm];E [GeV];Tracks", 600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer1_gen":TH2D("h_xE_layer1_gen","Layer 1;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer2_gen":TH2D("h_xE_layer2_gen","Layer 2;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer3_gen":TH2D("h_xE_layer3_gen","Layer 3;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer4_gen":TH2D("h_xE_layer4_gen","Layer 4;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )

   histos.update( {"h_nITSHits":TH1D("h_nITSHits",";nITSHits;Tracks", 10,0,10)} )
   histos.update( {"h_Snp_pull":TH1D("h_Snp_pull", ";(Snp_{rec}-Snp_{gen})/#DeltaSnp_{rec};Tracks",100,-5,+5)} ) ### new
   histos.update( {"h_dSnp_vs_Snp":TH2D("h_dSnp_vs_Snp", ";Snp_{gen};Snp_{rec}-Snp_{gen};Tracks",100,-0.65e-3,+0.65e-3, 100,-1.e-3,+1.e-3)} ) ### new
   histos.update( {"h_dinvpT_vs_invpT":TH2D("h_dinvpT_vs_invpT", ";Signed 1/p_{T}^{gen};Signed (1/p_{T}^{rec}-1/p_{T}^{gen});Tracks",100,-7.e-1,+7.e-1, 100,-5.e-3,+5.e-3)} ) ### new
   histos.update( {"h_dSnp_vs_dinvpT":TH2D("h_dSnp_vs_dinvpT", ";Signed (1/p_{T}^{rec}-1/p_{T}^{gen});Snp_{rec}-Snp_{gen};Tracks",100,-5.e-3,+5.e-3, 100,-1.5e-3,+1.5e-3)} ) ### new
   histos.update( {"h_dx_vs_x_z0"  :TH2D("h_dx_vs_x_z0",  "At z=0 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks",   180, -0.01, +0.01, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z100":TH2D("h_dx_vs_x_z100","At z=100 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -0.075, +0.075, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z200":TH2D("h_dx_vs_x_z200","At z=200 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -40, +40, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z300":TH2D("h_dx_vs_x_z300","At z=300 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -40, +40, 100,-0.3,+0.3)} )

   histos.update( {"h_dx_vs_E":TH2D("h_dx_vs_E",";E_{gen} [GeV];x_{rec}-x_{gen} [cm];Tracks", 34,0,17, 100,-0.5,+0.5)} )
   histos.update( {"h_dx_vs_x":TH2D("h_dx_vs_x",";x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -40, +40, 100,-0.5,+0.5)} )
   histos.update( {"h_dE_vs_x":TH2D("h_dE_vs_x",";x_{gen} [cm];E_{rec}-E_{gen} [GeV];Tracks", 180, -40, +40, 100,-5,+5)} )



#############################################

def GetTracks(event):
   for i in range(event.poll.size()):
      tracks.append( event.poll[i].Clone() )
   for i in range(event.polm.size()):
      points.append( event.polm[i].Clone() )

def GetTruthTracks(event):
   for i in range(event.poll_gen.size()):
      truthtracks.append( event.poll_gen[i].Clone() )
   for i in range(event.polm_gen.size()):
      truthpoints.append( event.polm_gen[i].Clone() )

def FillHistos(event):
   ntrk_gen = 0
   ntrk_rec = 0
   ntrk_acc = 0
   filledgen = []
   filledrec = []
   ###############################
   ### loop on generated particles
   for i in range(event.pgen.size()):
      wgt = event.wgtgen[i]
      ntrk_gen += 1
      Egen = event.pgen[i].E()
      pxgen = event.pgen[i].Px()
      pygen = event.pgen[i].Py()
      pzgen = event.pgen[i].Pz()
      pTgen = event.pgen[i].Pt()
      thetagen = event.pgen[i].Theta()
      phigen = event.pgen[i].Phi()
      histos["h_E_gen"].Fill(Egen,wgt)
      ### electron-positron quantities (generated)
      jpair = event.jpair[i]
      if(jpair>=0): ## has to have a pair
         if(i not in filledgen and jpair not in filledgen): ## avoid double counting
            filledgen.append(i)
            filledgen.append(jpair)
            wgt12 = wgt*event.wgtgen[jpair]
            peegen = event.pgen[i]+event.pgen[jpair]
            dEee_gen = abs(event.pgen[i].E()-event.pgen[jpair].E())
            histos["h_Mee_gen"].Fill(peegen.M(),wgt12)
            histos["h_Eee_gen"].Fill(peegen.E(),wgt12)
            histos["h_dEee_gen"].Fill(dEee_gen,wgt12)
      ### loop over points along track (generated)
      for jxy in range(event.polm_gen[i].GetN()):
         xgen = ROOT.Double()
         ygen = ROOT.Double()
         zgen = ROOT.Double()
         event.polm_gen[i].GetPoint(jxy,xgen,ygen,zgen)
         
         ### noam for tracking
         if(zgen==200): histos["h_xE_layer0_gen"].Fill(xgen,Egen,wgt)
         if(zgen==300): histos["h_xE_layer1_gen"].Fill(xgen,Egen,wgt)
         if(zgen==310): histos["h_xE_layer2_gen"].Fill(xgen,Egen,wgt)
         if(zgen==320): histos["h_xE_layer3_gen"].Fill(xgen,Egen,wgt)
         if(zgen==330): histos["h_xE_layer4_gen"].Fill(xgen,Egen,wgt)
         
         if(zgen<300): continue
         if(zgen==300): histos["h_xy_layer1_gen"].Fill(xgen,ygen,wgt)
         if(zgen==310): histos["h_xy_layer2_gen"].Fill(xgen,ygen,wgt)
         if(zgen==320): histos["h_xy_layer3_gen"].Fill(xgen,ygen,wgt)
         if(zgen==330): histos["h_xy_layer4_gen"].Fill(xgen,ygen,wgt)

   ###################################
   ### loop on reconstructed particles
   for i in range(event.prec.size()):
	  ### not cutting on the acceptance yet
      igen = event.igentrk[i]
      wgt  = event.wgtgen[igen]
      ntrk_rec += 1
      for jxy in range(event.polm[i].GetN()):
         xrec = ROOT.Double()
         yrec = ROOT.Double()
         zrec = ROOT.Double()
         event.polm[i].GetPoint(jxy,xrec,yrec,zrec)
         xgen = ROOT.Double()
         ygen = ROOT.Double()
         zgen = ROOT.Double()
         event.polm_gen[igen].GetPoint(jxy,xgen,ygen,zgen)
         dxrel = (xrec-xgen)/xgen if(xgen!=0) else -999
         dyrel = (yrec-ygen)/ygen if(ygen!=0) else -999
         if(zrec<=0.1): histos["h_dx_vs_x_z0"].Fill(xgen, xrec-xgen,wgt)
         if(zrec==100): histos["h_dx_vs_x_z100"].Fill(xgen, xrec-xgen,wgt)
         if(zrec==200): histos["h_dx_vs_x_z200"].Fill(xgen, xrec-xgen,wgt)
         if(zrec==300): histos["h_dx_vs_x_z300"].Fill(xgen, xrec-xgen,wgt)
         if(zrec>=300 and zgen>=300):
            histos["h_res_x"].Fill(dxrel,wgt)
            histos["h_res_y"].Fill(dyrel,wgt)
            histos["h_res_x_vs_x"].Fill(xgen,dxrel,wgt)
            histos["h_dx_vs_E"].Fill(event.pgen[igen].E(),xrec-xgen,wgt)
            histos["h_dx_vs_x"].Fill(xgen,xrec-xgen,wgt)
            histos["h_dE_vs_x"].Fill(xgen,event.prec[i].E()-event.pgen[igen].E(),wgt)
         if(zrec==300): histos["h_res_y_vs_x_1"].Fill(dxrel,dyrel,wgt)
         if(zrec==310): histos["h_res_y_vs_x_2"].Fill(dxrel,dyrel,wgt)
         if(zrec==320): histos["h_res_y_vs_x_3"].Fill(dxrel,dyrel,wgt)
         if(zrec==330): histos["h_res_y_vs_x_4"].Fill(dxrel,dyrel,wgt)
         if(zrec==300): histos["h_xy_layer1_rec"].Fill(xrec,yrec,wgt)
         if(zrec==310): histos["h_xy_layer2_rec"].Fill(xrec,yrec,wgt)
         if(zrec==320): histos["h_xy_layer3_rec"].Fill(xrec,yrec,wgt)
         if(zrec==330): histos["h_xy_layer4_rec"].Fill(xrec,yrec,wgt)

      histos["h_nITSHits"].Fill(event.nITSHits[i],wgt)
      histos["h_Snp_pull"].Fill((event.trk_rec_Snp[i]-event.trk_gen_Snp[igen])/math.sqrt(event.trk_rec_sigmaSnp2[i]))
      histos["h_dSnp_vs_Snp"].Fill(event.trk_gen_Snp[igen], (event.trk_rec_Snp[i]-event.trk_gen_Snp[igen]))
      histos["h_dinvpT_vs_invpT"].Fill(event.trk_gen_signedinvpT[igen], (event.trk_rec_signedinvpT[i]-event.trk_gen_signedinvpT[igen]))
      histos["h_dSnp_vs_dinvpT"].Fill(event.trk_rec_signedinvpT[i]-event.trk_gen_signedinvpT[igen], event.trk_rec_Snp[i]-event.trk_gen_Snp[igen])

      ### cutting on the acceptance
      if(event.acctrkrec[i]==1):
         ntrk_acc += 1
         Erec = event.prec[i].E()
         Egen = event.pgen[igen].E()
         pxrec = event.prec[i].Px()
         pxgen = event.pgen[igen].Px()
         pyrec = event.prec[i].Py()
         pygen = event.pgen[igen].Py()
         pzrec = event.prec[i].Pz()
         pzgen = event.pgen[igen].Pz()
         pTrec = event.prec[i].Pt()
         pTgen = event.pgen[igen].Pt()
         thetarec = event.prec[i].Theta()
         thetagen = event.pgen[igen].Theta()
         phirec = event.prec[i].Phi()
         phigen = event.pgen[igen].Phi()
         dErel = (Erec-Egen)/Egen
         dpxrel = (pxrec-pxgen)/pxgen
         dpyrel = (pyrec-pygen)/pygen
         dpzrel = (pzrec-pzgen)/pzgen
         dpTrel = (pTrec-pTgen)/pTgen
         dthetarel = (thetarec-thetagen)/thetagen
         dphirel = (phirec-phigen)/phigen
         ### individual leptons
         histos["h_E_rec"].Fill(Erec,wgt)
         histos["h_res_E"].Fill(dErel,wgt)
         histos["h_res_E_fine"].Fill(dErel,wgt)
         histos["h_res_E_vs_E"].Fill(Erec,dErel,wgt)
         histos["h_res_px"].Fill(dpxrel,wgt)
         histos["h_res_px_vs_px"].Fill(pxgen,dpxrel,wgt)
         histos["h_res_py"].Fill(dpyrel,wgt)
         histos["h_res_pz"].Fill(dpzrel,wgt)
         histos["h_res_pz_fine"].Fill(dpzrel,wgt)
         histos["h_res_pT"].Fill(dpTrel,wgt)
         histos["h_res_phi"].Fill(dphirel,wgt)
         histos["h_res_theta"].Fill(dthetarel,wgt)
         histos["h_res_theta_vs_E"].Fill(Erec,dthetarel,wgt)
         histos["h_res_phi_vs_E"].Fill(Erec,dphirel,wgt)
         
         ### electron-positron quantities
         igenpair = event.jpair[igen] ### the truth partner index
         if(igenpair>=0): ### pair index should be valid
            jpair = event.jrec[igenpair] ### the reconstructed partner index
            if(jpair>=0): ### truth partner was reconstructed
               if(i not in filledrec and jpair not in filledrec): ## avoid double counting
                  filledrec.append(i)
                  filledrec.append(jpair)
                  wgt12 = wgt*event.wgtgen[igenpair]
                  peerec = event.prec[i]+event.prec[jpair]
                  peegen = event.pgen[igen]+event.pgen[igenpair]
                  dMrel  = (peerec.M()-peegen.M())/peegen.M()
                  dEee_rec = abs(event.prec[i].E()-event.prec[jpair].E())
                  dEee_gen = abs(event.pgen[igen].E()-event.pgen[igenpair].E())
                  histos["h_res_M"].Fill(dMrel,wgt12)
                  histos["h_Mee_rec"].Fill(peerec.M(),wgt12)
                  histos["h_Eee_rec"].Fill(peerec.E(),wgt12)
                  histos["h_dEee_rec"].Fill(dEee_rec,wgt12)

   ### summary filling
   histos["h_ntrks_acc"].Fill(ntrk_acc)
   histos["h_ntrks_gen"].Fill(ntrk_gen)
   histos["h_ntrks_rec"].Fill(ntrk_rec)

#############################################

## analysis
def Analyze(n,event):
   if(n==1): GetTracks(event)
   if(n==1): GetTruthTracks(event)
   FillHistos(event)
   
## event loop
def EventLoop(tree):
   nevents = tree.GetEntries()
   print("with %d events" % nevents)
   n=0 ### init n
   for event in tree:
      if(n%100==0 and n>0): print("  processed %d events" % n)
      Analyze(n,event)
      n+=1
   print("Total events processed: ",n)
   return nevents

## file and tree
def Run(tfilename,ttreename):
   print("getting tree from ",tfilename)
   tfile = TFile(tfilename,"READ")
   tree = tfile.Get(ttreename)
   nevents = EventLoop(tree)
   return nevents

#############################################

## actually run
process = proc
beamenergy = ebeam
tfilename = "../data/root/rec_"+process+".root"
BookHistos(process)
nevents = Run(tfilename,"res")

#############################################

## summarise
allpdf = "../output/pdf/all_"+process+".pdf"
fn = "../output/pdf/"+process+"_"


tfile = TFile("../data/root/geometry.root","READ")
lines = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
          tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
          tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
          tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2"),
          tfile.Get("TPolyLine3D;1")]
leg_rec_gen = tfile.Get("TPave;1")

cnv_pl3d_rec = TCanvas("cnv_pl3d_rec","",500,500)
cnv_pl3d_rec.cd()
view_pl3d_rec = TView.CreateView(1)
view_pl3d_rec.ShowAxis()
view_pl3d_rec.SetRange(-80,-50,0, +80,+50,350)
cnv_pm3d_rec = TCanvas("cnv_pm3d_rec","",500,500)
cnv_pm3d_rec.cd()
view_pm3d_rec = TView.CreateView(1)
view_pm3d_rec.ShowAxis()
view_pm3d_rec.SetRange(-80,-50,0, +80,+50,350)
cnv_pl3d_rec.cd()
for line in lines: line.Draw()
for trk in tracks: trk.Draw()
leg_rec_gen.Draw("same")
label("Reconstructed (tracks)",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_rec.cd()
for line in lines: line.Draw()
for trk in points: trk.Draw()
leg_rec_gen.Draw("same")
label("Reconstructed (\"hits\")",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_rec.SaveAs(fn+"test_tracks.pdf")
cnv_pl3d_rec.SaveAs(allpdf+"(")
cnv_pm3d_rec.SaveAs(fn+"test_hits.pdf")
cnv_pm3d_rec.SaveAs(allpdf)

cnv_pl3d_gen = TCanvas("cnv_pl3d_gen","",500,500)
cnv_pl3d_gen.cd()
view_pl3d_gen = TView.CreateView(1)
view_pl3d_gen.ShowAxis()
view_pl3d_gen.SetRange(-80,-50,0, +80,+50,350)
cnv_pm3d_gen = TCanvas("cnv_pm3d_gen","",500,500)
cnv_pm3d_gen.cd()
view_pm3d_gen = TView.CreateView(1)
view_pm3d_gen.ShowAxis()
view_pm3d_gen.SetRange(-80,-50,0, +80,+50,350)
cnv_pl3d_gen.cd()
for line in lines: line.Draw()
for trk in truthtracks: trk.Draw()
leg_rec_gen.Draw("same")
label("Generated (tracks)",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_gen.cd()
for line in lines: line.Draw()
for trk in truthpoints: trk.Draw()
leg_rec_gen.Draw("same")
label("Generated (\"hits\")",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_gen.SaveAs(fn+"test_tracks_truth.pdf")
cnv_pl3d_gen.SaveAs(allpdf)
cnv_pm3d_gen.SaveAs(fn+"test_hits_truth.pdf")
cnv_pm3d_gen.SaveAs(allpdf)



cnv_xy = TCanvas("cnv_xy","",800,800)
cnv_xy.Divide(1,4)
p1 = cnv_xy.cd(1)
p2 = cnv_xy.cd(2)
p3 = cnv_xy.cd(3)
p4 = cnv_xy.cd(4)
np=5
### stave geometry
Hstave = 1.5  # cm
Lstave = 27   # cm
Rbeampipe = 4 # cm
x1L = -Rbeampipe-Lstave # = -31
x1R = -Rbeampipe        # = -4
x2L = +Rbeampipe        # = +4
x2R = +Rbeampipe+Lstave # = +31
yUp = +Hstave/2.        # = +0.75
yDn = -Hstave/2.        # = -0.75
xLlist = [x1L,x1L,x1R,x1R,x1L]
xRlist = [x2R,x2R,x2L,x2L,x2R]
ylist  = [yDn,yUp,yUp,yDn,yDn]
xL = array.array('d', xLlist)
xR = array.array('d', xRlist)
y  = array.array('d', ylist)
staveL = TPolyLine(np,xL,y)
staveR = TPolyLine(np,xR,y)
staveL.SetLineColor(ROOT.kGreen+2)
staveR.SetLineColor(ROOT.kGreen+2)
p1.cd()
histos["h_xy_layer4_rec"].Scale(1./nevents)
histos["h_xy_layer4_rec"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_rec"].Scale(1./nevents)
histos["h_xy_layer3_rec"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_rec"].Scale(1./nevents)
histos["h_xy_layer2_rec"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_rec"].Scale(1./nevents)
histos["h_xy_layer1_rec"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy.SaveAs(fn+"xy_rec.pdf")
cnv_xy.SaveAs(allpdf)


cnv_xy_gen = TCanvas("cnv_xy","",800,800)
cnv_xy_gen.Divide(1,4)
p1 = cnv_xy_gen.cd(1)
p2 = cnv_xy_gen.cd(2)
p3 = cnv_xy_gen.cd(3)
p4 = cnv_xy_gen.cd(4)
p1.cd()
histos["h_xy_layer4_gen"].Scale(1./nevents)
histos["h_xy_layer4_gen"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_gen"].Scale(1./nevents)
histos["h_xy_layer3_gen"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_gen"].Scale(1./nevents)
histos["h_xy_layer2_gen"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_gen"].Scale(1./nevents)
histos["h_xy_layer1_gen"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy_gen.SaveAs(fn+"xy_gen.pdf")
cnv_xy_gen.SaveAs(allpdf)


cnv_res_xy = TCanvas("cnv_res_xy","",800,800)
cnv_res_xy.Divide(2,2)
p1 = cnv_res_xy.cd(1)
p2 = cnv_res_xy.cd(2)
p3 = cnv_res_xy.cd(3)
p4 = cnv_res_xy.cd(4)
p1.cd()
histos["h_res_y_vs_x_1"].Draw("col")
p2.cd()
histos["h_res_y_vs_x_2"].Draw("col")
p3.cd()
histos["h_res_y_vs_x_3"].Draw("col")
p4.cd()
histos["h_res_y_vs_x_4"].Draw("col")
cnv_res_xy.SaveAs(fn+"res_xy.pdf")
cnv_res_xy.SaveAs(allpdf)



cnv_xE_gen = TCanvas("cnv_xE","",1000,1000)
cnv_xE_gen.Divide(1,5)
p1 = cnv_xE_gen.cd(1)
p2 = cnv_xE_gen.cd(2)
p3 = cnv_xE_gen.cd(3)
p4 = cnv_xE_gen.cd(4)
p5 = cnv_xE_gen.cd(5)
p1.cd()
# histos["h_xE_layer4_gen"].Scale(1./nevents)
histos["h_xE_layer4_gen"].Draw("col")
p2.cd()
# histos["h_xE_layer3_gen"].Scale(1./nevents)
histos["h_xE_layer3_gen"].Draw("col")
p3.cd()
# histos["h_xE_layer2_gen"].Scale(1./nevents)
histos["h_xE_layer2_gen"].Draw("col")
p4.cd()
# histos["h_xE_layer1_gen"].Scale(1./nevents)
histos["h_xE_layer1_gen"].Draw("col")
p5.cd()
# histos["h_xE_layer0_gen"].Scale(1./nevents)
histos["h_xE_layer0_gen"].Draw("col")
cnv_xE_gen.SaveAs(fn+"xE_gen.pdf")
cnv_xE_gen.SaveAs(allpdf)




leg_rec_gen = TLegend(0.60,0.73,0.87,0.87)
leg_rec_gen.SetFillStyle(4000) # will be transparent
leg_rec_gen.SetFillColor(0)
leg_rec_gen.SetTextFont(42)
leg_rec_gen.SetBorderSize(0)

leg_acc_rec_gen = TLegend(0.60,0.73,0.87,0.87)
leg_acc_rec_gen.SetFillStyle(4000) # will be transparent
leg_acc_rec_gen.SetFillColor(0)
leg_acc_rec_gen.SetTextFont(42)
leg_acc_rec_gen.SetBorderSize(0)

cnv_E = TCanvas("cnv_E","",500,500)
cnv_E.cd()
hmin,hmax = minmax(histos["h_E_gen"],histos["h_E_gen"])
histos["h_E_gen"].SetLineColor(ROOT.kRed)
histos["h_E_gen"].Draw("hist")
histos["h_E_rec"].Draw("hist same")
leg_rec_gen.AddEntry(histos["h_E_gen"],"Generated","l")
leg_rec_gen.AddEntry(histos["h_E_rec"],"Reconstructed","l")
leg_rec_gen.Draw("same")
cnv_E.SaveAs(fn+"Energy.pdf")
cnv_E.SaveAs(allpdf)

cnv_Eee = TCanvas("cnv_Eee","",500,500)
cnv_Eee.cd()
hmin,hmax = minmax(histos["h_Eee_gen"],histos["h_Eee_rec"])
histos["h_Eee_gen"].SetLineColor(ROOT.kRed)
histos["h_Eee_gen"].Draw("hist")
histos["h_Eee_rec"].Draw("hist same")
leg_rec_gen.Draw("same")
cnv_Eee.SaveAs(fn+"Eee.pdf")
cnv_Eee.SaveAs(allpdf)

cnv_Mee = TCanvas("cnv_Mee","",500,500)
cnv_Mee.cd()
hmin,hmax = minmax(histos["h_Mee_gen"],histos["h_Mee_rec"])
histos["h_Mee_gen"].SetLineColor(ROOT.kRed)
histos["h_Mee_gen"].Draw("hist")
histos["h_Mee_rec"].Draw("hist same")
leg_rec_gen.Draw("same")
cnv_Mee.SaveAs(fn+"Mee.pdf")
cnv_Mee.SaveAs(allpdf)

cnv_dEee = TCanvas("cnv_dEee","",500,500)
cnv_dEee.cd()
hmin,hmax = minmax(histos["h_dEee_gen"],histos["h_dEee_rec"])
histos["h_dEee_gen"].SetLineColor(ROOT.kRed)
histos["h_dEee_gen"].Draw("hist")
histos["h_dEee_rec"].Draw("hist same")
leg_rec_gen.Draw("same")
cnv_dEee.SaveAs(fn+"dEee.pdf")
cnv_dEee.SaveAs(allpdf)

cnv_ntrks = TCanvas("cnv_ntrks","",500,500)
cnv_ntrks.cd()
hmin,hmax = minmax(histos["h_ntrks_gen"],histos["h_ntrks_rec"])
histos["h_ntrks_gen"].SetLineColor(ROOT.kRed)
histos["h_ntrks_gen"].Draw("hist")
histos["h_ntrks_rec"].Draw("hist same")
leg_rec_gen.Draw("same")
cnv_ntrks.SaveAs(fn+"ntrks.pdf")
cnv_ntrks.SaveAs(allpdf)

cnv_nITSHits = TCanvas("cnv_nITSHits","",500,500)
cnv_nITSHits.cd()
histos["h_nITSHits"].Draw("hist")
cnv_nITSHits.SaveAs(fn+"nITSHits.pdf")
cnv_nITSHits.SaveAs(allpdf)

cnv_Snp_pull = TCanvas("cnv_Snp_pull","",500,500)
cnv_Snp_pull.cd()
histos["h_Snp_pull"].Draw("hist")
cnv_Snp_pull.SaveAs(fn+"Snp_pull.pdf")
cnv_Snp_pull.SaveAs(allpdf)

cnv_dSnp_vs_Snp = TCanvas("cnv_dSnp_vs_Snp","",500,500)
cnv_dSnp_vs_Snp.cd()
histos["h_dSnp_vs_Snp"].Draw("col")
cnv_dSnp_vs_Snp.SaveAs(fn+"dSnp_vs_Snp.pdf")
cnv_dSnp_vs_Snp.SaveAs(allpdf)

cnv_dinvpT_vs_invpT = TCanvas("cnv_dinvpT_vs_invpT","",500,500)
cnv_dinvpT_vs_invpT.cd()
histos["h_dinvpT_vs_invpT"].Draw("col")
cnv_dinvpT_vs_invpT.SaveAs(fn+"dinvpT_vs_invpT.pdf")
cnv_dinvpT_vs_invpT.SaveAs(allpdf)

cnv_dSnp_vs_dinvpT = TCanvas("cnv_dSnp_vs_dinvpT","",500,500)
cnv_dSnp_vs_dinvpT.cd()
histos["h_dSnp_vs_dinvpT"].Draw("col")
cnv_dSnp_vs_dinvpT.SaveAs(fn+"dSnp_vs_dinvpT.pdf")
cnv_dSnp_vs_dinvpT.SaveAs(allpdf)

cnv_res_px_vs_px = TCanvas("cnv_res_px_vs_px","",500,500)
cnv_res_px_vs_px.cd()
histos["h_res_px_vs_px"].Draw("col")
cnv_res_px_vs_px.SaveAs(fn+"res_px_vs_px.pdf")
cnv_res_px_vs_px.SaveAs(allpdf)

cnv_dx_vs_x_z = TCanvas("cnv_dx_vs_x_z","",800,800)
cnv_dx_vs_x_z.Divide(2,2)
p1 = cnv_dx_vs_x_z.cd(1)
p2 = cnv_dx_vs_x_z.cd(2)
p3 = cnv_dx_vs_x_z.cd(3)
p4 = cnv_dx_vs_x_z.cd(4)
p1.cd()
histos["h_dx_vs_x_z0"].Draw("col")
p2.cd()
histos["h_dx_vs_x_z100"].Draw("col")
p3.cd()
histos["h_dx_vs_x_z200"].Draw("col")
p4.cd()
histos["h_dx_vs_x_z300"].Draw("col")
cnv_dx_vs_x_z.SaveAs(fn+"dx_vs_x_atz.pdf")
cnv_dx_vs_x_z.SaveAs(allpdf)

cnv_res_x = TCanvas("cnv_res_x","",500,500)
cnv_res_x.cd()
# histos["h_res_x"].Fit("gaus","LEM")
histos["h_res_x"].Draw("hist")
cnv_res_x.SaveAs(fn+"res_x.pdf")
cnv_res_x.SaveAs(allpdf)

cnv_res_y = TCanvas("cnv_res_y","",500,500)
cnv_res_y.cd()
# histos["h_res_y"].Fit("gaus","LEM")
histos["h_res_y"].Draw("hist")
cnv_res_y.SaveAs(fn+"res_y.pdf")
cnv_res_y.SaveAs(allpdf)

cnv_res_M = TCanvas("cnv_res_M","",500,500)
cnv_res_M.cd()
# histos["h_res_M"].Fit("gaus","LEM")
histos["h_res_M"].Draw("hist")
cnv_res_M.SaveAs(fn+"res_M.pdf")
cnv_res_M.SaveAs(allpdf)

cnv_res_E = TCanvas("cnv_res_E","",1000,500)
cnv_res_E.Divide(2,1)
cnv_res_E.cd(1)
# histos["h_res_E"].Fit("gaus","LEM")
histos["h_res_E"].Draw("hist")
cnv_res_E.cd(2)
# histos["h_res_E_fine"].Fit("gaus","LEM")
histos["h_res_E_fine"].Draw("hist")
cnv_res_E.SaveAs(fn+"res_E.pdf")
cnv_res_E.SaveAs(allpdf)

cnv_res_pz = TCanvas("cnv_res_pz","",1000,500)
cnv_res_pz.Divide(2,1)
cnv_res_pz.cd(1)
# histos["h_res_pz"].Fit("gaus","LEM")
histos["h_res_pz"].Draw("hist")
cnv_res_pz.cd(2)
# histos["h_res_pz_fine"].Fit("gaus","LEM")
histos["h_res_pz_fine"].Draw("hist")
cnv_res_pz.SaveAs(fn+"res_pz.pdf")
cnv_res_pz.SaveAs(allpdf)

cnv_res_px = TCanvas("cnv_res_px","",500,500)
cnv_res_px.cd()
# histos["h_res_px"].Fit("gaus","LEM")
histos["h_res_px"].Draw("hist")
cnv_res_px.SaveAs(fn+"res_px.pdf")
cnv_res_px.SaveAs(allpdf)

cnv_res_py = TCanvas("cnv_res_py","",500,500)
cnv_res_py.cd()
# histos["h_res_py"].Fit("gaus","LEM")
histos["h_res_py"].Draw("hist")
cnv_res_py.SaveAs(fn+"res_py.pdf")
cnv_res_py.SaveAs(allpdf)

cnv_res_pT = TCanvas("cnv_res_pT","",500,500)
cnv_res_pT.cd()
# histos["h_res_pT"].Fit("gaus","LEM")
histos["h_res_pT"].Draw("hist")
cnv_res_pT.SaveAs(fn+"res_pT.pdf")
cnv_res_pT.SaveAs(allpdf)

cnv_res_theta = TCanvas("cnv_res_theta","",500,500)
cnv_res_theta.cd()
# histos["h_res_theta"].Fit("gaus","LEM")
histos["h_res_theta"].Draw("hist")
cnv_res_theta.SaveAs(fn+"res_theta.pdf")
cnv_res_theta.SaveAs(allpdf)

cnv_res_phi = TCanvas("cnv_res_phi","",500,500)
cnv_res_phi.cd()
# histos["h_res_phi"].Fit("gaus","LEM")
histos["h_res_phi"].Draw("hist")
cnv_res_phi.SaveAs(fn+"res_phi.pdf")
cnv_res_phi.SaveAs(allpdf)



cnv_res_E_vs_E = TCanvas("cnv_res_E_vs_E","",500,500)
cnv_res_E_vs_E.cd()
histos["h_res_E_vs_E"].Draw("col")
cnv_res_E_vs_E.SaveAs(fn+"res_E_vs_E.pdf")
cnv_res_E_vs_E.SaveAs(allpdf)

cnv_res_theta_vs_E = TCanvas("cnv_res_theta_vs_E","",500,500)
cnv_res_theta_vs_E.cd()
histos["h_res_theta_vs_E"].Draw("col")
cnv_res_theta_vs_E.SaveAs(fn+"res_theta_vs_E.pdf")
cnv_res_theta_vs_E.SaveAs(allpdf)

cnv_res_phi_vs_E = TCanvas("cnv_res_phi_vs_E","",500,500)
cnv_res_phi_vs_E.cd()
histos["h_res_phi_vs_E"].Draw("col")
cnv_res_phi_vs_E.SaveAs(fn+"res_phi_vs_E.pdf")
cnv_res_phi_vs_E.SaveAs(allpdf)

cnv_res_x_vs_x = TCanvas("cnv_res_x_vs_x","",500,500)
cnv_res_x_vs_x.cd()
histos["h_res_x_vs_x"].Draw("col")
cnv_res_x_vs_x.SaveAs(fn+"res_x_vs_x.pdf")
cnv_res_x_vs_x.SaveAs(allpdf)

cnv_dx_vs_E = TCanvas("cnv_dx_vs_E","",500,500)
cnv_dx_vs_E.cd()
histos["h_dx_vs_E"].Draw("col")
cnv_dx_vs_E.SaveAs(fn+"dx_vs_E.pdf")
cnv_dx_vs_E.SaveAs(allpdf)

cnv_dx_vs_x = TCanvas("cnv_dx_vs_x","",500,500)
cnv_dx_vs_x.cd()
histos["h_dx_vs_x"].Draw("col")
cnv_dx_vs_x.SaveAs(fn+"dx_vs_x.pdf")
cnv_dx_vs_x.SaveAs(allpdf)

cnv_dE_vs_x = TCanvas("cnv_dE_vs_x","",500,500)
cnv_dE_vs_x.cd()
histos["h_dE_vs_x"].Draw("col")
cnv_dE_vs_x.SaveAs(fn+"dE_vs_x.pdf")
cnv_dE_vs_x.SaveAs(allpdf+")")



########### write all histos to a root file
allroot = "../output/root/all_"+process+".root"
tfileout = TFile(allroot,"RECREATE")
tfileout.cd()
for hname,hist in histos.items(): hist.Write()
tfileout.Write()
tfileout.Close()