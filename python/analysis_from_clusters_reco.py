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

process = proc
beamenergy = ebeam

### stave geometry
np=5
Hstave = 1.5  # cm
Lstave = 27 if(process=="bppp") else 50 # cm
Rbeampipe = 4 # cm
RoffsetBfield22BPPP = 7.0  # cm for BPPP in B=2.2T
RoffsetBfield20BPPP = 5.7  # cm for BPPP in B=2.0T
RoffsetBfield14BPPP = 4.0  # cm for BPPP in B=1.4T
RoffsetBfield = RoffsetBfield20BPPP if(process=="bppp") else 14 # cm
x1L = -RoffsetBfield-Lstave 
x1R = -RoffsetBfield        
x2L = +RoffsetBfield        
x2R = +RoffsetBfield+Lstave 
yUp = +Hstave/2.        
yDn = -Hstave/2.        
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

#############################################

def mapid2index(all_clusters_id):
   id2index = {}
   for i in range(all_clusters_id.size()):
      id2index.update({all_clusters_id[i]:i})
   return id2index

def acceptcls(x,y,z):
   failx = (abs(x)>x2R or abs(x)<x2L)
   if(failx): return False
   faily = (abs(y)>yUp)
   if(faily): return False
   failz = (z!=300 and z!=310 and z!=320 and z!=330)
   if(failz): return False
   return True

def accepttrk(clusters, fullacc=True):
   nlayers = 4;
   acc = 0;
   x = ROOT.Double()
   y = ROOT.Double()
   z = ROOT.Double()
   clusters[0].GetPoint(0,x,y,z)
   acc += acceptcls(x,y,z)
   clusters[1].GetPoint(0,x,y,z)
   acc += acceptcls(x,y,z)
   clusters[2].GetPoint(0,x,y,z)
   acc += acceptcls(x,y,z)
   clusters[3].GetPoint(0,x,y,z)
   acc += acceptcls(x,y,z)
   return (acc==nlayers) if(fullacc) else (acc>0)

tracks = []
points = []
sigtracks = []
sigpoints = []
bkgtracks = []
bkgpoints = []

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
   histos.update( {"h_Eee_sig":TH1D("h_Eee_sig",";E_{ee} [GeV];Tracks", nEeebins,0,Eeemax)} )
   histos.update( {"h_Eee_bkg":TH1D("h_Eee_bkg",";E_{ee} [GeV];Tracks", nEeebins,0,Eeemax)} )
   histos.update( {"h_Mee_rec":TH1D("h_Mee_rec",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   histos.update( {"h_Mee_sig":TH1D("h_Mee_sig",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   histos.update( {"h_Mee_bkg":TH1D("h_Mee_bkg",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   dEeemax = 15 if(process=="bppp") else 20
   histos.update( {"h_dEee_rec":TH1D("h_dEee_rec",";|E_{1}-E_{2}| [GeV];Tracks", 100,0,dEeemax)} )
   histos.update( {"h_dEee_sig":TH1D("h_dEee_sig",";|E_{1}-E_{2}| [GeV];Tracks", 100,0,dEeemax)} )
   histos.update( {"h_px_rec":TH1D("h_px_rec",";p_{x} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_px_sig":TH1D("h_px_sig",";p_{x} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_px_bkg":TH1D("h_px_bkg",";p_{x} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_py_rec":TH1D("h_py_rec",";p_{y} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_py_sig":TH1D("h_py_sig",";p_{y} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_py_bkg":TH1D("h_py_bkg",";p_{y} [GeV];Tracks", 200,-1,+1)} )
   histos.update( {"h_pz_rec":TH1D("h_pz_rec",";p_{z} [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_pz_sig":TH1D("h_pz_sig",";p_{z} [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_pz_bkg":TH1D("h_pz_bkg",";p_{z} [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_E_rec":TH1D("h_E_rec",";E [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_E_sig":TH1D("h_E_sig",";E [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_E_bkg":TH1D("h_E_bkg",";E [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_chi2dof_rec":TH1D("h_chi2dof_rec",";#chi^{2}/N_{DOF};Tracks",150,0,15)} )
   ntrkmax = 200  if(process=="bppp") else 1250
   ntrkmin = 0    if(process=="bppp") else 1000
   ntrkbins = 100 if(process=="bppp") else 50
   histos.update( {"h_ntrks_rec":TH1D("h_ntrks_rec",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_sig":TH1D("h_ntrks_sig",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_bkg":TH1D("h_ntrks_bkg",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_acc":TH1D("h_ntrks_acc",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )

   histos.update( {"h_res_E":TH1D("h_res_E",";(E_{rec}-E_{gen})/E_{gen};Tracks", 200,-0.05,+0.05)} )
   histos.update( {"h_res_x":TH1D("h_res_x",";(x_{rec}-x_{gen})/x_{gen};Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_res_y":TH1D("h_res_y",";(y_{rec}-y_{gen})/y_{gen};Tracks", 200,-0.25,+0.25)} )
   histos.update( {"h_res_px":TH1D("h_res_px",";(p_{x}^{rec}-p_{x}^{gen})/p_{x}^{gen};Tracks", 200,-5,+5)} )
   histos.update( {"h_res_py":TH1D("h_res_py",";(p_{y}^{rec}-p_{y}^{gen})/p_{y}^{gen};Tracks", 200,-0.5,+0.5)} )
   histos.update( {"h_res_pz":TH1D("h_res_pz",";(p_{z}^{rec}-p_{z}^{gen})/p_{z}^{gen};Tracks", 200,-0.05,+0.05)} )
   histos.update( {"h_res_pT":TH1D("h_res_pT",";(p_{T}^{rec}-p_{T}^{gen})/p_{T}^{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_theta":TH1D("h_res_theta",";(#theta_{rec}-#theta_{gen})/#theta_{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_phi":TH1D("h_res_phi",";(#phi_{rec}-#phi_{gen})/#phi_{gen};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_E_vs_E":TH2D("h_res_E_vs_E",";E_{gen} [GeV];(E_{rec}-E_{gen})/E_{gen};Tracks", 34, 0, 17, 100,-0.05,+0.05)} )
   histos.update( {"h_res_x_vs_x":TH2D("h_res_x_vs_x",";x_{gen} [cm];(x_{rec}-x_{gen})/x_{gen};Tracks", 180, -40, +40, 100,0,+0.01)} )

   histos.update( {"h_res_y_vs_x_1":TH2D("h_res_y_vs_x_1","Layer1;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_2":TH2D("h_res_y_vs_x_2","Layer2;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_3":TH2D("h_res_y_vs_x_3","Layer3;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )
   histos.update( {"h_res_y_vs_x_4":TH2D("h_res_y_vs_x_4","Layer4;(x_{rec}-x_{gen})/x_{gen};(y_{rec}-y_{gen})/y_{gen};Tracks", 100,-0.01, +0.01, 100,-0.25,+0.25)} )

   histos.update( {"h_xy_layer1_rec":TH2D("h_xy_layer1_rec","Layer 1: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_rec":TH2D("h_xy_layer2_rec","Layer 2: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_rec":TH2D("h_xy_layer3_rec","Layer 3: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_rec":TH2D("h_xy_layer4_rec","Layer 4: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_sig":TH2D("h_xy_layer1_sig","Layer 1: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_sig":TH2D("h_xy_layer2_sig","Layer 2: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_sig":TH2D("h_xy_layer3_sig","Layer 3: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_sig":TH2D("h_xy_layer4_sig","Layer 4: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_bkg":TH2D("h_xy_layer1_bkg","Layer 1: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_bkg":TH2D("h_xy_layer2_bkg","Layer 2: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_bkg":TH2D("h_xy_layer3_bkg","Layer 3: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_bkg":TH2D("h_xy_layer4_bkg","Layer 4: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer0_rec":TH2D("h_xy_layer0_rec","Dipole exit;x [cm];y [cm];Tracks", 200,-40,+40, 200,-2.,+2.)} )
   histos.update( {"h_xy_layer0_sig":TH2D("h_xy_layer0_sig","Dipole exit;x [cm];y [cm];Tracks", 200,-40,+40, 200,-2.,+2.)} )
   histos.update( {"h_xy_layer0_bkg":TH2D("h_xy_layer0_bkg","Dipole exit;x [cm];y [cm];Tracks", 200,-40,+40, 200,-2.,+2.)} )
   
   histos.update( {"h_xE_layer0_sig":TH2D("h_xE_layer0_sig","Dipole exit;x [cm];E [GeV];Tracks", 600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer1_sig":TH2D("h_xE_layer1_sig","Layer 1;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer2_sig":TH2D("h_xE_layer2_sig","Layer 2;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer3_sig":TH2D("h_xE_layer3_sig","Layer 3;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer4_sig":TH2D("h_xE_layer4_sig","Layer 4;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )

   histos.update( {"h_xE_layer0_bkg":TH2D("h_xE_layer0_bkg","Dipole exit;x [cm];E [GeV];Tracks", 600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer1_bkg":TH2D("h_xE_layer1_bkg","Layer 1;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer2_bkg":TH2D("h_xE_layer2_bkg","Layer 2;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer3_bkg":TH2D("h_xE_layer3_bkg","Layer 3;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )
   histos.update( {"h_xE_layer4_bkg":TH2D("h_xE_layer4_bkg","Layer 4;x [cm];E [GeV];Tracks",     600,-75,+75, 100,0,+20)} )

   histos.update( {"h_dx_vs_x_z0"  :TH2D("h_dx_vs_x_z0",  "At z=0 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks",   180, -0.01, +0.01, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z100":TH2D("h_dx_vs_x_z100","At z=100 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -0.075, +0.075, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z200":TH2D("h_dx_vs_x_z200","At z=200 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -40, +40, 100,-0.3,+0.3)} )
   histos.update( {"h_dx_vs_x_z300":TH2D("h_dx_vs_x_z300","At z=300 cm;x_{gen} [cm];x_{rec}-x_{gen} [cm];Tracks", 180, -40, +40, 100,-0.3,+0.3)} )

   histos.update( {"h_dx_vs_E":TH2D("h_dx_vs_E",";E_{sig} [GeV];x_{rec}-x_{sig} [cm];Tracks", 34,0,17, 100,-0.5,+0.5)} )
   histos.update( {"h_dx_vs_x":TH2D("h_dx_vs_x",";x_{sig} [cm];x_{rec}-x_{sig} [cm];Tracks", 180, -40, +40, 100,-0.5,+0.5)} )
   histos.update( {"h_dE_vs_x":TH2D("h_dE_vs_x",";x_{sig} [cm];E_{rec}-E_{sig} [GeV];Tracks", 180, -40, +40, 100,-5,+5)} )



#############################################

def GetTracks(event):
   for i in range(event.reco_trcklin.size()):
      tracks.append( event.reco_trcklin[i].Clone() )
   for i in range(event.reco_trckmar.size()):
      points.append( event.reco_trckmar[i].Clone() )


def GetSigTracks(event):
   for i in range(event.true_trcklin.size()):
      sigtracks.append( event.true_trcklin[i].Clone() )
   for i in range(event.true_trckmar.size()):
      sigpoints.append( event.true_trckmar[i].Clone() )


def GetBkgTracks(event):
   for i in range(event.bkgr_trcklin.size()):
      bkgtracks.append( event.bkgr_trcklin[i].Clone() )
   for i in range(event.bkgr_trckmar.size()):
      bkgpoints.append( event.bkgr_trckmar[i].Clone() )


def FillHistos(event):
   ntrk_sig = 0
   ntrk_bkg = 0
   ntrk_rec = 0
   ntrk_acc = 0
   sigoffset = 100000
   bkgoffset = 10000
   
   id2index = mapid2index(event.all_clusters_id)
   
   ###############################
   ### loop on signal particles
   for i in range(event.true_p.size()):
      wgt = 1 # event.true_wgt[i]
      ntrk_sig += event.true_wgt[i]
      E = event.true_p[i].E()
      px = event.true_p[i].Px()
      py = event.true_p[i].Py()
      pz = event.true_p[i].Pz()
      pT = event.true_p[i].Pt()
      theta = event.true_p[i].Theta()
      phi = event.true_p[i].Phi()
      histos["h_E_sig"].Fill(E,wgt)
      histos["h_px_sig"].Fill(px,wgt)
      histos["h_py_sig"].Fill(py,wgt)
      histos["h_pz_sig"].Fill(pz,wgt)
      ### loop over points along track (generated)
      for jxy in range(event.true_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.true_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200): histos["h_xy_layer0_sig"].Fill(x,y,wgt)
         
         if(z==200): histos["h_xE_layer0_sig"].Fill(x,E,wgt)
         if(z==300): histos["h_xE_layer1_sig"].Fill(x,E,wgt)
         if(z==310): histos["h_xE_layer2_sig"].Fill(x,E,wgt)
         if(z==320): histos["h_xE_layer3_sig"].Fill(x,E,wgt)
         if(z==330): histos["h_xE_layer4_sig"].Fill(x,E,wgt)
         
         if(z<300): continue
         if(z==300): histos["h_xy_layer1_sig"].Fill(x,y,wgt)
         if(z==310): histos["h_xy_layer2_sig"].Fill(x,y,wgt)
         if(z==320): histos["h_xy_layer3_sig"].Fill(x,y,wgt)
         if(z==330): histos["h_xy_layer4_sig"].Fill(x,y,wgt)
         
   ###############################
   ### loop on backgrouns particles
   for i in range(event.bkgr_p.size()):
      ## only in acceptance
      id1 = event.bkgr_clusters_id[i][0]
      id2 = event.bkgr_clusters_id[i][1]
      id3 = event.bkgr_clusters_id[i][2]
      id4 = event.bkgr_clusters_id[i][3]
      ix1 = id2index[id1]
      ix2 = id2index[id2]
      ix3 = id2index[id3]
      ix4 = id2index[id4]
      clusters = [event.all_clusters_xyz[ix1],event.all_clusters_xyz[ix2],event.all_clusters_xyz[ix3],event.all_clusters_xyz[ix4]]
      if(not accepttrk(clusters)): continue
      # if(not event.bkgr_acc[i]): continue
      
      wgt = 1 # event.bkgr_wgt[i]
      ntrk_bkg += event.bkgr_wgt[i]
      E = event.bkgr_p[i].E()
      px = event.bkgr_p[i].Px()
      py = event.bkgr_p[i].Py()
      pz = event.bkgr_p[i].Pz()
      pT = event.bkgr_p[i].Pt()
      theta = event.bkgr_p[i].Theta()
      phi = event.bkgr_p[i].Phi()
      histos["h_E_bkg"].Fill(E,wgt)
      histos["h_px_bkg"].Fill(px,wgt)
      histos["h_py_bkg"].Fill(py,wgt)
      histos["h_pz_bkg"].Fill(pz,wgt)
      ### loop over points along track (generated)
      for jxy in range(event.bkgr_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.bkgr_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200): histos["h_xy_layer0_bkg"].Fill(x,y,wgt)
         
         if(z==200): histos["h_xE_layer0_bkg"].Fill(x,E,wgt)
         if(z==300): histos["h_xE_layer1_bkg"].Fill(x,E,wgt)
         if(z==310): histos["h_xE_layer2_bkg"].Fill(x,E,wgt)
         if(z==320): histos["h_xE_layer3_bkg"].Fill(x,E,wgt)
         if(z==330): histos["h_xE_layer4_bkg"].Fill(x,E,wgt)
         
         if(z<300): continue
         if(z==300): histos["h_xy_layer1_bkg"].Fill(x,y,wgt)
         if(z==310): histos["h_xy_layer2_bkg"].Fill(x,y,wgt)
         if(z==320): histos["h_xy_layer3_bkg"].Fill(x,y,wgt)
         if(z==330): histos["h_xy_layer4_bkg"].Fill(x,y,wgt)

   ###################################
   ### loop on reconstructed particles
   for i in range(event.reco_p.size()):
	   ### not cutting on the acceptance yet
      wgt = 1
      ntrk_rec += 1
      Erec     = event.reco_p[i].E()
      pxrec    = event.reco_p[i].Px()
      pyrec    = event.reco_p[i].Py()
      pzrec    = event.reco_p[i].Pz()
      pTrec    = event.reco_p[i].Pt()
      thetarec = event.reco_p[i].Theta()
      phirec   = event.reco_p[i].Phi()
      chi2dof  = event.reco_chi2dof[i]
      # Tgl      = event.reco_Tgl[i]
      # Snp      = event.reco_Snp[i]
      # alpha    = event.reco_alpha[i]
      # signedinvpT = event.reco_signedinvpT[i]
      # sigmaY2     = event.reco_sigmaY2[i]
      # sigmaZY     = event.reco_sigmaZY[i]
      # sigmaZ2     = event.reco_sigmaZ2[i]
      # sigmaSnpY   = event.reco_sigmaSnpY[i]
      # sigmaSnpZ   = event.reco_sigmaSnpZ[i]
      # sigmaSnp2   = event.reco_sigmaSnp2[i]
      # sigmaTglY   = event.reco_sigmaTglY[i]
      # sigmaTglZ   = event.reco_sigmaTglZ[i]
      # sigmaTglSnp = event.reco_sigmaTglSnp[i]
      # sigmaTgl2   = event.reco_sigmaTgl2[i]
      # sigma1PtY   = event.reco_sigma1PtY[i]
      # sigma1PtZ   = event.reco_sigma1PtZ[i]
      # sigma1PtSnp = event.reco_sigma1PtSnp[i]
      # sigma1PtTgl = event.reco_sigma1PtTgl[i]
      # sigma1Pt2   = event.reco_sigma1Pt2[i]
      # invpT       = event.reco_invpT[i]
      
      
      
      
      ### loop over points along track (generated)
      for jxy in range(event.reco_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.reco_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200): histos["h_xy_layer0_rec"].Fill(x,y,wgt)
         
         if(z<300): continue
         if(z==300): histos["h_xy_layer1_rec"].Fill(x,y,wgt)
         if(z==310): histos["h_xy_layer2_rec"].Fill(x,y,wgt)
         if(z==320): histos["h_xy_layer3_rec"].Fill(x,y,wgt)
         if(z==330): histos["h_xy_layer4_rec"].Fill(x,y,wgt)
      
      
      ### fill histos
      histos["h_E_rec"].Fill(Erec,wgt)
      histos["h_px_rec"].Fill(pxrec,wgt)
      histos["h_py_rec"].Fill(pyrec,wgt)
      histos["h_pz_rec"].Fill(pzrec,wgt)
      histos["h_chi2dof_rec"].Fill(chi2dof,wgt)
      
      ismatched = event.reco_ismtchd[i]
      isig = event.reco_ixmtchd[i]
      if(not ismatched or isig<0): continue
      
      Esig     = event.true_p[isig].E()
      pxsig    = event.true_p[isig].Px()
      pysig    = event.true_p[isig].Py()
      pzsig    = event.true_p[isig].Pz()
      pTsig    = event.true_p[isig].Pt()
      thetasig = event.true_p[isig].Theta()
      phisig   = event.true_p[isig].Phi()
      
      dErel  = (Erec-Esig)/Esig
      dpxrel = (pxrec-pxsig)/pxsig
      dpyrel = (pyrec-pysig)/pysig
      dpzrel = (pzrec-pzsig)/pzsig
      dpTrel = (pTrec-pTsig)/pTsig
      dthetarel = (thetarec-thetasig)/thetasig
      dphirel = (phirec-phisig)/phisig
      
      ### fill histos
      histos["h_res_E"].Fill(dErel,wgt)
      histos["h_res_E_vs_E"].Fill(Erec,dErel,wgt)
      histos["h_res_px"].Fill(dpxrel,wgt)
      histos["h_res_py"].Fill(dpyrel,wgt)
      histos["h_res_pz"].Fill(dpzrel,wgt)
      histos["h_res_pT"].Fill(dpTrel,wgt)
      histos["h_res_phi"].Fill(dphirel,wgt)
      histos["h_res_theta"].Fill(dthetarel,wgt)
         
   ### summary filling
   histos["h_ntrks_sig"].Fill(ntrk_sig)
   histos["h_ntrks_bkg"].Fill(ntrk_bkg)
   histos["h_ntrks_rec"].Fill(ntrk_rec)

#############################################

## analysis
def Analyze(n,event):
   if(n==1): GetTracks(event)
   if(n==1): GetSigTracks(event)
   if(n==1): GetBkgTracks(event)
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
tfilename = "../data/root/rec_from_clusters_"+process+".root"
BookHistos(process)
nevents = Run(tfilename,"reco")

#############################################

## summarise
allpdf = "../output/pdf/analysis_reco_"+process+".pdf"
fn = "../output/pdf/analysis_reco_"+process+"_"


tfile = TFile("../data/root/"+process+"_geometry_truth.root","READ")
lines = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
          tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
          tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
          tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2"),
          tfile.Get("TPolyLine3D;1")]
leg_rec_sig = tfile.Get("TPave;1")

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
leg_rec_sig.Draw("same")
label("Reconstructed (tracks)",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_rec.cd()
for line in lines: line.Draw()
for trk in points: trk.Draw()
leg_rec_sig.Draw("same")
label("Reconstructed (\"clusters\")",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_rec.SaveAs(fn+"tracks.pdf")
cnv_pl3d_rec.SaveAs(allpdf+"(")
cnv_pm3d_rec.SaveAs(fn+"hits.pdf")
cnv_pm3d_rec.SaveAs(allpdf)

cnv_pl3d_sig = TCanvas("cnv_pl3d_sig","",500,500)
cnv_pl3d_sig.cd()
view_pl3d_sig = TView.CreateView(1)
view_pl3d_sig.ShowAxis()
view_pl3d_sig.SetRange(-80,-50,0, +80,+50,350)
cnv_pm3d_sig = TCanvas("cnv_pm3d_sig","",500,500)
cnv_pm3d_sig.cd()
view_pm3d_sig = TView.CreateView(1)
view_pm3d_sig.ShowAxis()
view_pm3d_sig.SetRange(-80,-50,0, +80,+50,350)
cnv_pl3d_sig.cd()
for line in lines: line.Draw()
for trk in sigtracks: trk.Draw()
for trk in bkgtracks: trk.Draw()
leg_rec_sig.Draw("same")
label("Generated (tracks)",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_sig.cd()
for line in lines: line.Draw()
for trk in sigpoints: trk.Draw()
for trk in bkgpoints: trk.Draw()
leg_rec_sig.Draw("same")
label("Generated (clusters)",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_sig.SaveAs(fn+"tracks_truth.pdf")
cnv_pl3d_sig.SaveAs(allpdf)
cnv_pm3d_sig.SaveAs(fn+"hits_truth.pdf")
cnv_pm3d_sig.SaveAs(allpdf)



cnv_xy = TCanvas("cnv_xy","",800,800)
cnv_xy.Divide(1,4)
p1 = cnv_xy.cd(1)
p2 = cnv_xy.cd(2)
p3 = cnv_xy.cd(3)
p4 = cnv_xy.cd(4)
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


cnv_xy_sig = TCanvas("cnv_xy","",800,800)
cnv_xy_sig.Divide(1,4)
p1 = cnv_xy_sig.cd(1)
p2 = cnv_xy_sig.cd(2)
p3 = cnv_xy_sig.cd(3)
p4 = cnv_xy_sig.cd(4)
p1.cd()
histos["h_xy_layer4_sig"].Scale(1./nevents)
histos["h_xy_layer4_sig"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_sig"].Scale(1./nevents)
histos["h_xy_layer3_sig"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_sig"].Scale(1./nevents)
histos["h_xy_layer2_sig"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_sig"].Scale(1./nevents)
histos["h_xy_layer1_sig"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy_sig.SaveAs(fn+"xy_sig.pdf")
cnv_xy_sig.SaveAs(allpdf)




cnv_xE_sig = TCanvas("cnv_xE","",1000,1000)
cnv_xE_sig.Divide(1,5)
p1 = cnv_xE_sig.cd(1)
p2 = cnv_xE_sig.cd(2)
p3 = cnv_xE_sig.cd(3)
p4 = cnv_xE_sig.cd(4)
p5 = cnv_xE_sig.cd(5)
p1.cd()
# histos["h_xE_layer4_sig"].Scale(1./nevents)
histos["h_xE_layer4_sig"].Draw("col")
p2.cd()
# histos["h_xE_layer3_sig"].Scale(1./nevents)
histos["h_xE_layer3_sig"].Draw("col")
p3.cd()
# histos["h_xE_layer2_sig"].Scale(1./nevents)
histos["h_xE_layer2_sig"].Draw("col")
p4.cd()
# histos["h_xE_layer1_sig"].Scale(1./nevents)
histos["h_xE_layer1_sig"].Draw("col")
p5.cd()
# histos["h_xE_layer0_sig"].Scale(1./nevents)
histos["h_xE_layer0_sig"].Draw("col")
cnv_xE_sig.SaveAs(fn+"xE_sig.pdf")
cnv_xE_sig.SaveAs(allpdf)

cnv_xy_sig = TCanvas("cnv_xy","",1000,1000)
histos["h_xy_layer0_sig"].Draw("col")
cnv_xy_sig.SaveAs(fn+"xy_sig.pdf")
cnv_xy_sig.SaveAs(allpdf)



leg_rec_sig = TLegend(0.60,0.73,0.87,0.87)
leg_rec_sig.SetFillStyle(4000) # will be transparent
leg_rec_sig.SetFillColor(0)
leg_rec_sig.SetTextFont(42)
leg_rec_sig.SetBorderSize(0)

leg_acc_rec_sig = TLegend(0.60,0.73,0.87,0.87)
leg_acc_rec_sig.SetFillStyle(4000) # will be transparent
leg_acc_rec_sig.SetFillColor(0)
leg_acc_rec_sig.SetTextFont(42)
leg_acc_rec_sig.SetBorderSize(0)

cnv_E = TCanvas("cnv_E","",500,500)
cnv_E.cd()
hmin,hmax = minmax(histos["h_E_sig"],histos["h_E_bkg"])
hmin,hmax = minmax(histos["h_E_sig"],histos["h_E_rec"])
histos["h_E_sig"].SetLineColor(ROOT.kRed)
histos["h_E_bkg"].SetLineColor(ROOT.kBlue)
histos["h_E_rec"].SetLineColor(ROOT.kBlack)
histos["h_E_sig"].Draw("hist")
histos["h_E_bkg"].Draw("hist same")
histos["h_E_rec"].Draw("hist same")
leg_rec_sig.AddEntry(histos["h_E_sig"],"Signal","l")
leg_rec_sig.AddEntry(histos["h_E_bkg"],"Background","l")
leg_rec_sig.AddEntry(histos["h_E_rec"],"Reconstructed","l")
leg_rec_sig.Draw("same")
cnv_E.SaveAs(fn+"Energy.pdf")
cnv_E.SaveAs(allpdf)


cnv_pxyz = TCanvas("cnv_E","",1500,500)
cnv_pxyz.Divide(3,1)
cnv_pxyz.cd(1)
hmin,hmax = minmax(histos["h_px_sig"],histos["h_px_bkg"])
hmin,hmax = minmax(histos["h_px_sig"],histos["h_px_rec"])
histos["h_px_sig"].SetLineColor(ROOT.kRed)
histos["h_px_bkg"].SetLineColor(ROOT.kBlue)
histos["h_px_rec"].SetLineColor(ROOT.kBlack)
histos["h_px_sig"].Draw("hist")
histos["h_px_bkg"].Draw("hist same")
histos["h_px_rec"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.cd(2)
hmin,hmax = minmax(histos["h_py_sig"],histos["h_py_bkg"])
hmin,hmax = minmax(histos["h_py_sig"],histos["h_py_rec"])
histos["h_py_sig"].SetLineColor(ROOT.kRed)
histos["h_py_bkg"].SetLineColor(ROOT.kBlue)
histos["h_py_rec"].SetLineColor(ROOT.kBlack)
histos["h_py_sig"].Draw("hist")
histos["h_py_bkg"].Draw("hist same")
histos["h_py_rec"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.cd(3)
hmin,hmax = minmax(histos["h_pz_sig"],histos["h_pz_bkg"])
hmin,hmax = minmax(histos["h_pz_sig"],histos["h_pz_rec"])
histos["h_pz_sig"].SetLineColor(ROOT.kRed)
histos["h_pz_bkg"].SetLineColor(ROOT.kBlue)
histos["h_pz_rec"].SetLineColor(ROOT.kBlack)
histos["h_pz_sig"].Draw("hist")
histos["h_pz_bkg"].Draw("hist same")
histos["h_pz_rec"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.SaveAs(fn+"Momentum.pdf")
cnv_pxyz.SaveAs(allpdf)


cnv_ntrks = TCanvas("cnv_ntrks","",500,500)
cnv_ntrks.cd()
hmin,hmax = minmax(histos["h_ntrks_sig"],histos["h_ntrks_bkg"])
hmin,hmax = minmax(histos["h_ntrks_sig"],histos["h_ntrks_rec"])
histos["h_ntrks_sig"].SetLineColor(ROOT.kRed)
histos["h_ntrks_bkg"].SetLineColor(ROOT.kBlue)
histos["h_ntrks_rec"].SetLineColor(ROOT.kBlack)
histos["h_ntrks_sig"].Draw("hist")
histos["h_ntrks_bkg"].Draw("hist same")
histos["h_ntrks_rec"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_ntrks.SaveAs(fn+"ntrks.pdf")
cnv_ntrks.SaveAs(allpdf)


# cnv_res_x = TCanvas("cnv_res_x","",500,500)
# cnv_res_x.cd()
# # histos["h_res_x"].Fit("gaus","LEM")
# histos["h_res_x"].Draw("hist")
# cnv_res_x.SaveAs(fn+"res_x.pdf")
# cnv_res_x.SaveAs(allpdf)
#
# cnv_res_y = TCanvas("cnv_res_y","",500,500)
# cnv_res_y.cd()
# # histos["h_res_y"].Fit("gaus","LEM")
# histos["h_res_y"].Draw("hist")
# cnv_res_y.SaveAs(fn+"res_y.pdf")
# cnv_res_y.SaveAs(allpdf)

cnv_chi2 = TCanvas("cnv_res_E","",500,500)
cnv_chi2.cd()
histos["h_chi2dof_rec"].Draw("hist")
cnv_chi2.SaveAs(fn+"chi2dof.pdf")
cnv_chi2.SaveAs(allpdf)

cnv_res_E = TCanvas("cnv_res_E","",500,500)
cnv_res_E.cd()
# histos["h_res_E"].Fit("gaus","LEM")
histos["h_res_E"].Draw("hist")
cnv_res_E.SaveAs(fn+"res_E.pdf")
cnv_res_E.SaveAs(allpdf)

cnv_res_pz = TCanvas("cnv_res_pz","",500,500)
cnv_res_pz.cd()
# histos["h_res_pz"].Fit("gaus","LEM")
histos["h_res_pz"].Draw("hist")
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
cnv_res_E_vs_E.SaveAs(allpdf+")")

# cnv_res_x_vs_x = TCanvas("cnv_res_x_vs_x","",500,500)
# cnv_res_x_vs_x.cd()
# histos["h_res_x_vs_x"].Draw("col")
# cnv_res_x_vs_x.SaveAs(fn+"res_x_vs_x.pdf")
# cnv_res_x_vs_x.SaveAs(allpdf)
# 
# cnv_dx_vs_E = TCanvas("cnv_dx_vs_E","",500,500)
# cnv_dx_vs_E.cd()
# histos["h_dx_vs_E"].Draw("col")
# cnv_dx_vs_E.SaveAs(fn+"dx_vs_E.pdf")
# cnv_dx_vs_E.SaveAs(allpdf)
#
# cnv_dx_vs_x = TCanvas("cnv_dx_vs_x","",500,500)
# cnv_dx_vs_x.cd()
# histos["h_dx_vs_x"].Draw("col")
# cnv_dx_vs_x.SaveAs(fn+"dx_vs_x.pdf")
# cnv_dx_vs_x.SaveAs(allpdf)
#
# cnv_dE_vs_x = TCanvas("cnv_dE_vs_x","",500,500)
# cnv_dE_vs_x.cd()
# histos["h_dE_vs_x"].Draw("col")
# cnv_dE_vs_x.SaveAs(fn+"dE_vs_x.pdf")
# cnv_dE_vs_x.SaveAs(allpdf+")")



########### write all histos to a root file
allroot = "../output/root/analysis_reco_"+process+".root"
tfileout = TFile(allroot,"RECREATE")
tfileout.cd()
for hname,hist in histos.items(): hist.Write()
tfileout.Write()
tfileout.Close()