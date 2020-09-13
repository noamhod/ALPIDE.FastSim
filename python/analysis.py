#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ctypes
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
storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")

#############################################
process = proc

### stave geometry
np=5
Hstave = 1.5  # cm
Lstave = 50 #27.12 if(process=="bppp") else 50 # cm
Rbeampipe = 2.413 # cm
RoffsetBfield = 5.7 if(process=="bppp") else 14 # cm
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

tracks = []
points = []
sigtracks = []
sigpoints = []
bkgtracks = []
bkgpoints = []
histos = {}
id2index = {}

#############################################

def Buildid2indexDict(all_clusters_id):
   for i in range(all_clusters_id.size()):
      id2index.update({all_clusters_id[i]:i})

def acceptcls(x,y,z):
   if(abs(x)>x2R):                              return False
   if(abs(x)<x2L):                              return False
   if(abs(y)>yUp):                              return False
   if(z!=300 and z!=310 and z!=320 and z!=330): return False
   return True

def accepttrk(itrk,clusters_id,clusters_xyz,fullacc=False):
   id1 = clusters_id[itrk][0]
   id2 = clusters_id[itrk][1]
   id3 = clusters_id[itrk][2]
   id4 = clusters_id[itrk][3]
   ix1 = id2index[id1] if(id1>0) else -1
   ix2 = id2index[id2] if(id2>0) else -1
   ix3 = id2index[id3] if(id3>0) else -1
   ix4 = id2index[id4] if(id4>0) else -1
   nlayers = 4;
   acc = 0;
   x = ROOT.Double()
   y = ROOT.Double()
   z = ROOT.Double()
   if(ix1>=0 and ix1<clusters_xyz.size()):
      clusters_xyz[ix1].GetPoint(0,x,y,z)
      acc += acceptcls(x,y,z)
   if(ix2>=0 and ix2<clusters_xyz.size()):
      clusters_xyz[ix2].GetPoint(0,x,y,z)
      acc += acceptcls(x,y,z)
   if(ix3>=0 and ix3<clusters_xyz.size()):
      clusters_xyz[ix3].GetPoint(0,x,y,z)
      acc += acceptcls(x,y,z)
   if(ix4>=0 and ix4<clusters_xyz.size()):
      clusters_xyz[ix4].GetPoint(0,x,y,z)
      acc += acceptcls(x,y,z)
   return (acc==nlayers) if(fullacc) else (acc>0)
   
def accepttrkpts(points,fullacc=False,isbkg=False):
   nlayers = 4;
   acc = 0;
   x = ROOT.Double()
   y = ROOT.Double()
   z = ROOT.Double()
   for i in range(points.GetN()):
      points.GetPoint(i,x,y,z)
      if(abs(x)<x2L): continue
      if(isbkg): z = ROOT.Double(ROOT.TMath.Nint(z))
      if(z!=300 and z!=310 and z!=320 and z!=330): continue
      acc += acceptcls(x,y,z)
      if(z>230): break
   return (acc==nlayers) if(fullacc) else (acc>0)

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

def GetTracks(event):
   for i in range(event.reco_trcklin.size()):
      tracks.append( event.reco_trcklin[i].Clone() )
      points.append( event.reco_trckmar[i].Clone() )

def GetSigTracks(event):
   for i in range(event.true_trcklin.size()):
      sigtracks.append( event.true_trcklin[i].Clone() )
      sigpoints.append( event.true_trckmar[i].Clone() )

def GetBkgTracks(event):
   print("Nbtrks=",event.bkgr_trcklin.size())
   for i in range(event.bkgr_trcklin.size()):
      if(not accepttrkpts(event.bkgr_trckmar[i],fullacc=False,isbkg=True)): continue
      bkgtracks.append( event.bkgr_trcklin[i].Clone() )
      bkgpoints.append( event.bkgr_trckmar[i].Clone() )
      
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

def GetSideFromX(x):
   side = ""
   if(x<0): side = "Pside"
   if(x>0): side = "Eside"
   return side
   
def GetSideFromQ(q):
   side = ""
   if(q>0): side = "Pside"
   if(q<0): side = "Eside"
   return side

def BookHistos(process):
   histos.update( {"h_px_sed":TH1D("h_px_sed",";p_{x} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_px_rec":TH1D("h_px_rec",";p_{x} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_px_sel":TH1D("h_px_sel",";p_{x} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_px_sig":TH1D("h_px_sig",";p_{x} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_px_bkg":TH1D("h_px_bkg",";p_{x} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_px_zoom_sed":TH1D("h_px_zoom_sed",";p_{x} [GeV];Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_px_zoom_rec":TH1D("h_px_zoom_rec",";p_{x} [GeV];Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_px_zoom_sel":TH1D("h_px_zoom_sel",";p_{x} [GeV];Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_px_zoom_sig":TH1D("h_px_zoom_sig",";p_{x} [GeV];Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_px_zoom_bkg":TH1D("h_px_zoom_bkg",";p_{x} [GeV];Tracks", 200,-0.01,+0.01)} )
   histos.update( {"h_py_sed":TH1D("h_py_sed",";p_{y} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_py_rec":TH1D("h_py_rec",";p_{y} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_py_sel":TH1D("h_py_sel",";p_{y} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_py_sig":TH1D("h_py_sig",";p_{y} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_py_bkg":TH1D("h_py_bkg",";p_{y} [GeV];Tracks", 200,-0.1,+0.1)} )
   histos.update( {"h_py_zoom_sed":TH1D("h_py_zoom_sed",";p_{y} [GeV];Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_py_zoom_rec":TH1D("h_py_zoom_rec",";p_{y} [GeV];Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_py_zoom_sel":TH1D("h_py_zoom_sel",";p_{y} [GeV];Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_py_zoom_sig":TH1D("h_py_zoom_sig",";p_{y} [GeV];Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_py_zoom_bkg":TH1D("h_py_zoom_bkg",";p_{y} [GeV];Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_pz_sed":TH1D("h_pz_sed",";p_{z} [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_pz_rec":TH1D("h_pz_rec",";p_{z} [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_pz_sel":TH1D("h_pz_sel",";p_{z} [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_pz_sig":TH1D("h_pz_sig",";p_{z} [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_pz_bkg":TH1D("h_pz_bkg",";p_{z} [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_E_sed":TH1D("h_E_sed",";E [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_E_rec":TH1D("h_E_rec",";E [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_E_sel":TH1D("h_E_sel",";E [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_E_sig":TH1D("h_E_sig",";E [GeV];Tracks", 170,0,17)} )
   histos.update( {"h_E_bkg":TH1D("h_E_bkg",";E [GeV];Tracks", 170,0,17)} )
   
   histos.update( {"h_chi2dof_rec"        :TH1D("h_chi2dof_rec","Inclusive;#chi^{2}/N_{DOF};Tracks",150,0,15)} )
   histos.update( {"h_chi2dof_sel"        :TH1D("h_chi2dof_sel","Selected;#chi^{2}/N_{DOF};Tracks",150,0,15)} )
   histos.update( {"h_chi2dof_rec_sig_mat":TH1D("h_chi2dof_rec_sig_mat","Signal matched;#chi^{2}/N_{DOF};Tracks",150,0,15)} )
   histos.update( {"h_chi2dof_rec_sig_non":TH1D("h_chi2dof_rec_sig_non","Non-signal mathed;#chi^{2}/N_{DOF};Tracks",150,0,15)} )
   histos.update( {"h_SnpSig_rec"         :TH1D("h_SnpSig_rec", "Inclusive;Snp/#sigma(Snp);Tracks",  150,-500,+500)} )
   histos.update( {"h_SnpSig_sel"         :TH1D("h_SnpSig_sel", "Selected;Snp/#sigma(Snp);Tracks",  150,-500,+500)} )
   histos.update( {"h_SnpSig_rec_sig_mat" :TH1D("h_SnpSig_rec_sig_mat", "Signal matched;Snp/#sigma(Snp);Tracks",  150,-500,+500)} )
   histos.update( {"h_SnpSig_rec_sig_non" :TH1D("h_SnpSig_rec_sig_non", "Non-signal matched;Snp/#sigma(Snp);Tracks",  150,-500,+500)} )
   histos.update( {"h_TglSig_rec"         :TH1D("h_TglSig_rec", "Inclusive;Tgl/#sigma(Tgl);Tracks",  150,-1000,+1000)} )   
   histos.update( {"h_TglSig_sel"         :TH1D("h_TglSig_sel", "Selected;Tgl/#sigma(Tgl);Tracks",  150,-1000,+1000)} )   
   histos.update( {"h_TglSig_rec_sig_mat" :TH1D("h_TglSig_rec_sig_mat", "Signal matched;Tgl/#sigma(Tgl);Tracks",  150,-1000,+1000)} )   
   histos.update( {"h_TglSig_rec_sig_non" :TH1D("h_TglSig_rec_sig_non", "Non-signal matched;Tgl/#sigma(Tgl);Tracks",  150,-1000,+1000)} )   
   histos.update( {"h_xVtxSig_rec"         :TH1D("h_xVtxSig_rec",         "Inclusive;x_{vtx}/#sigma(x_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_xVtxSig_sel"         :TH1D("h_xVtxSig_sel",         "Selected;x_{vtx}/#sigma(x_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_xVtxSig_rec_sig_mat" :TH1D("h_xVtxSig_rec_sig_mat", "Signal matched;x_{vtx}/#sigma(x_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_xVtxSig_rec_sig_non" :TH1D("h_xVtxSig_rec_sig_non", "Non-signal matched;x_{vtx}/#sigma(x_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_yVtxSig_rec"         :TH1D("h_yVtxSig_rec",         "Inclusive;y_{vtx}/#sigma(y_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_yVtxSig_sel"         :TH1D("h_yVtxSig_sel",         "Selected;y_{vtx}/#sigma(y_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_yVtxSig_rec_sig_mat" :TH1D("h_yVtxSig_rec_sig_mat", "Signal matched;y_{vtx}/#sigma(y_{vtx});Tracks",  150,-0.5,+0.5)} )   
   histos.update( {"h_yVtxSig_rec_sig_non" :TH1D("h_yVtxSig_rec_sig_non", "Non-signal matched;y_{vtx}/#sigma(y_{vtx});Tracks",  150,-0.5,+0.5)} )   

   ntrkmax = 450  if(process=="bppp") else 1250
   ntrkmin = 0    if(process=="bppp") else 1000
   ntrkbins = 100 if(process=="bppp") else 50
   histos.update( {"h_ntrks_rec_zoom":TH1D("h_ntrks_rec_zoom",";Track multiplicity;Events", 15,50,170)} )
   histos.update( {"h_ntrks_sel_zoom":TH1D("h_ntrks_sel_zoom",";Track multiplicity;Events", 15,50,170)} )
   histos.update( {"h_ntrks_rec":TH1D("h_ntrks_rec",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_sel":TH1D("h_ntrks_sel",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_sed":TH1D("h_ntrks_sed",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_sig":TH1D("h_ntrks_sig",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_bkg":TH1D("h_ntrks_bkg",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_acc":TH1D("h_ntrks_acc",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )

   histos.update( {"h_rat_E"     :TH1D("h_rat_E", ";E_{rec}/E_{tru};Tracks",                       200,0.98,1.02)} )
   histos.update( {"h_res_E"     :TH1D("h_res_E", ";(E_{rec}-E_{tru})/E_{tru};Tracks",             200,-0.015,+0.015)} )
   histos.update( {"h_rat_sed_E" :TH1D("h_rat_sed_E", ";E_{rec}/E_{tru};Tracks",                   200,0.98,1.02)} )
   histos.update( {"h_res_sed_E" :TH1D("h_res_sed_E", ";(E_{rec}-E_{tru})/E_{tru};Tracks",             200,-0.015,+0.015)} )
   histos.update( {"h_res_sed_pz":TH1D("h_res_sed_pz",";(p_{z}^{rec}-p_{z}^{tru})/p_{z}^{tru};Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_res_pz":TH1D("h_res_pz",";(p_{z}^{rec}-p_{z}^{tru})/p_{z}^{tru};Tracks", 200,-0.015,+0.015)} )
   histos.update( {"h_res_px":TH1D("h_res_px",";(p_{x}^{rec}-p_{x}^{tru})/p_{x}^{tru};Tracks", 200,-8,+8)} )
   histos.update( {"h_res_sed_px":TH1D("h_res_sed_px",";(p_{x}^{rec}-p_{x}^{tru})/p_{x}^{tru};Tracks", 200,-8,+8)} )
   histos.update( {"h_res_py":TH1D("h_res_py",";(p_{y}^{rec}-p_{y}^{tru})/p_{y}^{tru};Tracks", 200,-0.12,+0.12)} )
   histos.update( {"h_res_py_zoomout":TH1D("h_res_py_zoomout",";(p_{y}^{rec}-p_{y}^{tru})/p_{y}^{tru};Tracks", 100,-5,+5)} )
   histos.update( {"h_res_sed_py":TH1D("h_res_sed_py",";(p_{y}^{rec}-p_{y}^{tru})/p_{y}^{tru};Tracks", 200,-0.12,+0.12)} )
   histos.update( {"h_res_sed_py_zoomout":TH1D("h_res_sed_py_zoomout",";(p_{y}^{rec}-p_{y}^{tru})/p_{y}^{tru};Tracks", 100,-5,+5)} )
   histos.update( {"h_res_pT":TH1D("h_res_pT",";(p_{T}^{rec}-p_{T}^{tru})/p_{T}^{tru};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_sed_pT":TH1D("h_res_sed_pT",";(p_{T}^{rec}-p_{T}^{tru})/p_{T}^{tru};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_theta":TH1D("h_res_theta",";(#theta_{rec}-#theta_{tru})/#theta_{tru};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_sed_theta":TH1D("h_res_sed_theta",";(#theta_{rec}-#theta_{tru})/#theta_{tru};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_phi":TH1D("h_res_phi",";(#phi_{rec}-#phi_{tru})/#phi_{tru};Tracks", 200,-1.5,+2)} )
   histos.update( {"h_res_sed_phi":TH1D("h_res_sed_phi",";(#phi_{rec}-#phi_{tru})/#phi_{tru};Tracks", 200,-1.5,+2)} )
      
   histos.update( {"h_res_x_truvscls":TH1D("h_res_x_truvscls",";(x_{cls}-x_{tru})/x_{tru};Tracks", 200,-0.0005,+0.0005)} )
   histos.update( {"h_res_y_truvscls":TH1D("h_res_y_truvscls",";(y_{cls}-y_{tru})/y_{tru};Tracks", 200,-0.10,+0.10)} )
   histos.update( {"h_res_x_recvscls":TH1D("h_res_x_recvscls",";(x_{rec}-x_{cls})/x_{cls};Tracks", 200,-0.00015,+0.00015)} )
   histos.update( {"h_res_y_recvscls":TH1D("h_res_y_recvscls",";(y_{rec}-y_{cls})/y_{cls};Tracks", 200,-0.05,+0.05)} )
   histos.update( {"h_res_x_recvstru":TH1D("h_res_x_recvstru",";(x_{rec}-x_{tru})/x_{tru};Tracks", 200,-0.0004,+0.0004)} )
   histos.update( {"h_res_y_recvstru":TH1D("h_res_y_recvstru",";(y_{rec}-y_{tru})/y_{tru};Tracks", 200,-0.03,+0.03)} )
   
   histos.update( {"h_E_tru_sig"        :TH1D("h_E_tru_sig",        ";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_sig_acc"    :TH1D("h_E_tru_sig_acc",    ";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_rec_mat"    :TH1D("h_E_tru_rec_mat",    ";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_sel_mat"    :TH1D("h_E_tru_sel_mat",    ";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_sed_mat"    :TH1D("h_E_tru_sed_mat",    ";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_rec_mat_acc":TH1D("h_E_tru_rec_mat_acc",";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_sel_mat_acc":TH1D("h_E_tru_sel_mat_acc",";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_eff_sed"    :TH1D("h_E_tru_eff_sed",    "Seeding efficiency;#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_eff"        :TH1D("h_E_tru_eff",        "Reconstruction efficiency;#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_eff_sel"    :TH1D("h_E_tru_eff_sel",    "Reconstruction+Selection efficiency;#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_eff_acc"    :TH1D("h_E_tru_eff_acc",    "Reconstruction efficiency in-acceptance;#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   histos.update( {"h_E_tru_eff_acc_sel":TH1D("h_E_tru_eff_acc_sel","Reconstruction+Selection efficiency in acceptance;#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)})
   
   histos.update( {"h_res_x_vs_xtru_recvstru":TH2D("h_res_x_vs_xtru_recvstru",";x_{tru} [cm];(x_{rec}-x_{tru})/x_{tru};Tracks", 200,-70,+70,   200,-0.001,+0.001)} )
   histos.update( {"h_res_y_vs_ytru_recvstru":TH2D("h_res_y_vs_ytru_recvstru",";y_{tru} [cm];(y_{rec}-y_{tru})/y_{tru};Tracks", 200,-0.4,+0.4, 200,-0.075,+0.075)} )
   
   nEbins,Ebins = GetLogBinning(50,0.1,18)
   histos.update( {"h_rat_E_vs_Etru_recvstru":TH2D("h_rat_E_vs_Etru_recvstru",";E_{tru} [GeV];E_{rec}/E_{tru};Tracks",          len(Ebins)-1,array.array("d",Ebins), 200,0.8,+1.2)} )
   histos.update( {"h_res_E_vs_Etru_recvstru":TH2D("h_res_E_vs_Etru_recvstru",";E_{tru} [GeV];(E_{rec}-E_{tru})/E_{tru};Tracks",len(Ebins)-1,array.array("d",Ebins), 200,-0.02,+0.02)} )
   histos.update( {"h_dE_vs_Etru_recvstru":TH2D("h_dE_vs_Etru_recvstru",";E_{tru} [GeV];E_{rec}-E_{tru};Tracks",                len(Ebins)-1,array.array("d",Ebins), 200,-0.4,+0.4)} )
   
   ninvpbinstmp,invpbinstmp = GetLogBinning(50,0.05,1)
   invpbins = []
   for n in reversed(range(len(invpbinstmp))): invpbins.append(-invpbinstmp[n])
   for n in range(len(invpbinstmp)):           invpbins.append(+invpbinstmp[n])
   histos.update( {"h_dx_vs_invptru_recvstru_logbins":TH2D("h_dx_vs_invptru_recvstru_logbins",";1/p_{tru} [cm];x_{rec}-x_{tru} [cm];Tracks", len(invpbins)-1,array.array("d",invpbins), 200,-0.02,+0.02)} )
   
   nxbinstmp,xbinstmp = GetLogBinning(50,1,70)
   xbins = []
   for n in reversed(range(len(xbinstmp))): xbins.append(-xbinstmp[n])
   for n in range(len(xbinstmp)):           xbins.append(+xbinstmp[n])
   histos.update( {"h_res_x_vs_xtru_clsvstru_logbins":TH2D("h_res_x_vs_xtru_clsvstru_logbins",";x_{tru} [cm];(x_{cls}-x_{tru})/x_{tru};Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.001,+0.001)} )
   histos.update( {"h_res_x_vs_xcls_recvscls_logbins":TH2D("h_res_x_vs_xcls_recvscls_logbins",";x_{cls} [cm];(x_{rec}-x_{cls})/x_{cls};Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.001,+0.001)} )
   histos.update( {"h_res_x_vs_xtru_recvstru_logbins":TH2D("h_res_x_vs_xtru_recvstru_logbins",";x_{tru} [cm];(x_{rec}-x_{tru})/x_{tru};Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.001,+0.001)} )
   histos.update( {"h_dx_vs_xtru_clsvstru_logbins":TH2D("h_dx_vs_xtru_clsvstru_logbins",";x_{tru} [cm];x_{cls}-x_{tru} [cm];Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.02,+0.02)} )
   histos.update( {"h_dx_vs_xcls_recvscls_logbins":TH2D("h_dx_vs_xcls_recvscls_logbins",";x_{cls} [cm];x_{rec}-x_{cls} [cm];Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.02,+0.02)} )
   histos.update( {"h_dx_vs_xtru_recvstru_logbins":TH2D("h_dx_vs_xtru_recvstru_logbins",";x_{tru} [cm];x_{rec}-x_{tru} [cm];Tracks", len(xbins)-1,array.array("d",xbins), 200,-0.02,+0.02)} )
   
   nybinstmp,ybinstmp = GetLogBinning(50,0.001,0.4)
   ybins = []
   for n in reversed(range(len(ybinstmp))): ybins.append(-ybinstmp[n])
   for n in range(len(ybinstmp)):           ybins.append(+ybinstmp[n])
   histos.update( {"h_res_y_vs_ytru_clsvstru_logbins":TH2D("h_res_y_vs_ytru_clsvstru_logbins",";y_{tru} [cm];(y_{cls}-y_{tru})/y_{tru};Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.075,+0.075)} )  
   histos.update( {"h_res_y_vs_ycls_recvscls_logbins":TH2D("h_res_y_vs_ycls_recvscls_logbins",";y_{cls} [cm];(y_{rec}-y_{cls})/y_{cls};Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.075,+0.075)} )
   histos.update( {"h_res_y_vs_ytru_recvstru_logbins":TH2D("h_res_y_vs_ytru_recvstru_logbins",";y_{tru} [cm];(y_{rec}-y_{tru})/y_{tru};Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.075,+0.075)} )
   histos.update( {"h_dy_vs_ytru_clsvstru_logbins":TH2D("h_dy_vs_ytru_clsvstru_logbins",";y_{tru} [cm];y_{cls}-y_{tru} [cm];Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.0075,+0.0075)} )
   histos.update( {"h_dy_vs_ycls_recvscls_logbins":TH2D("h_dy_vs_ycls_recvscls_logbins",";y_{cls} [cm];y_{rec}-y_{cls} [cm];Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.0075,+0.0075)} )
   histos.update( {"h_dy_vs_ytru_recvstru_logbins":TH2D("h_dy_vs_ytru_recvstru_logbins",";y_{tru} [cm];y_{rec}-y_{tru} [cm];Tracks", len(ybins)-1,array.array("d",ybins), 200,-0.0075,+0.0075)} )

   histos.update( {"h_xy_layer1_rec":TH2D("h_xy_layer1_rec","Layer 1: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_rec":TH2D("h_xy_layer2_rec","Layer 2: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_rec":TH2D("h_xy_layer3_rec","Layer 3: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_rec":TH2D("h_xy_layer4_rec","Layer 4: reconstructed tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_sel":TH2D("h_xy_layer1_sel","Layer 1: selected tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_sel":TH2D("h_xy_layer2_sel","Layer 2: selected tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_sel":TH2D("h_xy_layer3_sel","Layer 3: selected tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_sel":TH2D("h_xy_layer4_sel","Layer 4: selected tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_sig":TH2D("h_xy_layer1_sig","Layer 1: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_sig":TH2D("h_xy_layer2_sig","Layer 2: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_sig":TH2D("h_xy_layer3_sig","Layer 3: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_sig":TH2D("h_xy_layer4_sig","Layer 4: generated sig tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_bkg":TH2D("h_xy_layer1_bkg","Layer 1: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_bkg":TH2D("h_xy_layer2_bkg","Layer 2: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_bkg":TH2D("h_xy_layer3_bkg","Layer 3: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_bkg":TH2D("h_xy_layer4_bkg","Layer 4: generated bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer1_sigbkg":TH2D("h_xy_layer1_sigbkg","Layer 1: generated sig+bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer2_sigbkg":TH2D("h_xy_layer2_sigbkg","Layer 2: generated sig+bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer3_sigbkg":TH2D("h_xy_layer3_sigbkg","Layer 3: generated sig+bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )
   histos.update( {"h_xy_layer4_sigbkg":TH2D("h_xy_layer4_sigbkg","Layer 4: generated sig+bkg tracks per event in bins of 250#times250 #mum^{2};x [cm];y [cm];Tracks", 6000,-75,+75, 80,-1,+1)} )

   histos.update( {"h_xy_layer0_rec":TH2D("h_xy_layer0_rec","Dipole exit;x [cm];y [cm];Tracks", 200,-30,+30, 200,-0.55,+0.55)} )
   histos.update( {"h_xy_layer0_sel":TH2D("h_xy_layer0_sel","Dipole exit;x [cm];y [cm];Tracks", 200,-30,+30, 200,-0.55,+0.55)} )
   histos.update( {"h_xy_layer0_sig":TH2D("h_xy_layer0_sig","Dipole exit;x [cm];y [cm];Tracks", 200,-30,+30, 200,-0.55,+0.55)} )
   histos.update( {"h_xy_layer0_bkg":TH2D("h_xy_layer0_bkg","Dipole exit;x [cm];y [cm];Tracks", 200,-30,+30, 200,-0.55,+0.55)} )
   histos.update( {"h_xy_layer0_sigbkg":TH2D("h_xy_layer0_sigbkg","Dipole exit;x [cm];y [cm];Tracks", 200,-30,+30, 200,-0.55,+0.55)} )
   
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

   ### duplicate for Eside and Pside
   hnames = list(histos.keys()).copy()
   for hname in hnames:
      for side in ["Eside","Pside"]:
         newhname = hname+"_"+side
         histos.update( {newhname: histos[hname].Clone(newhname) } )

   for name,h in histos.items(): h.Sumw2()
   


#############################################


def FillHistos(event):
   ntrk_sed = 0
   ntrk_sig = 0
   ntrk_bkg = 0
   ntrk_rec = 0
   ntrk_acc = 0
   ntrk_sel = 0
   sigoffset = 100000
   bkgoffset = 10000
   Buildid2indexDict(event.all_clusters_id)
   
   ###############################
   ### loop on seeds
   nseedsmat = 0
   for i in range(event.seed_p.size()):
      if(process=="trident" and event.seed_q[i]<0): continue
      wgt = 1 # event.seed_wgt[i]
      ntrk_sed   += wgt
      seedtype   = event.seed_type[i]
      clusterids = event.seed_clusters_id[i]
      seed_q     = event.seed_q[i]
      seed_p     = event.seed_p[i]
      side       = GetSideFromQ( seed_q )
      
      histos["h_E_sed"].Fill(seed_p.E(),wgt)
      histos["h_px_sed"].Fill(seed_p.Px(),wgt)
      histos["h_px_zoom_sed"].Fill(seed_p.Px(),wgt)
      histos["h_py_sed"].Fill(seed_p.Py(),wgt)
      histos["h_py_zoom_sed"].Fill(seed_p.Py(),wgt)
      histos["h_pz_sed"].Fill(seed_p.Pz(),wgt)
      
      histos["h_E_sed_"+side].Fill(seed_p.E(),wgt)
      histos["h_px_sed_"+side].Fill(seed_p.Px(),wgt)
      histos["h_px_zoom_sed_"+side].Fill(seed_p.Px(),wgt)
      histos["h_py_sed_"+side].Fill(seed_p.Py(),wgt)
      histos["h_py_zoom_sed_"+side].Fill(seed_p.Py(),wgt)
      histos["h_pz_sed_"+side].Fill(seed_p.Pz(),wgt)
      
      
      if(event.seed_type[i]==1):
         seed_true_id = event.seed_clusters_id[i][0] ## id of cluster at layer 1
         seed_true_ix = -1
         for j in range(event.true_clusters_id.size()):
            if(event.true_clusters_id[j][0]==seed_true_id): ## id of cluster at layer 1
               seed_true_ix = j
               break
         if(seed_true_ix>=0):
            nseedsmat += 1
            histos["h_E_tru_sed_mat"].Fill(event.true_p[seed_true_ix].E(),wgt)
            histos["h_E_tru_sed_mat_"+side].Fill(event.true_p[seed_true_ix].E(),wgt)
            
            
            E_sig     = event.true_p[seed_true_ix].E()
            px_sig    = event.true_p[seed_true_ix].Px()
            py_sig    = event.true_p[seed_true_ix].Py()
            pz_sig    = event.true_p[seed_true_ix].Pz()
            pT_sig    = event.true_p[seed_true_ix].Pt()
            theta_sig = event.true_p[seed_true_ix].Theta()
            phi_sig   = event.true_p[seed_true_ix].Phi()
      
            sed_Erat   = seed_p.E()/E_sig if(E_sig>0) else -9999
            sed_dErel  = (seed_p.E()-E_sig)/E_sig
            sed_dpxrel = (seed_p.Px()-px_sig)/px_sig
            sed_dpyrel = (seed_p.Py()-py_sig)/py_sig
            sed_dpzrel = (seed_p.Pz()-pz_sig)/pz_sig
            sed_dpTrel = (seed_p.Pt()-pT_sig)/pT_sig
            sed_dthetarel = (seed_p.Theta()-theta_sig)/theta_sig
            sed_dphirel = (seed_p.Phi()-phi_sig)/phi_sig
            
            # print("sed_dpyrel=%g, seed_p.Py()=%g, py_sig=%g" %(sed_dpyrel,seed_p.Py(),py_sig))
      
            histos["h_rat_sed_E"].Fill(sed_Erat,wgt)
            histos["h_res_sed_E"].Fill(sed_dErel,wgt)
            histos["h_res_sed_px"].Fill(sed_dpxrel,wgt)
            histos["h_res_sed_py"].Fill(sed_dpyrel,wgt)
            histos["h_res_sed_py_zoomout"].Fill(sed_dpyrel,wgt)
            histos["h_res_sed_pz"].Fill(sed_dpzrel,wgt)
            histos["h_res_sed_pT"].Fill(sed_dpTrel,wgt)
            histos["h_res_sed_phi"].Fill(sed_dphirel,wgt)
            histos["h_res_sed_theta"].Fill(sed_dthetarel,wgt)
            
            histos["h_rat_sed_E_"+side].Fill(sed_Erat,wgt)
            histos["h_res_sed_E_"+side].Fill(sed_dErel,wgt)
            histos["h_res_sed_px_"+side].Fill(sed_dpxrel,wgt)
            histos["h_res_sed_py_"+side].Fill(sed_dpyrel,wgt)
            histos["h_res_sed_py_zoomout_"+side].Fill(sed_dpyrel,wgt)
            histos["h_res_sed_pz_"+side].Fill(sed_dpzrel,wgt)
            histos["h_res_sed_pT_"+side].Fill(sed_dpTrel,wgt)
            histos["h_res_sed_phi_"+side].Fill(sed_dphirel,wgt)
            histos["h_res_sed_theta_"+side].Fill(sed_dthetarel,wgt)
   
   
   ###############################
   ### loop on signal particles
   for i in range(event.true_p.size()):
      if(process=="trident" and event.true_q[i]<0): continue
      wgt = 1 # event.true_wgt[i]
      ntrk_sig += wgt
      E = event.true_p[i].E()
      px = event.true_p[i].Px()
      py = event.true_p[i].Py()
      pz = event.true_p[i].Pz()
      pT = event.true_p[i].Pt()
      theta = event.true_p[i].Theta()
      phi = event.true_p[i].Phi()
      side = GetSideFromQ( event.true_q[i] )
      
      ### get the true clusters
      cls_tru_id1 = event.true_clusters_id[i][0]
      cls_tru_id2 = event.true_clusters_id[i][1]
      cls_tru_id3 = event.true_clusters_id[i][2]
      cls_tru_id4 = event.true_clusters_id[i][3]
      cls_tru_ix1 = id2index[cls_tru_id1] if(cls_tru_id1>=0) else -1
      cls_tru_ix2 = id2index[cls_tru_id2] if(cls_tru_id2>=0) else -1
      cls_tru_ix3 = id2index[cls_tru_id3] if(cls_tru_id3>=0) else -1
      cls_tru_ix4 = id2index[cls_tru_id4] if(cls_tru_id4>=0) else -1
      cls_tru_r1  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_tru_r2  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_tru_r3  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_tru_r4  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      if(cls_tru_ix1>=0 and cls_tru_ix1<event.all_clusters_xyz.size()): event.all_clusters_xyz[cls_tru_ix1].GetPoint(0,cls_tru_r1[0],cls_tru_r1[1],cls_tru_r1[2])
      else:               cls_tru_r1 = [-999,-999,-999]
      if(cls_tru_ix2>=0 and cls_tru_ix2<event.all_clusters_xyz.size()): event.all_clusters_xyz[cls_tru_ix2].GetPoint(0,cls_tru_r2[0],cls_tru_r2[1],cls_tru_r2[2])
      else:               cls_tru_r2 = [-999,-999,-999] 
      if(cls_tru_ix3>=0 and cls_tru_ix3<event.all_clusters_xyz.size()): event.all_clusters_xyz[cls_tru_ix3].GetPoint(0,cls_tru_r3[0],cls_tru_r3[1],cls_tru_r3[2])
      else:               cls_tru_r3 = [-999,-999,-999] 
      if(cls_tru_ix4>=0 and cls_tru_ix4<event.all_clusters_xyz.size()): event.all_clusters_xyz[cls_tru_ix4].GetPoint(0,cls_tru_r4[0],cls_tru_r4[1],cls_tru_r4[2])
      else:               cls_tru_r4 = [-999,-999,-999]
      
      histos["h_E_sig"].Fill(E,wgt)
      histos["h_px_sig"].Fill(px,wgt)
      histos["h_px_zoom_sig"].Fill(px,wgt)
      histos["h_py_sig"].Fill(py,wgt)
      histos["h_py_zoom_sig"].Fill(py,wgt)
      histos["h_pz_sig"].Fill(pz,wgt)
      histos["h_E_tru_sig"].Fill(E,wgt)
      if(accepttrk(i,event.true_clusters_id,event.all_clusters_xyz,fullacc=True)): histos["h_E_tru_sig_acc"].Fill(E,wgt)
      
      histos["h_E_sig_"+side].Fill(E,wgt)
      histos["h_px_sig_"+side].Fill(px,wgt)
      histos["h_px_zoom_sig_"+side].Fill(px,wgt)
      histos["h_py_sig_"+side].Fill(py,wgt)
      histos["h_py_zoom_sig_"+side].Fill(py,wgt)
      histos["h_pz_sig_"+side].Fill(pz,wgt)
      histos["h_E_tru_sig_"+side].Fill(E,wgt)
      if(accepttrk(i,event.true_clusters_id,event.all_clusters_xyz,fullacc=True)): histos["h_E_tru_sig_acc_"+side].Fill(E,wgt)
      
      
      ### loop over points along track (generated)
      for jxy in range(event.true_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.true_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200):
            histos["h_xy_layer0_sig"].Fill(x,y,wgt)
            histos["h_xy_layer0_sigbkg"].Fill(x,y,wgt)
            histos["h_xE_layer0_sig"].Fill(x,E,wgt)
            
            histos["h_xy_layer0_sig_"+side].Fill(x,y,wgt)
            histos["h_xy_layer0_sigbkg_"+side].Fill(x,y,wgt)
            histos["h_xE_layer0_sig_"+side].Fill(x,E,wgt)
         if(z==300):
            histos["h_xE_layer1_sig"].Fill(x,E,wgt)
            histos["h_xy_layer1_sig"].Fill(x,y,wgt)
            histos["h_xy_layer1_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xE_layer1_sig_"+side].Fill(x,E,wgt)
            histos["h_xy_layer1_sig_"+side].Fill(x,y,wgt)
            histos["h_xy_layer1_sigbkg_"+side].Fill(x,y,wgt)
            
            dx = (cls_tru_r1[0]-x)
            dy = (cls_tru_r1[1]-y)
            dxrel = dx/x if(cls_tru_r1[0]!=-999 and x!=0) else -999
            dyrel = dy/y if(cls_tru_r1[1]!=-999 and y!=0) else -999
            histos["h_res_x_truvscls"].Fill(dxrel,wgt)
            histos["h_res_y_truvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins"].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins"].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins"].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins"].Fill(y,dy,wgt)
            
            histos["h_res_x_truvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_truvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins_"+side].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins_"+side].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins_"+side].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins_"+side].Fill(y,dy,wgt)
            
         if(z==310):
            histos["h_xE_layer2_sig"].Fill(x,E,wgt)
            histos["h_xy_layer2_sig"].Fill(x,y,wgt)
            histos["h_xy_layer2_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xE_layer2_sig_"+side].Fill(x,E,wgt)
            histos["h_xy_layer2_sig_"+side].Fill(x,y,wgt)
            histos["h_xy_layer2_sigbkg_"+side].Fill(x,y,wgt)
            
            dx = (cls_tru_r2[0]-x)
            dy = (cls_tru_r2[1]-y)
            dxrel = dx/x if(cls_tru_r2[0]!=-999 and x!=0) else -999
            dyrel = dy/y if(cls_tru_r2[1]!=-999 and y!=0) else -999
            histos["h_res_x_truvscls"].Fill(dxrel,wgt)
            histos["h_res_y_truvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins"].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins"].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins"].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins"].Fill(y,dy,wgt)
            
            histos["h_res_x_truvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_truvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins_"+side].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins_"+side].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins_"+side].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins_"+side].Fill(y,dy,wgt)
            
         if(z==320):
            histos["h_xE_layer3_sig"].Fill(x,E,wgt)
            histos["h_xy_layer3_sig"].Fill(x,y,wgt)
            histos["h_xy_layer3_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xE_layer3_sig_"+side].Fill(x,E,wgt)
            histos["h_xy_layer3_sig_"+side].Fill(x,y,wgt)
            histos["h_xy_layer3_sigbkg_"+side].Fill(x,y,wgt)
            
            dx = (cls_tru_r3[0]-x)
            dy = (cls_tru_r3[1]-y)
            dxrel = dx/x if(cls_tru_r3[0]!=-999 and x!=0) else -999
            dyrel = dy/y if(cls_tru_r3[1]!=-999 and y!=0) else -999
            histos["h_res_x_truvscls"].Fill(dxrel,wgt)
            histos["h_res_y_truvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins"].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins"].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins"].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins"].Fill(y,dy,wgt)
            
            histos["h_res_x_truvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_truvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins_"+side].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins_"+side].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins_"+side].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins_"+side].Fill(y,dy,wgt)
            
         if(z==330):
            histos["h_xE_layer4_sig"].Fill(x,E,wgt)
            histos["h_xy_layer4_sig"].Fill(x,y,wgt)
            histos["h_xy_layer4_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xE_layer4_sig_"+side].Fill(x,E,wgt)
            histos["h_xy_layer4_sig_"+side].Fill(x,y,wgt)
            histos["h_xy_layer4_sigbkg_"+side].Fill(x,y,wgt)
            
            dx = (cls_tru_r4[0]-x)
            dy = (cls_tru_r4[1]-y)
            dxrel = dx/x if(cls_tru_r4[0]!=-999 and x!=0) else -999
            dyrel = dy/y if(cls_tru_r4[1]!=-999 and y!=0) else -999
            histos["h_res_x_truvscls"].Fill(dxrel,wgt)
            histos["h_res_y_truvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins"].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins"].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins"].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins"].Fill(y,dy,wgt)
            
            histos["h_res_x_truvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_truvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_clsvstru_logbins_"+side].Fill(x,dxrel,wgt)
            histos["h_res_y_vs_ytru_clsvstru_logbins_"+side].Fill(y,dyrel,wgt)
            histos["h_dx_vs_xtru_clsvstru_logbins_"+side].Fill(x,dx,wgt)
            histos["h_dy_vs_ytru_clsvstru_logbins_"+side].Fill(y,dy,wgt)
         
   ###############################
   ### loop on backgrouns particles
   for i in range(event.bkgr_p.size()):
      if(process=="trident"):
         ### loop over points along track (generated) ad check if in Eside (charge test is not good here!)
         inEside = False
         for jxy in range(event.bkgr_trckmar[i].GetN()):
            x = ROOT.Double()
            y = ROOT.Double()
            z = ROOT.Double()
            event.bkgr_trckmar[i].GetPoint(jxy,x,y,z)
            if(z==300 and x>0):
               inEside = True
               break
         if(inEside): continue
            
            
      ## only in acceptance
      if(not accepttrk(i,event.bkgr_clusters_id,event.all_clusters_xyz,fullacc=False)): continue
      if(not accepttrkpts(event.bkgr_trckmar[i],fullacc=False,isbkg=True)):             continue
      wgt = 1 # event.bkgr_wgt[i]
      ntrk_bkg += wgt
      E = event.bkgr_p[i].E()
      px = event.bkgr_p[i].Px()
      py = event.bkgr_p[i].Py()
      pz = event.bkgr_p[i].Pz()
      pT = event.bkgr_p[i].Pt()
      theta = event.bkgr_p[i].Theta()
      phi = event.bkgr_p[i].Phi()
      side = GetSideFromQ( event.bkgr_q[i] )
      
      histos["h_E_bkg"].Fill(E,wgt)
      histos["h_px_bkg"].Fill(px,wgt)
      histos["h_px_zoom_bkg"].Fill(px,wgt)
      histos["h_py_bkg"].Fill(py,wgt)
      histos["h_py_zoom_bkg"].Fill(py,wgt)
      histos["h_pz_bkg"].Fill(pz,wgt)
      
      histos["h_E_bkg_"+side].Fill(E,wgt)
      histos["h_px_bkg_"+side].Fill(px,wgt)
      histos["h_px_zoom_bkg_"+side].Fill(px,wgt)
      histos["h_py_bkg_"+side].Fill(py,wgt)
      histos["h_py_zoom_bkg_"+side].Fill(py,wgt)
      histos["h_pz_bkg_"+side].Fill(pz,wgt)
      
      ### loop over points along track (generated)
      for jxy in range(event.bkgr_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.bkgr_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200):
            histos["h_xy_layer0_bkg"].Fill(x,y,wgt)
            histos["h_xy_layer0_sigbkg"].Fill(x,y,wgt)
         
         if(z==200): histos["h_xE_layer0_bkg"].Fill(x,E,wgt)
         if(z==300): histos["h_xE_layer1_bkg"].Fill(x,E,wgt)
         if(z==310): histos["h_xE_layer2_bkg"].Fill(x,E,wgt)
         if(z==320): histos["h_xE_layer3_bkg"].Fill(x,E,wgt)
         if(z==330): histos["h_xE_layer4_bkg"].Fill(x,E,wgt)
         
         if(z==200): histos["h_xE_layer0_bkg_"+side].Fill(x,E,wgt)
         if(z==300): histos["h_xE_layer1_bkg_"+side].Fill(x,E,wgt)
         if(z==310): histos["h_xE_layer2_bkg_"+side].Fill(x,E,wgt)
         if(z==320): histos["h_xE_layer3_bkg_"+side].Fill(x,E,wgt)
         if(z==330): histos["h_xE_layer4_bkg_"+side].Fill(x,E,wgt)
         
         if(z<300): continue
         if(z==300):
            histos["h_xy_layer1_bkg"].Fill(x,y,wgt)
            histos["h_xy_layer1_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xy_layer1_bkg_"+side].Fill(x,y,wgt)
            histos["h_xy_layer1_sigbkg_"+side].Fill(x,y,wgt)
         if(z==310):
            histos["h_xy_layer2_bkg"].Fill(x,y,wgt)
            histos["h_xy_layer2_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xy_layer2_bkg_"+side].Fill(x,y,wgt)
            histos["h_xy_layer2_sigbkg_"+side].Fill(x,y,wgt)
         if(z==320):
            histos["h_xy_layer3_bkg"].Fill(x,y,wgt)
            histos["h_xy_layer3_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xy_layer3_bkg_"+side].Fill(x,y,wgt)
            histos["h_xy_layer3_sigbkg_"+side].Fill(x,y,wgt)
         if(z==330):
            histos["h_xy_layer4_bkg"].Fill(x,y,wgt)
            histos["h_xy_layer4_sigbkg"].Fill(x,y,wgt)
            
            histos["h_xy_layer4_bkg_"+side].Fill(x,y,wgt)
            histos["h_xy_layer4_sigbkg_"+side].Fill(x,y,wgt)


   ###################################
   ### loop on reconstructed particles
   for i in range(event.reco_p.size()):
      if(process=="trident" and event.reco_q[i]<0): continue
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
      Tgl         = event.reco_Tgl[i]
      Snp         = event.reco_Snp[i] ## the slope in X direction: probe->GetTrack()->GetSnp()
      alpha       = event.reco_alpha[i]
      signedinvpT = event.reco_signedinvpT[i] ## the curvature (q/Pyz): probe->GetTrack()->GetSigned1Pt()
      sigmaY2     = event.reco_sigmaY2[i]
      sigmaZY     = event.reco_sigmaZY[i]
      sigmaZ2     = event.reco_sigmaZ2[i]
      sigmaSnpY   = event.reco_sigmaSnpY[i]
      sigmaSnpZ   = event.reco_sigmaSnpZ[i]
      sigmaSnp2   = event.reco_sigmaSnp2[i] ## probe->GetTrack()->GetSigmaSnp2()
      sigmaTglY   = event.reco_sigmaTglY[i]
      sigmaTglZ   = event.reco_sigmaTglZ[i]
      sigmaTglSnp = event.reco_sigmaTglSnp[i]
      sigmaTgl2   = event.reco_sigmaTgl2[i]
      sigma1PtY   = event.reco_sigma1PtY[i]
      sigma1PtZ   = event.reco_sigma1PtZ[i]
      sigma1PtSnp = event.reco_sigma1PtSnp[i]
      sigma1PtTgl = event.reco_sigma1PtTgl[i]
      sigma1Pt2   = event.reco_sigma1Pt2[i]
      invpT       = event.reco_invpT[i]
      signedpT    = event.reco_signedpT[i]
      xVtx        = event.reco_x[i]
      yVtx        = event.reco_y[i]
      side        = GetSideFromQ( event.reco_q[i] )
      
      ### significances
      SnpSig  = Snp/math.sqrt(sigmaSnp2) if(sigmaSnp2>0) else -1e10
      TglSig  = Tgl/math.sqrt(sigmaTgl2) if(sigmaTgl2>0) else -1e10
      xVtxSig = xVtx/math.sqrt(sigmaY2)  if(sigmaY2>0)   else -1e10
      yVtxSig = yVtx/math.sqrt(sigmaZ2)  if(sigmaZ2>0)   else -1e10
      
      ### selection?
      selected = True
      if(chi2dof>6):         selected = False
      if(abs(xVtxSig)>0.01): selected = False
      if(abs(yVtxSig)>0.01): selected = False
      if(abs(SnpSig)>10):    selected = False
      if(abs(TglSig)>400):   selected = False
      if(Erec<1.5):          selected = False
      if(Erec>15.5):         selected = False
      if(selected): ntrk_sel += 1
      
      ### get the reco clusters
      cls_rec_id1 = event.reco_clusters_id[i][0]
      cls_rec_id2 = event.reco_clusters_id[i][1]
      cls_rec_id3 = event.reco_clusters_id[i][2]
      cls_rec_id4 = event.reco_clusters_id[i][3]
      cls_rec_ix1 = id2index[cls_rec_id1] if(cls_rec_id1>=0) else -1
      cls_rec_ix2 = id2index[cls_rec_id2] if(cls_rec_id2>=0) else -1
      cls_rec_ix3 = id2index[cls_rec_id3] if(cls_rec_id3>=0) else -1
      cls_rec_ix4 = id2index[cls_rec_id4] if(cls_rec_id4>=0) else -1
      cls_rec_r1  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_rec_r2  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_rec_r3  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      cls_rec_r4  = [ ROOT.Double(),ROOT.Double(),ROOT.Double() ]
      if(cls_rec_ix1>=0): event.all_clusters_xyz[cls_rec_ix1].GetPoint(0,cls_rec_r1[0],cls_rec_r1[1],cls_rec_r1[2])
      else:               cls_rec_r1 = [-999,-999,-999]
      if(cls_rec_ix2>=0): event.all_clusters_xyz[cls_rec_ix2].GetPoint(0,cls_rec_r2[0],cls_rec_r2[1],cls_rec_r2[2])
      else:               cls_rec_r2 = [-999,-999,-999] 
      if(cls_rec_ix3>=0): event.all_clusters_xyz[cls_rec_ix3].GetPoint(0,cls_rec_r3[0],cls_rec_r3[1],cls_rec_r3[2])
      else:               cls_rec_r3 = [-999,-999,-999] 
      if(cls_rec_ix4>=0): event.all_clusters_xyz[cls_rec_ix4].GetPoint(0,cls_rec_r4[0],cls_rec_r4[1],cls_rec_r4[2])
      else:               cls_rec_r4 = [-999,-999,-999]
      
      
      ### loop over points along the reco track
      for jxy in range(event.reco_trckmar[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.reco_trckmar[i].GetPoint(jxy,x,y,z)
         
         if(z==200):
            histos["h_xy_layer0_rec"].Fill(x,y,wgt)
            histos["h_xy_layer0_rec_"+side].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer0_sel"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer0_sel_"+side].Fill(x,y,wgt)
         if(z==300):
            histos["h_xy_layer1_rec"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer1_sel"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer1_sel_"+side].Fill(x,y,wgt)
            dx = (x-cls_rec_r1[0])
            dy = (y-cls_rec_r1[1])
            dxrel = dx/cls_rec_r1[0] if(cls_rec_r1[0]!=-999 and cls_rec_r1[0]!=0) else -999
            dyrel = dy/cls_rec_r1[1] if(cls_rec_r1[1]!=-999 and cls_rec_r1[1]!=0) else -999
            histos["h_res_x_recvscls"].Fill(dxrel,wgt)
            histos["h_res_y_recvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins"].Fill(cls_rec_r1[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins"].Fill(cls_rec_r1[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins"].Fill(cls_rec_r1[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins"].Fill(cls_rec_r1[1],dy,wgt)
            
            histos["h_res_x_recvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_recvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r1[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r1[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r1[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r1[1],dy,wgt)
         if(z==310):
            histos["h_xy_layer2_rec"].Fill(x,y,wgt)
            histos["h_xy_layer2_rec_"+side].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer2_sel"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer2_sel_"+side].Fill(x,y,wgt)
            dx = (x-cls_rec_r2[0])
            dy = (y-cls_rec_r2[1])
            dxrel = dx/cls_rec_r2[0] if(cls_rec_r2[0]!=-999 and cls_rec_r2[0]!=0) else -999
            dyrel = dy/cls_rec_r2[1] if(cls_rec_r2[1]!=-999 and cls_rec_r2[1]!=0) else -999
            histos["h_res_x_recvscls"].Fill(dxrel,wgt)
            histos["h_res_y_recvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins"].Fill(cls_rec_r2[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins"].Fill(cls_rec_r2[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins"].Fill(cls_rec_r2[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins"].Fill(cls_rec_r2[1],dy,wgt)
            
            histos["h_res_x_recvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_recvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r2[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r2[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r2[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r2[1],dy,wgt)
         if(z==320):
            histos["h_xy_layer3_rec"].Fill(x,y,wgt)
            histos["h_xy_layer3_rec_"+side].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer3_sel"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer3_sel_"+side].Fill(x,y,wgt)
            dx = (x-cls_rec_r3[0])
            dy = (y-cls_rec_r3[1])
            dxrel = dx/cls_rec_r3[0] if(cls_rec_r3[0]!=-999 and cls_rec_r3[0]!=0) else -999
            dyrel = dy/cls_rec_r3[1] if(cls_rec_r3[1]!=-999 and cls_rec_r3[1]!=0) else -999
            histos["h_res_x_recvscls"].Fill(dxrel,wgt)
            histos["h_res_y_recvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins"].Fill(cls_rec_r3[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins"].Fill(cls_rec_r3[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins"].Fill(cls_rec_r3[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins"].Fill(cls_rec_r3[1],dy,wgt)
            
            histos["h_res_x_recvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_recvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r3[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r3[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r3[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r3[1],dy,wgt)
         if(z==330):
            histos["h_xy_layer4_rec"].Fill(x,y,wgt)
            histos["h_xy_layer4_rec_"+side].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer4_sel"].Fill(x,y,wgt)
            if(selected): histos["h_xy_layer4_sel_"+side].Fill(x,y,wgt)
            dx = (x-cls_rec_r4[0])
            dy = (y-cls_rec_r4[1])
            dxrel = dx/cls_rec_r4[0] if(cls_rec_r4[0]!=-999 and cls_rec_r4[0]!=0) else -999
            dyrel = dy/cls_rec_r4[1] if(cls_rec_r4[1]!=-999 and cls_rec_r4[1]!=0) else -999
            histos["h_res_x_recvscls"].Fill(dxrel,wgt)
            histos["h_res_y_recvscls"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins"].Fill(cls_rec_r4[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins"].Fill(cls_rec_r4[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins"].Fill(cls_rec_r4[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins"].Fill(cls_rec_r4[1],dy,wgt)
            
            histos["h_res_x_recvscls_"+side].Fill(dxrel,wgt)
            histos["h_res_y_recvscls_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r4[0],dxrel,wgt)
            histos["h_res_y_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r4[1],dyrel,wgt)
            histos["h_dx_vs_xcls_recvscls_logbins_"+side].Fill(cls_rec_r4[0],dx,wgt)
            histos["h_dy_vs_ycls_recvscls_logbins_"+side].Fill(cls_rec_r4[1],dy,wgt)
      
      
      ### fill histos
      histos["h_E_rec"].Fill(Erec,wgt)
      histos["h_px_rec"].Fill(pxrec,wgt)
      histos["h_px_zoom_rec"].Fill(pxrec,wgt)
      histos["h_py_rec"].Fill(pyrec,wgt)
      histos["h_py_zoom_rec"].Fill(pyrec,wgt)
      histos["h_pz_rec"].Fill(pzrec,wgt)
      histos["h_chi2dof_rec"].Fill(chi2dof,wgt)
      histos["h_SnpSig_rec"].Fill(SnpSig,wgt)
      histos["h_TglSig_rec"].Fill(TglSig,wgt)
      histos["h_xVtxSig_rec"].Fill(xVtxSig,wgt)
      histos["h_yVtxSig_rec"].Fill(yVtxSig,wgt)
      
      histos["h_E_rec_"+side].Fill(Erec,wgt)
      histos["h_px_rec_"+side].Fill(pxrec,wgt)
      histos["h_px_zoom_rec_"+side].Fill(pxrec,wgt)
      histos["h_py_rec_"+side].Fill(pyrec,wgt)
      histos["h_py_zoom_rec_"+side].Fill(pyrec,wgt)
      histos["h_pz_rec_"+side].Fill(pzrec,wgt)
      histos["h_chi2dof_rec_"+side].Fill(chi2dof,wgt)
      histos["h_SnpSig_rec_"+side].Fill(SnpSig,wgt)
      histos["h_TglSig_rec_"+side].Fill(TglSig,wgt)
      histos["h_xVtxSig_rec_"+side].Fill(xVtxSig,wgt)
      histos["h_yVtxSig_rec_"+side].Fill(yVtxSig,wgt)
      
      if(selected):
         histos["h_E_sel"].Fill(Erec,wgt)
         histos["h_px_sel"].Fill(pxrec,wgt)
         histos["h_px_zoom_sel"].Fill(pxrec,wgt)
         histos["h_py_sel"].Fill(pyrec,wgt)
         histos["h_py_zoom_sel"].Fill(pyrec,wgt)
         histos["h_pz_sel"].Fill(pzrec,wgt)
         histos["h_chi2dof_sel"].Fill(chi2dof,wgt)
         histos["h_SnpSig_sel"].Fill(SnpSig,wgt)
         histos["h_TglSig_sel"].Fill(TglSig,wgt)
         histos["h_xVtxSig_sel"].Fill(xVtxSig,wgt)
         histos["h_yVtxSig_sel"].Fill(yVtxSig,wgt)
         
         histos["h_E_sel_"+side].Fill(Erec,wgt)
         histos["h_px_sel_"+side].Fill(pxrec,wgt)
         histos["h_px_zoom_sel_"+side].Fill(pxrec,wgt)
         histos["h_py_sel_"+side].Fill(pyrec,wgt)
         histos["h_py_zoom_sel_"+side].Fill(pyrec,wgt)
         histos["h_pz_sel_"+side].Fill(pzrec,wgt)
         histos["h_chi2dof_sel_"+side].Fill(chi2dof,wgt)
         histos["h_SnpSig_sel_"+side].Fill(SnpSig,wgt)
         histos["h_TglSig_sel_"+side].Fill(TglSig,wgt)
         histos["h_xVtxSig_sel_"+side].Fill(xVtxSig,wgt)
         histos["h_yVtxSig_sel_"+side].Fill(yVtxSig,wgt)
      
      ismatched = event.reco_ismtchd[i]
      isig = event.reco_ixmtchd[i]
      if(not ismatched or isig<0):
         histos["h_chi2dof_rec_sig_non"].Fill(chi2dof,wgt)
         histos["h_SnpSig_rec_sig_non"].Fill(SnpSig,wgt)
         histos["h_TglSig_rec_sig_non"].Fill(TglSig,wgt)
         histos["h_xVtxSig_rec_sig_non"].Fill(xVtxSig,wgt)
         histos["h_yVtxSig_rec_sig_non"].Fill(yVtxSig,wgt)
         
         histos["h_chi2dof_rec_sig_non_"+side].Fill(chi2dof,wgt)
         histos["h_SnpSig_rec_sig_non_"+side].Fill(SnpSig,wgt)
         histos["h_TglSig_rec_sig_non_"+side].Fill(TglSig,wgt)
         histos["h_xVtxSig_rec_sig_non_"+side].Fill(xVtxSig,wgt)
         histos["h_yVtxSig_rec_sig_non_"+side].Fill(yVtxSig,wgt)
         continue
      
      ### loop over points along the reco track
      for jxy_rec in range(event.reco_trckmar[i].GetN()):
         xrec = ROOT.Double()
         yrec = ROOT.Double()
         zrec = ROOT.Double()
         event.reco_trckmar[i].GetPoint(jxy_rec,xrec,yrec,zrec)
         if(zrec!=300 and zrec!=310 and zrec!=320 and zrec!=330): continue
         for jxy_tru in range(event.true_trckmar[isig].GetN()):
            xtru = ROOT.Double()
            ytru = ROOT.Double()
            ztru = ROOT.Double()
            event.true_trckmar[isig].GetPoint(jxy_tru,xtru,ytru,ztru)
            if(ztru!=zrec): continue
            invptru = 1./event.true_p[isig].P()*(-1) if(xtru>0) else 1./event.true_p[isig].P()*(+1)
            dx = (xrec-xtru)
            dy = (yrec-ytru)
            dxrel = dx/xtru if(xtru!=0) else -999
            dyrel = dy/ytru if(ytru!=0) else -999
            histos["h_res_x_recvstru"].Fill(dxrel,wgt)
            histos["h_res_y_recvstru"].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_recvstru"].Fill(xtru,dxrel,wgt)
            histos["h_res_y_vs_ytru_recvstru"].Fill(ytru,dyrel,wgt)
            histos["h_res_x_vs_xtru_recvstru_logbins"].Fill(xtru,dxrel,wgt)
            histos["h_res_y_vs_ytru_recvstru_logbins"].Fill(ytru,dyrel,wgt)
            histos["h_dx_vs_invptru_recvstru_logbins"].Fill(invptru,dx,wgt)
            histos["h_dx_vs_xtru_recvstru_logbins"].Fill(xtru,dx,wgt)
            histos["h_dy_vs_ytru_recvstru_logbins"].Fill(ytru,dy,wgt)
            
            histos["h_res_x_recvstru_"+side].Fill(dxrel,wgt)
            histos["h_res_y_recvstru_"+side].Fill(dyrel,wgt)
            histos["h_res_x_vs_xtru_recvstru_"+side].Fill(xtru,dxrel,wgt)
            histos["h_res_y_vs_ytru_recvstru_"+side].Fill(ytru,dyrel,wgt)
            histos["h_res_x_vs_xtru_recvstru_logbins_"+side].Fill(xtru,dxrel,wgt)
            histos["h_res_y_vs_ytru_recvstru_logbins_"+side].Fill(ytru,dyrel,wgt)
            histos["h_dx_vs_invptru_recvstru_logbins_"+side].Fill(invptru,dx,wgt)
            histos["h_dx_vs_xtru_recvstru_logbins_"+side].Fill(xtru,dx,wgt)
            histos["h_dy_vs_ytru_recvstru_logbins_"+side].Fill(ytru,dy,wgt)
         
      
      Esig     = event.true_p[isig].E()
      pxsig    = event.true_p[isig].Px()
      pysig    = event.true_p[isig].Py()
      pzsig    = event.true_p[isig].Pz()
      pTsig    = event.true_p[isig].Pt()
      thetasig = event.true_p[isig].Theta()
      phisig   = event.true_p[isig].Phi()
      
      Erat   = Erec/Esig if(Esig>0) else -9999
      dE     = Erec-Esig
      dErel  = dE/Esig
      dpxrel = (pxrec-pxsig)/pxsig
      dpyrel = (pyrec-pysig)/pysig
      dpzrel = (pzrec-pzsig)/pzsig
      dpTrel = (pTrec-pTsig)/pTsig
      dthetarel = (thetarec-thetasig)/thetasig
      dphirel = (phirec-phisig)/phisig
      
      ### fill histos
      histos["h_chi2dof_rec_sig_mat"].Fill(chi2dof,wgt)
      histos["h_SnpSig_rec_sig_mat"].Fill(SnpSig,wgt)
      histos["h_TglSig_rec_sig_mat"].Fill(TglSig,wgt)
      histos["h_xVtxSig_rec_sig_mat"].Fill(xVtxSig,wgt)
      histos["h_yVtxSig_rec_sig_mat"].Fill(yVtxSig,wgt)
      histos["h_E_tru_rec_mat"].Fill(Esig)
      
      histos["h_chi2dof_rec_sig_mat_"+side].Fill(chi2dof,wgt)
      histos["h_SnpSig_rec_sig_mat_"+side].Fill(SnpSig,wgt)
      histos["h_TglSig_rec_sig_mat_"+side].Fill(TglSig,wgt)
      histos["h_xVtxSig_rec_sig_mat_"+side].Fill(xVtxSig,wgt)
      histos["h_yVtxSig_rec_sig_mat_"+side].Fill(yVtxSig,wgt)
      histos["h_E_tru_rec_mat_"+side].Fill(Esig)
      
      if(selected): histos["h_E_tru_sel_mat"].Fill(Esig)
      if(selected): histos["h_E_tru_sel_mat_"+side].Fill(Esig)
      if(accepttrk(isig,event.true_clusters_id,event.all_clusters_xyz,fullacc=True)):
         histos["h_E_tru_rec_mat_acc"].Fill(Esig)
         histos["h_E_tru_rec_mat_acc_"+side].Fill(Esig)
         if(selected): histos["h_E_tru_sel_mat_acc"].Fill(Esig)
         if(selected): histos["h_E_tru_sel_mat_acc_"+side].Fill(Esig)
      histos["h_rat_E"].Fill(Erat,wgt)
      histos["h_res_E"].Fill(dErel,wgt)
      histos["h_rat_E_vs_Etru_recvstru"].Fill(Esig,Erat,wgt)
      histos["h_res_E_vs_Etru_recvstru"].Fill(Esig,dErel,wgt)
      histos["h_dE_vs_Etru_recvstru"].Fill(Esig,dE,wgt)
      histos["h_res_px"].Fill(dpxrel,wgt)
      histos["h_res_py"].Fill(dpyrel,wgt)
      histos["h_res_py_zoomout"].Fill(dpyrel,wgt)
      histos["h_res_pz"].Fill(dpzrel,wgt)
      histos["h_res_pT"].Fill(dpTrel,wgt)
      histos["h_res_phi"].Fill(dphirel,wgt)
      histos["h_res_theta"].Fill(dthetarel,wgt)
      
      histos["h_rat_E_"+side].Fill(Erat,wgt)
      histos["h_res_E_"+side].Fill(dErel,wgt)
      histos["h_rat_E_vs_Etru_recvstru_"+side].Fill(Esig,Erat,wgt)
      histos["h_res_E_vs_Etru_recvstru_"+side].Fill(Esig,dErel,wgt)
      histos["h_dE_vs_Etru_recvstru_"+side].Fill(Esig,dE,wgt)
      histos["h_res_px_"+side].Fill(dpxrel,wgt)
      histos["h_res_py_"+side].Fill(dpyrel,wgt)
      histos["h_res_py_zoomout_"+side].Fill(dpyrel,wgt)
      histos["h_res_pz_"+side].Fill(dpzrel,wgt)
      histos["h_res_pT_"+side].Fill(dpTrel,wgt)
      histos["h_res_phi_"+side].Fill(dphirel,wgt)
      histos["h_res_theta_"+side].Fill(dthetarel,wgt)
         
   ### summary filling
   histos["h_ntrks_sig"].Fill(ntrk_sig)
   histos["h_ntrks_bkg"].Fill(ntrk_bkg)
   histos["h_ntrks_sed"].Fill(ntrk_sed)
   histos["h_ntrks_rec"].Fill(ntrk_rec)
   histos["h_ntrks_sel"].Fill(ntrk_sel)
   histos["h_ntrks_sel_zoom"].Fill(ntrk_sel)
   histos["h_ntrks_rec_zoom"].Fill(ntrk_rec)
   
   histos["h_ntrks_sig_"+side].Fill(ntrk_sig)
   histos["h_ntrks_bkg_"+side].Fill(ntrk_bkg)
   histos["h_ntrks_sed_"+side].Fill(ntrk_sed)
   histos["h_ntrks_rec_"+side].Fill(ntrk_rec)
   histos["h_ntrks_sel_"+side].Fill(ntrk_sel)
   histos["h_ntrks_sel_zoom_"+side].Fill(ntrk_sel)
   histos["h_ntrks_rec_zoom_"+side].Fill(ntrk_rec)

#############################################

## analysis
def Analyze(n,event):
   if(n==0): GetTracks(event)
   if(n==0): GetSigTracks(event)
   if(n==0): GetBkgTracks(event)
   FillHistos(event)
   
## event loop
def EventLoop(tree):
   nevents = tree.GetEntries()
   print("with %d events" % nevents)
   n=0 ### init n
   for event in tree:
      if(n%10==0 and n>0): print("  processed %d events" % n)
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
tfilename = storage+"/data/root/rec_"+process+".root"
BookHistos(process)
nevents = Run(tfilename,"reco")

#############################################

## summarise
allpdf = storage+"/output/pdf/analysis_rec_"+process+".pdf"
fn = storage+"/output/pdf/analysis_rec_"+process+"_"


tfile = TFile(storage+"/data/root/"+process+"_geometry_truth.root","READ")
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
cnv_pl3d_rec.SaveAs(fn+"rec_tracks.pdf")
cnv_pl3d_rec.SaveAs(allpdf+"(")
cnv_pm3d_rec.SaveAs(fn+"rec_clusters.pdf")
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
leg_rec_sig.Draw("same")
label("Signal tracks",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_sig.cd()
for line in lines: line.Draw()
for trk in sigpoints: trk.Draw()
leg_rec_sig.Draw("same")
label("Signal clusters",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_sig.SaveAs(fn+"sig_tracks.pdf")
cnv_pl3d_sig.SaveAs(allpdf)
cnv_pm3d_sig.SaveAs(fn+"sig_clusters.pdf")
cnv_pm3d_sig.SaveAs(allpdf)

cnv_pl3d_bkg = TCanvas("cnv_pl3d_bkg","",500,500)
cnv_pl3d_bkg.cd()
view_pl3d_bkg = TView.CreateView(1)
view_pl3d_bkg.ShowAxis()
view_pl3d_bkg.SetRange(-80,-50,0, +80,+50,350)
cnv_pm3d_bkg = TCanvas("cnv_pm3d_bkg","",500,500)
cnv_pm3d_bkg.cd()
view_pm3d_bkg = TView.CreateView(1)
view_pm3d_bkg.ShowAxis()
view_pm3d_bkg.SetRange(-80,-50,0, +80,+50,350)
cnv_pl3d_bkg.cd()
for line in lines: line.Draw()
for trk in bkgtracks: trk.Draw()
leg_rec_sig.Draw("same")
label("Background tracks",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_bkg.cd()
for line in lines: line.Draw()
for trk in bkgpoints: trk.Draw()
leg_rec_sig.Draw("same")
label("Background clusters",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_bkg.SaveAs(fn+"bkg_tracks.pdf")
cnv_pl3d_bkg.SaveAs(allpdf)
cnv_pm3d_bkg.SaveAs(fn+"bkg_clusters.pdf")
cnv_pm3d_bkg.SaveAs(allpdf)

cnv_pl3d_sigbkg = TCanvas("cnv_pl3d_sigbkg","",500,500)
cnv_pl3d_sigbkg.cd()
view_pl3d_sigbkg = TView.CreateView(1)
view_pl3d_sigbkg.ShowAxis()
view_pl3d_sigbkg.SetRange(-80,-50,0, +80,+50,350)
cnv_pm3d_sigbkg = TCanvas("cnv_pm3d_sigbkg","",500,500)
cnv_pm3d_sigbkg.cd()
view_pm3d_sigbkg = TView.CreateView(1)
view_pm3d_sigbkg.ShowAxis()
view_pm3d_sigbkg.SetRange(-80,-50,0, +80,+50,350)
cnv_pl3d_sigbkg.cd()
for line in lines: line.Draw()
for trk in sigtracks: trk.Draw()
for trk in bkgtracks: trk.Draw()
leg_rec_sig.Draw("same")
label("Signal+Background tracks",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pm3d_sigbkg.cd()
for line in lines: line.Draw()
for trk in sigpoints: trk.Draw()
for trk in bkgpoints: trk.Draw()
leg_rec_sig.Draw("same")
label("Signal+Background clusters",0.4,0.8)
label("Staves",0.62,0.74,ROOT.kGreen+2)
label("Dipole",0.8,0.45,ROOT.kGray)
cnv_pl3d_sigbkg.SaveAs(fn+"sigbkg_tracks.pdf")
cnv_pl3d_sigbkg.SaveAs(allpdf)
cnv_pm3d_sigbkg.SaveAs(fn+"sigbkg_clusters.pdf")
cnv_pm3d_sigbkg.SaveAs(allpdf)



cnv_xy = TCanvas("cnv_xy_rec","",800,800)
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

cnv_xy = TCanvas("cnv_xy_sel","",800,800)
cnv_xy.Divide(1,4)
p1 = cnv_xy.cd(1)
p2 = cnv_xy.cd(2)
p3 = cnv_xy.cd(3)
p4 = cnv_xy.cd(4)
p1.cd()
histos["h_xy_layer4_sel"].Scale(1./nevents)
histos["h_xy_layer4_sel"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_sel"].Scale(1./nevents)
histos["h_xy_layer3_sel"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_sel"].Scale(1./nevents)
histos["h_xy_layer2_sel"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_sel"].Scale(1./nevents)
histos["h_xy_layer1_sel"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy.SaveAs(fn+"xy_sel.pdf")
cnv_xy.SaveAs(allpdf)

cnv_xy_sig = TCanvas("cnv_xy_sig","",800,800)
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


cnv_xy_bkg = TCanvas("cnv_xy_bkg","",800,800)
cnv_xy_bkg.Divide(1,4)
p1 = cnv_xy_bkg.cd(1)
p2 = cnv_xy_bkg.cd(2)
p3 = cnv_xy_bkg.cd(3)
p4 = cnv_xy_bkg.cd(4)
p1.cd()
histos["h_xy_layer4_bkg"].Scale(1./nevents)
histos["h_xy_layer4_bkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_bkg"].Scale(1./nevents)
histos["h_xy_layer3_bkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_bkg"].Scale(1./nevents)
histos["h_xy_layer2_bkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_bkg"].Scale(1./nevents)
histos["h_xy_layer1_bkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy_bkg.SaveAs(fn+"xy_bkg.pdf")
cnv_xy_bkg.SaveAs(allpdf)


cnv_xy_sigbkg = TCanvas("cnv_xy_sigbkg","",800,800)
cnv_xy_sigbkg.Divide(1,4)
p1 = cnv_xy_sigbkg.cd(1)
p2 = cnv_xy_sigbkg.cd(2)
p3 = cnv_xy_sigbkg.cd(3)
p4 = cnv_xy_sigbkg.cd(4)
p1.cd()
histos["h_xy_layer4_sigbkg"].Scale(1./nevents)
histos["h_xy_layer4_sigbkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p2.cd()
histos["h_xy_layer3_sigbkg"].Scale(1./nevents)
histos["h_xy_layer3_sigbkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p3.cd()
histos["h_xy_layer2_sigbkg"].Scale(1./nevents)
histos["h_xy_layer2_sigbkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
p4.cd()
histos["h_xy_layer1_sigbkg"].Scale(1./nevents)
histos["h_xy_layer1_sigbkg"].Draw("colz")
staveL.Draw()
staveR.Draw()
cnv_xy_sigbkg.SaveAs(fn+"xy_sigbkg.pdf")
cnv_xy_sigbkg.SaveAs(allpdf)


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


cnv_xy_sig = TCanvas("cnv_xy0_sig","",1000,1000)
histos["h_xy_layer0_sig"].Draw("col")
cnv_xy_sig.SaveAs(fn+"xy0_sig.pdf")
cnv_xy_sig.SaveAs(allpdf)

cnv_xy_bkg = TCanvas("cnv_xy0_bkg","",1000,1000)
histos["h_xy_layer0_bkg"].Draw("col")
cnv_xy_bkg.SaveAs(fn+"xy0_bkg.pdf")
cnv_xy_bkg.SaveAs(allpdf)

cnv_xy_sigbkg = TCanvas("cnv_xy0_sigbkg","",1000,1000)
histos["h_xy_layer0_sigbkg"].Draw("col")
cnv_xy_sigbkg.SaveAs(fn+"xy0_sigbkg.pdf")
cnv_xy_sigbkg.SaveAs(allpdf)


leg_rec_sig = TLegend(0.60,0.73,0.87,0.87)
leg_rec_sig.SetFillStyle(4000) # will be transparent
leg_rec_sig.SetFillColor(0)
leg_rec_sig.SetTextFont(42)
leg_rec_sig.SetBorderSize(0)

leg_bkg_sig = TLegend(0.60,0.73,0.87,0.87)
leg_bkg_sig.SetFillStyle(4000) # will be transparent
leg_bkg_sig.SetFillColor(0)
leg_bkg_sig.SetTextFont(42)
leg_bkg_sig.SetBorderSize(0)

leg_acc_rec_sig = TLegend(0.60,0.73,0.87,0.87)
leg_acc_rec_sig.SetFillStyle(4000) # will be transparent
leg_acc_rec_sig.SetFillColor(0)
leg_acc_rec_sig.SetTextFont(42)
leg_acc_rec_sig.SetBorderSize(0)



cnv_E = TCanvas("cnv_E","",500,500)
cnv_E.cd()
cnv_E.SetLogy()
hmin,hmax = minmax(histos["h_E_sig"],histos["h_E_bkg"])
hmin,hmax = minmax(histos["h_E_sig"],histos["h_E_rec"])
hmin,hmax = minmax(histos["h_E_rec"],histos["h_E_sed"])
hmin,hmax = minmax(histos["h_E_sed"],histos["h_E_sig"])

hmin,hmax = minmax(histos["h_E_sig_Eside"],histos["h_E_bkg_Eside"])
hmin,hmax = minmax(histos["h_E_sig_Eside"],histos["h_E_rec_Eside"])
hmin,hmax = minmax(histos["h_E_rec_Eside"],histos["h_E_sed_Eside"])
hmin,hmax = minmax(histos["h_E_sed_Eside"],histos["h_E_sig_Eside"])

hmin,hmax = minmax(histos["h_E_sig_Pside"],histos["h_E_bkg_Pside"])
hmin,hmax = minmax(histos["h_E_sig_Pside"],histos["h_E_rec_Pside"])
hmin,hmax = minmax(histos["h_E_rec_Pside"],histos["h_E_sed_Pside"])
hmin,hmax = minmax(histos["h_E_sed_Pside"],histos["h_E_sig_Pside"])

histos["h_E_sig"].SetLineColor(ROOT.kRed)
histos["h_E_bkg"].SetLineColor(ROOT.kBlue)
histos["h_E_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_E_rec"].SetLineColor(ROOT.kBlack)
histos["h_E_sig"].Draw("hist")
histos["h_E_bkg"].Draw("hist same")
histos["h_E_rec"].Draw("hist same")
histos["h_E_sed"].Draw("hist same")
leg_rec_sig.AddEntry(histos["h_E_sig"],"Signal","l")
leg_rec_sig.AddEntry(histos["h_E_bkg"],"Background","l")
leg_rec_sig.AddEntry(histos["h_E_sed"],"Seeds","l")
leg_rec_sig.AddEntry(histos["h_E_rec"],"Reconstructed","l")
leg_rec_sig.Draw("same")
cnv_E.SaveAs(fn+"Energy.pdf")
cnv_E.SaveAs(allpdf)

cnv_E_norec = TCanvas("cnv_E_norec","",500,500)
cnv_E_norec.cd()
cnv_E_norec.SetLogy()
hmin,hmax = minmax(histos["h_E_sig"],histos["h_E_bkg"])
hmin,hmax = minmax(histos["h_E_sig_Eside"],histos["h_E_bkg_Eside"])
hmin,hmax = minmax(histos["h_E_sig_Pside"],histos["h_E_bkg_Pside"])
histos["h_E_sig"].SetLineColor(ROOT.kRed)
histos["h_E_bkg"].SetLineColor(ROOT.kBlue)
histos["h_E_sig"].Draw("hist")
histos["h_E_bkg"].Draw("hist same")
leg_bkg_sig.AddEntry(histos["h_E_sig"],"Signal","l")
leg_bkg_sig.AddEntry(histos["h_E_bkg"],"Background","l")
leg_bkg_sig.Draw("same")
cnv_E_norec.SaveAs(fn+"Energy_norec.pdf")
cnv_E_norec.SaveAs(allpdf)

cnv_pxyz = TCanvas("cnv_pxyz","",1500,500)
cnv_pxyz.Divide(3,1)
cnv_pxyz.cd(1)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_px_sig"],histos["h_px_bkg"])
hmin,hmax = minmax(histos["h_px_sig"],histos["h_px_rec"])
hmin,hmax = minmax(histos["h_px_rec"],histos["h_px_sed"])
hmin,hmax = minmax(histos["h_px_sed"],histos["h_px_sig"])

hmin,hmax = minmax(histos["h_px_sig_Eside"],histos["h_px_bkg_Eside"])
hmin,hmax = minmax(histos["h_px_sig_Eside"],histos["h_px_rec_Eside"])
hmin,hmax = minmax(histos["h_px_rec_Eside"],histos["h_px_sed_Eside"])
hmin,hmax = minmax(histos["h_px_sed_Eside"],histos["h_px_sig_Eside"])

hmin,hmax = minmax(histos["h_px_sig_Pside"],histos["h_px_bkg_Pside"])
hmin,hmax = minmax(histos["h_px_sig_Pside"],histos["h_px_rec_Pside"])
hmin,hmax = minmax(histos["h_px_rec_Pside"],histos["h_px_sed_Pside"])
hmin,hmax = minmax(histos["h_px_sed_Pside"],histos["h_px_sig_Pside"])

histos["h_px_sig"].SetLineColor(ROOT.kRed)
histos["h_px_bkg"].SetLineColor(ROOT.kBlue)
histos["h_px_rec"].SetLineColor(ROOT.kBlack)
histos["h_px_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_px_sig"].Draw("hist")
histos["h_px_bkg"].Draw("hist same")
histos["h_px_rec"].Draw("hist same")
histos["h_px_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.cd(2)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_py_sig"],histos["h_py_bkg"])
hmin,hmax = minmax(histos["h_py_sig"],histos["h_py_rec"])
hmin,hmax = minmax(histos["h_py_rec"],histos["h_py_sed"])
hmin,hmax = minmax(histos["h_py_sed"],histos["h_py_sig"])

hmin,hmax = minmax(histos["h_py_sig_Eside"],histos["h_py_bkg_Eside"])
hmin,hmax = minmax(histos["h_py_sig_Eside"],histos["h_py_rec_Eside"])
hmin,hmax = minmax(histos["h_py_rec_Eside"],histos["h_py_sed_Eside"])
hmin,hmax = minmax(histos["h_py_sed_Eside"],histos["h_py_sig_Eside"])

hmin,hmax = minmax(histos["h_py_sig_Pside"],histos["h_py_bkg_Pside"])
hmin,hmax = minmax(histos["h_py_sig_Pside"],histos["h_py_rec_Pside"])
hmin,hmax = minmax(histos["h_py_rec_Pside"],histos["h_py_sed_Pside"])
hmin,hmax = minmax(histos["h_py_sed_Pside"],histos["h_py_sig_Pside"])

histos["h_py_sig"].SetLineColor(ROOT.kRed)
histos["h_py_bkg"].SetLineColor(ROOT.kBlue)
histos["h_py_rec"].SetLineColor(ROOT.kBlack)
histos["h_py_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_py_sig"].Draw("hist")
histos["h_py_bkg"].Draw("hist same")
histos["h_py_rec"].Draw("hist same")
histos["h_py_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.cd(3)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_pz_sig"],histos["h_pz_bkg"])
hmin,hmax = minmax(histos["h_pz_sig"],histos["h_pz_rec"])
hmin,hmax = minmax(histos["h_pz_rec"],histos["h_pz_sed"])
hmin,hmax = minmax(histos["h_pz_sed"],histos["h_pz_sig"])

hmin,hmax = minmax(histos["h_pz_sig_Eside"],histos["h_pz_bkg_Eside"])
hmin,hmax = minmax(histos["h_pz_sig_Eside"],histos["h_pz_rec_Eside"])
hmin,hmax = minmax(histos["h_pz_rec_Eside"],histos["h_pz_sed_Eside"])
hmin,hmax = minmax(histos["h_pz_sed_Eside"],histos["h_pz_sig_Eside"])

hmin,hmax = minmax(histos["h_pz_sig_Pside"],histos["h_pz_bkg_Pside"])
hmin,hmax = minmax(histos["h_pz_sig_Pside"],histos["h_pz_rec_Pside"])
hmin,hmax = minmax(histos["h_pz_rec_Pside"],histos["h_pz_sed_Pside"])
hmin,hmax = minmax(histos["h_pz_sed_Pside"],histos["h_pz_sig_Pside"])

histos["h_pz_sig"].SetLineColor(ROOT.kRed)
histos["h_pz_bkg"].SetLineColor(ROOT.kBlue)
histos["h_pz_rec"].SetLineColor(ROOT.kBlack)
histos["h_pz_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_pz_sig"].Draw("hist")
histos["h_pz_bkg"].Draw("hist same")
histos["h_pz_rec"].Draw("hist same")
histos["h_pz_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxyz.SaveAs(fn+"Momentum.pdf")
cnv_pxyz.SaveAs(allpdf)


cnv_pxyz = TCanvas("cnv_pxyz_norec","",1500,500)
cnv_pxyz.Divide(3,1)
cnv_pxyz.cd(1)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_px_sig"],histos["h_px_bkg"])
hmin,hmax = minmax(histos["h_px_sig_Eside"],histos["h_px_bkg_Eside"])
hmin,hmax = minmax(histos["h_px_sig_Pside"],histos["h_px_bkg_Pside"])
histos["h_px_sig"].SetLineColor(ROOT.kRed)
histos["h_px_bkg"].SetLineColor(ROOT.kBlue)
histos["h_px_sig"].Draw("hist")
histos["h_px_bkg"].Draw("hist same")
leg_bkg_sig.Draw("same")
cnv_pxyz.cd(2)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_py_sig"],histos["h_py_bkg"])
hmin,hmax = minmax(histos["h_py_sig_Eside"],histos["h_py_bkg_Eside"])
hmin,hmax = minmax(histos["h_py_sig_Pside"],histos["h_py_bkg_Pside"])
histos["h_py_sig"].SetLineColor(ROOT.kRed)
histos["h_py_bkg"].SetLineColor(ROOT.kBlue)
histos["h_py_sig"].Draw("hist")
histos["h_py_bkg"].Draw("hist same")
leg_bkg_sig.Draw("same")
cnv_pxyz.cd(3)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_pz_sig"],histos["h_pz_bkg"])
hmin,hmax = minmax(histos["h_pz_sig_Eside"],histos["h_pz_bkg_Eside"])
hmin,hmax = minmax(histos["h_pz_sig_Pside"],histos["h_pz_bkg_Pside"])
histos["h_pz_sig"].SetLineColor(ROOT.kRed)
histos["h_pz_bkg"].SetLineColor(ROOT.kBlue)
histos["h_pz_sig"].Draw("hist")
histos["h_pz_bkg"].Draw("hist same")
leg_bkg_sig.Draw("same")
cnv_pxyz.SaveAs(fn+"Momentum_norec.pdf")
cnv_pxyz.SaveAs(allpdf)


cnv_pxy_zoom = TCanvas("cnv_pxyz_zoom","",1000,500)
cnv_pxy_zoom.Divide(2,1)
cnv_pxy_zoom.cd(1)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_px_zoom_sig"],histos["h_px_zoom_bkg"])
hmin,hmax = minmax(histos["h_px_zoom_sig"],histos["h_px_zoom_rec"])
hmin,hmax = minmax(histos["h_px_zoom_rec"],histos["h_px_zoom_sed"])
hmin,hmax = minmax(histos["h_px_zoom_sed"],histos["h_px_zoom_sig"])

hmin,hmax = minmax(histos["h_px_zoom_sig_Eside"],histos["h_px_zoom_bkg_Eside"])
hmin,hmax = minmax(histos["h_px_zoom_sig_Eside"],histos["h_px_zoom_rec_Eside"])
hmin,hmax = minmax(histos["h_px_zoom_rec_Eside"],histos["h_px_zoom_sed_Eside"])
hmin,hmax = minmax(histos["h_px_zoom_sed_Eside"],histos["h_px_zoom_sig_Eside"])

hmin,hmax = minmax(histos["h_px_zoom_sig_Pside"],histos["h_px_zoom_bkg_Pside"])
hmin,hmax = minmax(histos["h_px_zoom_sig_Pside"],histos["h_px_zoom_rec_Pside"])
hmin,hmax = minmax(histos["h_px_zoom_rec_Pside"],histos["h_px_zoom_sed_Pside"])
hmin,hmax = minmax(histos["h_px_zoom_sed_Pside"],histos["h_px_zoom_sig_Pside"])

histos["h_px_zoom_sig"].SetLineColor(ROOT.kRed)
histos["h_px_zoom_bkg"].SetLineColor(ROOT.kBlue)
histos["h_px_zoom_rec"].SetLineColor(ROOT.kBlack)
histos["h_px_zoom_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_px_zoom_sig"].Draw("hist")
histos["h_px_zoom_bkg"].Draw("hist same")
histos["h_px_zoom_rec"].Draw("hist same")
histos["h_px_zoom_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxy_zoom.cd(2)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_py_zoom_sig"],histos["h_py_zoom_bkg"])
hmin,hmax = minmax(histos["h_py_zoom_sig"],histos["h_py_zoom_rec"])
hmin,hmax = minmax(histos["h_py_zoom_rec"],histos["h_py_zoom_sed"])
hmin,hmax = minmax(histos["h_py_zoom_sed"],histos["h_py_zoom_sig"])

hmin,hmax = minmax(histos["h_py_zoom_sig_Eside"],histos["h_py_zoom_bkg_Eside"])
hmin,hmax = minmax(histos["h_py_zoom_sig_Eside"],histos["h_py_zoom_rec_Eside"])
hmin,hmax = minmax(histos["h_py_zoom_rec_Eside"],histos["h_py_zoom_sed_Eside"])
hmin,hmax = minmax(histos["h_py_zoom_sed_Eside"],histos["h_py_zoom_sig_Eside"])

hmin,hmax = minmax(histos["h_py_zoom_sig_Pside"],histos["h_py_zoom_bkg_Pside"])
hmin,hmax = minmax(histos["h_py_zoom_sig_Pside"],histos["h_py_zoom_rec_Pside"])
hmin,hmax = minmax(histos["h_py_zoom_rec_Pside"],histos["h_py_zoom_sed_Pside"])
hmin,hmax = minmax(histos["h_py_zoom_sed_Pside"],histos["h_py_zoom_sig_Pside"])

histos["h_py_zoom_sig"].SetLineColor(ROOT.kRed)
histos["h_py_zoom_bkg"].SetLineColor(ROOT.kBlue)
histos["h_py_zoom_rec"].SetLineColor(ROOT.kBlack)
histos["h_py_zoom_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_py_zoom_sig"].Draw("hist")
histos["h_py_zoom_bkg"].Draw("hist same")
histos["h_py_zoom_rec"].Draw("hist same")
histos["h_py_zoom_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_pxy_zoom.SaveAs(fn+"pxy_zoom.pdf")
cnv_pxy_zoom.SaveAs(allpdf)


cnv_ntrks = TCanvas("cnv_ntrks","",500,500)
cnv_ntrks.cd()
hmin,hmax = minmax(histos["h_ntrks_sig"],histos["h_ntrks_bkg"])
hmin,hmax = minmax(histos["h_ntrks_sig"],histos["h_ntrks_rec"])
hmin,hmax = minmax(histos["h_ntrks_rec"],histos["h_ntrks_sed"])
hmin,hmax = minmax(histos["h_ntrks_sed"],histos["h_ntrks_sig"])

hmin,hmax = minmax(histos["h_ntrks_sig_Eside"],histos["h_ntrks_bkg_Eside"])
hmin,hmax = minmax(histos["h_ntrks_sig_Eside"],histos["h_ntrks_rec_Eside"])
hmin,hmax = minmax(histos["h_ntrks_rec_Eside"],histos["h_ntrks_sed_Eside"])
hmin,hmax = minmax(histos["h_ntrks_sed_Eside"],histos["h_ntrks_sig_Eside"])

hmin,hmax = minmax(histos["h_ntrks_sig_Pside"],histos["h_ntrks_bkg_Pside"])
hmin,hmax = minmax(histos["h_ntrks_sig_Pside"],histos["h_ntrks_rec_Pside"])
hmin,hmax = minmax(histos["h_ntrks_rec_Pside"],histos["h_ntrks_sed_Pside"])
hmin,hmax = minmax(histos["h_ntrks_sed_Pside"],histos["h_ntrks_sig_Pside"])

histos["h_ntrks_sig"].SetLineColor(ROOT.kRed)
histos["h_ntrks_bkg"].SetLineColor(ROOT.kBlue)
histos["h_ntrks_rec"].SetLineColor(ROOT.kBlack)
histos["h_ntrks_sed"].SetLineColor(ROOT.kGreen+2)
histos["h_ntrks_sig"].Draw("hist")
histos["h_ntrks_bkg"].Draw("hist same")
histos["h_ntrks_rec"].Draw("hist same")
histos["h_ntrks_sed"].Draw("hist same")
leg_rec_sig.Draw("same")
cnv_ntrks.SaveAs(fn+"ntrks.pdf")
cnv_ntrks.SaveAs(allpdf)


cnv_ntrks_norec = TCanvas("cnv_ntrks_norec","",500,500)
cnv_ntrks_norec.cd()
hmin,hmax = minmax(histos["h_ntrks_sig"],histos["h_ntrks_bkg"])
hmin,hmax = minmax(histos["h_ntrks_sig_Eside"],histos["h_ntrks_bkg_Eside"])
hmin,hmax = minmax(histos["h_ntrks_sig_Pside"],histos["h_ntrks_bkg_Pside"])
histos["h_ntrks_sig"].SetLineColor(ROOT.kRed)
histos["h_ntrks_bkg"].SetLineColor(ROOT.kBlue)
histos["h_ntrks_sig"].Draw("hist")
histos["h_ntrks_bkg"].Draw("hist same")
leg_bkg_sig.Draw("same")
cnv_ntrks_norec.SaveAs(fn+"ntrks_norec.pdf")
cnv_ntrks_norec.SaveAs(allpdf)


histos["h_E_tru_eff_sed"].Divide(histos["h_E_tru_sed_mat"],histos["h_E_tru_sig"]);
histos["h_E_tru_eff_sed_Eside"].Divide(histos["h_E_tru_sed_mat_Eside"],histos["h_E_tru_sig_Eside"]);
histos["h_E_tru_eff_sed_Pside"].Divide(histos["h_E_tru_sed_mat_Pside"],histos["h_E_tru_sig_Pside"]);
cnv_eff_seeds = TCanvas("cnv_eff_seeds","",500,500)
cnv_eff_seeds.cd()
histos["h_E_tru_eff_sed"].SetMinimum(0)
histos["h_E_tru_eff_sed"].SetMaximum(1.05)
histos["h_E_tru_eff_sed"].SetLineColor(ROOT.kBlack)
histos["h_E_tru_eff_sed"].SetMarkerColor(ROOT.kBlack)
histos["h_E_tru_eff_sed"].SetMarkerStyle(20)
histos["h_E_tru_eff_sed"].Draw("ep")
cnv_eff_seeds.SaveAs(fn+"eff_seeds.pdf")
cnv_eff_seeds.SaveAs(allpdf)

histos["h_E_tru_eff"].Divide(histos["h_E_tru_rec_mat"],histos["h_E_tru_sig"]);
histos["h_E_tru_eff_Eside"].Divide(histos["h_E_tru_rec_mat_Eside"],histos["h_E_tru_sig_Eside"]);
histos["h_E_tru_eff_Pside"].Divide(histos["h_E_tru_rec_mat_Pside"],histos["h_E_tru_sig_Pside"]);
histos["h_E_tru_eff_acc"].Divide(histos["h_E_tru_rec_mat_acc"],histos["h_E_tru_sig_acc"]);
histos["h_E_tru_eff_acc_Eside"].Divide(histos["h_E_tru_rec_mat_acc_Eside"],histos["h_E_tru_sig_acc_Eside"]);
histos["h_E_tru_eff_acc_Pside"].Divide(histos["h_E_tru_rec_mat_acc_Pside"],histos["h_E_tru_sig_acc_Pside"]);
cnv_eff = TCanvas("cnv_eff","",1000,500)
cnv_eff.Divide(2,1)
cnv_eff.cd(1)
histos["h_E_tru_eff"].SetMinimum(0)
histos["h_E_tru_eff"].SetMaximum(1.05)
histos["h_E_tru_eff"].SetLineColor(ROOT.kBlack)
histos["h_E_tru_eff"].SetMarkerColor(ROOT.kBlack)
histos["h_E_tru_eff"].SetMarkerStyle(20)
histos["h_E_tru_eff"].Draw("ep")
cnv_eff.cd(2)
histos["h_E_tru_eff_acc"].SetMinimum(0)
histos["h_E_tru_eff_acc"].SetMaximum(1.05)
histos["h_E_tru_eff_acc"].SetLineColor(ROOT.kBlack)
histos["h_E_tru_eff_acc"].SetMarkerColor(ROOT.kBlack)
histos["h_E_tru_eff_acc"].SetMarkerStyle(20)
histos["h_E_tru_eff_acc"].Draw("ep")
cnv_eff.SaveAs(fn+"eff.pdf")
cnv_eff.SaveAs(allpdf)

histos["h_E_tru_eff_sel"].Divide(histos["h_E_tru_sel_mat"],histos["h_E_tru_sig"]);
histos["h_E_tru_eff_sel_Eside"].Divide(histos["h_E_tru_sel_mat_Eside"],histos["h_E_tru_sig_Eside"]);
histos["h_E_tru_eff_sel_Pside"].Divide(histos["h_E_tru_sel_mat_Pside"],histos["h_E_tru_sig_Pside"]);
histos["h_E_tru_eff_acc_sel"].Divide(histos["h_E_tru_sel_mat_acc"],histos["h_E_tru_sig_acc"]);
histos["h_E_tru_eff_acc_sel_Eside"].Divide(histos["h_E_tru_sel_mat_acc_Eside"],histos["h_E_tru_sig_acc_Eside"]);
histos["h_E_tru_eff_acc_sel_Pside"].Divide(histos["h_E_tru_sel_mat_acc_Pside"],histos["h_E_tru_sig_acc_Pside"]);
cnv_eff_sel = TCanvas("cnv_eff_sel","",1000,500)
cnv_eff_sel.Divide(2,1)
cnv_eff_sel.cd(1)
histos["h_E_tru_eff_sel"].SetMinimum(0)
histos["h_E_tru_eff_sel"].SetMaximum(1.05)
histos["h_E_tru_eff_sel"].SetLineColor(ROOT.kBlack)
histos["h_E_tru_eff_sel"].SetMarkerColor(ROOT.kBlack)
histos["h_E_tru_eff_sel"].SetMarkerStyle(20)
histos["h_E_tru_eff_sel"].Draw("ep")
cnv_eff_sel.cd(2)
histos["h_E_tru_eff_acc_sel"].SetMinimum(0)
histos["h_E_tru_eff_acc_sel"].SetMaximum(1.05)
histos["h_E_tru_eff_acc_sel"].SetLineColor(ROOT.kBlack)
histos["h_E_tru_eff_acc_sel"].SetMarkerColor(ROOT.kBlack)
histos["h_E_tru_eff_acc_sel"].SetMarkerStyle(20)
histos["h_E_tru_eff_acc_sel"].Draw("ep")
cnv_eff_sel.SaveAs(fn+"eff_sel.pdf")
cnv_eff_sel.SaveAs(allpdf)


cnv_res_xy = TCanvas("cnv_res_xy","",1000,1500)
cnv_res_xy.Divide(2,3)
cnv_res_xy.cd(1)
gaus_res_x_truvscls = TF1("gaus_res_x_truvscls","gaus",-0.00005,+0.00005)
res = histos["h_res_x_truvscls"].Fit(gaus_res_x_truvscls,"EMRS")
chi2dof = gaus_res_x_truvscls.GetChisquare()/gaus_res_x_truvscls.GetNDF() if(gaus_res_x_truvscls.GetNDF()>0) else -1
print("Res(x tru:cls) chi2/Ndof=",chi2dof)
histos["h_res_x_truvscls"].Draw("hist")
gaus_res_x_truvscls.Draw("same")
cnv_res_xy.cd(2);
gaus_res_y_truvscls = TF1("gaus_res_y_truvscls","gaus",-0.01,+0.01)
res = histos["h_res_y_truvscls"].Fit(gaus_res_y_truvscls,"EMRS")
chi2dof = gaus_res_y_truvscls.GetChisquare()/gaus_res_y_truvscls.GetNDF() if(gaus_res_y_truvscls.GetNDF()>0) else -1
print("Res(y tru:cls) chi2/Ndof=",chi2dof)
histos["h_res_y_truvscls"].Draw("hist")
gaus_res_y_truvscls.Draw("same")
cnv_res_xy.cd(3);
gaus_res_x_recvscls = TF1("gaus_res_x_recvscls","gaus",-0.00001,+0.00001)
res = histos["h_res_x_recvscls"].Fit(gaus_res_x_recvscls,"EMRS")
chi2dof = gaus_res_x_recvscls.GetChisquare()/gaus_res_x_recvscls.GetNDF() if(gaus_res_x_recvscls.GetNDF()>0) else -1
print("Res(x rec:cls) chi2/Ndof=",chi2dof)
histos["h_res_x_recvscls"].Draw("hist")
gaus_res_x_recvscls.Draw("same")
cnv_res_xy.cd(4)
gaus_res_y_recvscls = TF1("gaus_res_y_recvscls","gaus",-0.0035,+0.0035)
res = histos["h_res_y_recvscls"].Fit(gaus_res_y_recvscls,"EMRS")
chi2dof = gaus_res_y_recvscls.GetChisquare()/gaus_res_y_recvscls.GetNDF() if(gaus_res_y_recvscls.GetNDF()>0) else -1
print("Res(y rec:cls) chi2/Ndof=",chi2dof)
histos["h_res_y_recvscls"].Draw("hist")
gaus_res_y_recvscls.Draw("same")
cnv_res_xy.cd(5)
gaus_res_x_recvstru = TF1("gaus_res_x_recvstru","gaus",-0.00006,+0.00006)
res = histos["h_res_x_recvstru"].Fit(gaus_res_x_recvstru,"EMRS")
chi2dof = gaus_res_x_recvstru.GetChisquare()/gaus_res_x_recvstru.GetNDF() if(gaus_res_x_recvstru.GetNDF()>0) else -1
print("Res(x rec:tru) chi2/Ndof=",chi2dof)
histos["h_res_x_recvstru"].Draw("hist")
gaus_res_x_recvstru.Draw("same")
cnv_res_xy.cd(6)
gaus_res_y_recvstru = TF1("gaus_res_y_recvstru","gaus",-0.01,+0.01)
res = histos["h_res_y_recvstru"].Fit(gaus_res_y_recvstru,"EMRS")
chi2dof = gaus_res_y_recvstru.GetChisquare()/gaus_res_y_recvstru.GetNDF() if(gaus_res_y_recvstru.GetNDF()>0) else -1
print("Res(y rec:tru) chi2/Ndof=",chi2dof)
histos["h_res_y_recvstru"].Draw("hist")
gaus_res_y_recvstru.Draw("same")
cnv_res_xy.SaveAs(fn+"res_xy.pdf")
cnv_res_xy.SaveAs(allpdf)


cnv_chi2 = TCanvas("cnv_chi2","",1500,500)
cnv_chi2.Divide(3,1)
hmin,hmax = minmax(histos["h_chi2dof_rec"],histos["h_chi2dof_rec_sig_mat"])
hmin,hmax = minmax(histos["h_chi2dof_rec_Eside"],histos["h_chi2dof_rec_sig_mat_Eside"])
hmin,hmax = minmax(histos["h_chi2dof_rec_Pside"],histos["h_chi2dof_rec_sig_mat_Pside"])
hmin,hmax = minmax(histos["h_chi2dof_rec_sig_mat"],histos["h_chi2dof_rec_sig_non"])
hmin,hmax = minmax(histos["h_chi2dof_rec_sig_mat_Eside"],histos["h_chi2dof_rec_sig_non_Eside"])
hmin,hmax = minmax(histos["h_chi2dof_rec_sig_mat_Pside"],histos["h_chi2dof_rec_sig_non_Pside"])
histos["h_chi2dof_rec"].SetMinimum(0.5)
histos["h_chi2dof_rec_sig_mat"].SetMinimum(0.5)
histos["h_chi2dof_rec_sig_non"].SetMinimum(0.5)
cnv_chi2.cd(1)
ROOT.gPad.SetLogy()
histos["h_chi2dof_rec"].Draw("hist")
cnv_chi2.cd(2)
ROOT.gPad.SetLogy()
histos["h_chi2dof_rec_sig_mat"].Draw("hist")
cnv_chi2.cd(3)
ROOT.gPad.SetLogy()
histos["h_chi2dof_rec_sig_non"].Draw("hist")
cnv_chi2.SaveAs(fn+"chi2dof.pdf")
cnv_chi2.SaveAs(allpdf)

cnv_SnpSig = TCanvas("cnv_SnpSig","",1500,500)
cnv_SnpSig.Divide(3,1)
hmin,hmax = minmax(histos["h_SnpSig_rec"],histos["h_SnpSig_rec_sig_mat"])
hmin,hmax = minmax(histos["h_SnpSig_rec_Eside"],histos["h_SnpSig_rec_sig_mat_Eside"])
hmin,hmax = minmax(histos["h_SnpSig_rec_Pside"],histos["h_SnpSig_rec_sig_mat_Pside"])
hmin,hmax = minmax(histos["h_SnpSig_rec_sig_mat"],histos["h_SnpSig_rec_sig_non"])
hmin,hmax = minmax(histos["h_SnpSig_rec_sig_mat_Eside"],histos["h_SnpSig_rec_sig_non_Eside"])
hmin,hmax = minmax(histos["h_SnpSig_rec_sig_mat_Pside"],histos["h_SnpSig_rec_sig_non_Pside"])
histos["h_SnpSig_rec"].SetMinimum(0.5)
histos["h_SnpSig_rec_sig_mat"].SetMinimum(0.5)
histos["h_SnpSig_rec_sig_non"].SetMinimum(0.5)
cnv_SnpSig.cd(1)
ROOT.gPad.SetLogy()
histos["h_SnpSig_rec"].Draw("hist")
cnv_SnpSig.cd(2)
ROOT.gPad.SetLogy()
histos["h_SnpSig_rec_sig_mat"].Draw("hist")
cnv_SnpSig.cd(3)
ROOT.gPad.SetLogy()
histos["h_SnpSig_rec_sig_non"].Draw("hist")
cnv_SnpSig.SaveAs(fn+"SnpSig.pdf")
cnv_SnpSig.SaveAs(allpdf)

cnv_TglSig = TCanvas("cnv_TglSig","",1500,500)
cnv_TglSig.Divide(3,1)
hmin,hmax = minmax(histos["h_TglSig_rec"],histos["h_TglSig_rec_sig_mat"])
hmin,hmax = minmax(histos["h_TglSig_rec_Eside"],histos["h_TglSig_rec_sig_mat_Eside"])
hmin,hmax = minmax(histos["h_TglSig_rec_Pside"],histos["h_TglSig_rec_sig_mat_Pside"])
hmin,hmax = minmax(histos["h_TglSig_rec_sig_mat"],histos["h_TglSig_rec_sig_non"])
hmin,hmax = minmax(histos["h_TglSig_rec_sig_mat_Eside"],histos["h_TglSig_rec_sig_non_Eside"])
hmin,hmax = minmax(histos["h_TglSig_rec_sig_mat_Pside"],histos["h_TglSig_rec_sig_non_Pside"])
histos["h_TglSig_rec"].SetMinimum(0.5)
histos["h_TglSig_rec_sig_mat"].SetMinimum(0.5)
histos["h_TglSig_rec_sig_non"].SetMinimum(0.5)
cnv_TglSig.cd(1)
ROOT.gPad.SetLogy()
histos["h_TglSig_rec"].Draw("hist")
cnv_TglSig.cd(2)
ROOT.gPad.SetLogy()
histos["h_TglSig_rec_sig_mat"].Draw("hist")
cnv_TglSig.cd(3)
ROOT.gPad.SetLogy()
histos["h_TglSig_rec_sig_non"].Draw("hist")
cnv_TglSig.SaveAs(fn+"TglSig.pdf")
cnv_TglSig.SaveAs(allpdf)

cnv_vtxSig = TCanvas("cnv_vtxSig","",1500,1000)
cnv_vtxSig.Divide(3,2)
cnv_vtxSig.cd(1)
ROOT.gPad.SetLogy()
hmin,hmax = minmax(histos["h_xVtxSig_rec"],histos["h_xVtxSig_rec_sig_mat"])
hmin,hmax = minmax(histos["h_xVtxSig_rec_Eside"],histos["h_xVtxSig_rec_sig_mat_Eside"])
hmin,hmax = minmax(histos["h_xVtxSig_rec_Pside"],histos["h_xVtxSig_rec_sig_mat_Pside"])
hmin,hmax = minmax(histos["h_xVtxSig_rec_sig_mat"],histos["h_xVtxSig_rec_sig_non"])
hmin,hmax = minmax(histos["h_xVtxSig_rec_sig_mat_Eside"],histos["h_xVtxSig_rec_sig_non_Eside"])
hmin,hmax = minmax(histos["h_xVtxSig_rec_sig_mat_Pside"],histos["h_xVtxSig_rec_sig_non_Pside"])
histos["h_xVtxSig_rec"].SetMinimum(0.5)
histos["h_xVtxSig_rec_sig_mat"].SetMinimum(0.5)
histos["h_xVtxSig_rec_sig_non"].SetMinimum(0.5)
histos["h_xVtxSig_rec"].Draw("hist")
cnv_vtxSig.cd(2)
ROOT.gPad.SetLogy()
histos["h_xVtxSig_rec_sig_mat"].Draw("hist")
cnv_vtxSig.cd(3)
ROOT.gPad.SetLogy()
histos["h_xVtxSig_rec_sig_non"].Draw("hist")
hmin,hmax = minmax(histos["h_yVtxSig_rec"],histos["h_yVtxSig_rec_sig_mat"])
hmin,hmax = minmax(histos["h_yVtxSig_rec_Eside"],histos["h_yVtxSig_rec_sig_mat_Eside"])
hmin,hmax = minmax(histos["h_yVtxSig_rec_Pside"],histos["h_yVtxSig_rec_sig_mat_Pside"])
hmin,hmax = minmax(histos["h_yVtxSig_rec_sig_mat"],histos["h_yVtxSig_rec_sig_non"])
hmin,hmax = minmax(histos["h_yVtxSig_rec_sig_mat_Eside"],histos["h_yVtxSig_rec_sig_non_Eside"])
hmin,hmax = minmax(histos["h_yVtxSig_rec_sig_mat_Pside"],histos["h_yVtxSig_rec_sig_non_Pside"])
histos["h_yVtxSig_rec"].SetMinimum(0.5)
histos["h_yVtxSig_rec_sig_mat"].SetMinimum(0.5)
histos["h_yVtxSig_rec_sig_non"].SetMinimum(0.5)
cnv_vtxSig.cd(4)
ROOT.gPad.SetLogy()
histos["h_yVtxSig_rec"].Draw("hist")
cnv_vtxSig.cd(5)
ROOT.gPad.SetLogy()
histos["h_yVtxSig_rec_sig_mat"].Draw("hist")
cnv_vtxSig.cd(6)
ROOT.gPad.SetLogy()
histos["h_yVtxSig_rec_sig_non"].Draw("hist")
cnv_vtxSig.SaveAs(fn+"vtxSig.pdf")
cnv_vtxSig.SaveAs(allpdf)


cnv_rat_E = TCanvas("cnv_rat_E","",500,500)
cnv_rat_E.cd()
hmin,hmax = minmax(histos["h_rat_E"],histos["h_rat_E"])
hmin,hmax = minmax(histos["h_rat_E_Eside"],histos["h_rat_E_Eside"])
hmin,hmax = minmax(histos["h_rat_E_Pside"],histos["h_rat_E_Pside"])
gaus_ratE = TF1("gaus_ratE","gaus",0.995,+1.004)
rat = histos["h_rat_E"].Fit(gaus_ratE,"EMRS")
chi2dof = gaus_ratE.GetChisquare()/gaus_ratE.GetNDF() if(gaus_ratE.GetNDF()>0) else -1
print("Res(E) chi2/Ndof=",chi2dof)
histos["h_rat_E"].Draw("hist")
gaus_ratE.Draw("same")
histos["h_rat_sed_E"].SetLineColor(ROOT.kGreen+2)
histos["h_rat_sed_E"].Draw("hist same")
leg_rec_sed = TLegend(0.60,0.73,0.87,0.87)
leg_rec_sed.SetFillStyle(4000) # will be transparent
leg_rec_sed.SetFillColor(0)
leg_rec_sed.SetTextFont(42)
leg_rec_sed.SetBorderSize(0)
leg_rec_sed.AddEntry(histos["h_rat_sed_E"],"Seeding","l")
leg_rec_sed.AddEntry(histos["h_rat_E"],"Reconstruction","l")
leg_rec_sed.Draw("same")
cnv_rat_E.SaveAs(fn+"rat_E.pdf")
cnv_rat_E.SaveAs(allpdf)


cnv_res_E = TCanvas("cnv_res_E","",500,500)
cnv_res_E.cd()
hmin,hmax = minmax(histos["h_res_E"],histos["h_res_E"])
hmin,hmax = minmax(histos["h_res_E_Eside"],histos["h_res_E_Eside"])
hmin,hmax = minmax(histos["h_res_E_Pside"],histos["h_res_E_Pside"])
gaus_resE = TF1("gaus_resE","gaus",-0.004,+0.004)
res = histos["h_res_E"].Fit(gaus_resE,"EMRS")
chi2dof = gaus_resE.GetChisquare()/gaus_resE.GetNDF() if(gaus_resE.GetNDF()>0) else -1
print("Res(E) chi2/Ndof=",chi2dof)
histos["h_res_E"].Draw("hist")
gaus_resE.Draw("same")
histos["h_res_sed_E"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_E"].Draw("hist same")
leg_rec_sed = TLegend(0.60,0.73,0.87,0.87)
leg_rec_sed.SetFillStyle(4000) # will be transparent
leg_rec_sed.SetFillColor(0)
leg_rec_sed.SetTextFont(42)
leg_rec_sed.SetBorderSize(0)
leg_rec_sed.AddEntry(histos["h_res_sed_E"],"Seeding","l")
leg_rec_sed.AddEntry(histos["h_res_E"],"Reconstruction","l")
leg_rec_sed.Draw("same")
cnv_res_E.SaveAs(fn+"res_E.pdf")
cnv_res_E.SaveAs(allpdf)

cnv_res_pz = TCanvas("cnv_res_pz","",500,500)
cnv_res_pz.cd()
hmin,hmax = minmax(histos["h_res_pz"],histos["h_res_pz"])
hmin,hmax = minmax(histos["h_res_pz_Eside"],histos["h_res_pz_Eside"])
hmin,hmax = minmax(histos["h_res_pz_Pside"],histos["h_res_pz_Pside"])
gaus_respz = TF1("gaus_respz","gaus",-0.003,+0.003)
res = histos["h_res_pz"].Fit(gaus_respz,"EMRS")
chi2dof = gaus_respz.GetChisquare()/gaus_respz.GetNDF() if(gaus_respz.GetNDF()>0) else -1
print("Res(pz) chi2/Ndof=",chi2dof)
histos["h_res_pz"].Draw("hist")
gaus_respz.Draw("same")
histos["h_res_sed_pz"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_pz"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_pz.SaveAs(fn+"res_pz.pdf")
cnv_res_pz.SaveAs(allpdf)

cnv_res_px = TCanvas("cnv_res_px","",500,500)
cnv_res_px.cd()
hmin,hmax = minmax(histos["h_res_px"],histos["h_res_px"])
hmin,hmax = minmax(histos["h_res_px_Eside"],histos["h_res_px_Eside"])
hmin,hmax = minmax(histos["h_res_px_Pside"],histos["h_res_px_Pside"])
gaus_respx = TF1("gaus_respx","gaus",-0.75,+0.75)
res = histos["h_res_px"].Fit(gaus_respx,"EMRS")
chi2dof = gaus_respx.GetChisquare()/gaus_respx.GetNDF() if(gaus_respx.GetNDF()>0) else -1
print("Res(px) chi2/Ndof=",chi2dof)
histos["h_res_px"].Draw("hist")
gaus_respx.Draw("same")
histos["h_res_sed_px"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_px"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_px.SaveAs(fn+"res_px.pdf")
cnv_res_px.SaveAs(allpdf)

cnv_res_py = TCanvas("cnv_res_py","",500,500)
cnv_res_py.cd()
hmin,hmax = minmax(histos["h_res_py"],histos["h_res_sed_py"])
hmin,hmax = minmax(histos["h_res_py_Eside"],histos["h_res_sed_py_Eside"])
hmin,hmax = minmax(histos["h_res_py_Pside"],histos["h_res_sed_py_Pside"])
gaus_respy = TF1("gaus_respy","gaus",-0.012,+0.012)
res = histos["h_res_py"].Fit(gaus_respy,"EMRS")
chi2dof = gaus_respy.GetChisquare()/gaus_respy.GetNDF() if(gaus_respy.GetNDF()>0) else -1
print("Res(py) chi2/Ndof=",chi2dof)
histos["h_res_py"].Draw("hist")
gaus_respy.Draw("same")
histos["h_res_sed_py"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_py"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_py.SaveAs(fn+"res_py.pdf")
cnv_res_py.SaveAs(allpdf)

cnv_res_py_zoomout = TCanvas("cnv_res_py_zoomout","",500,500)
cnv_res_py_zoomout.cd()
hmin,hmax = minmax(histos["h_res_py_zoomout"],histos["h_res_sed_py_zoomout"])
hmin,hmax = minmax(histos["h_res_py_zoomout_Eside"],histos["h_res_sed_py_zoomout_Eside"])
hmin,hmax = minmax(histos["h_res_py_zoomout_Pside"],histos["h_res_sed_py_zoomout_Pside"])
histos["h_res_py_zoomout"].Draw("hist")
histos["h_res_sed_py_zoomout"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_py_zoomout"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_py_zoomout.SaveAs(fn+"res_py_zoomout.pdf")
cnv_res_py_zoomout.SaveAs(allpdf)

cnv_res_pT = TCanvas("cnv_res_pT","",500,500)
cnv_res_pT.cd()
hmin,hmax = minmax(histos["h_res_pT"],histos["h_res_pT"])
hmin,hmax = minmax(histos["h_res_pT_Eside"],histos["h_res_pT_Eside"])
hmin,hmax = minmax(histos["h_res_pT_Pside"],histos["h_res_pT_Pside"])
histos["h_res_pT"].Draw("hist")
histos["h_res_sed_pT"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_pT"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_pT.SaveAs(fn+"res_pT.pdf")
cnv_res_pT.SaveAs(allpdf)

cnv_res_theta = TCanvas("cnv_res_theta","",500,500)
cnv_res_theta.cd()
hmin,hmax = minmax(histos["h_res_theta"],histos["h_res_theta"])
hmin,hmax = minmax(histos["h_res_theta_Eside"],histos["h_res_theta_Eside"])
hmin,hmax = minmax(histos["h_res_theta_Pside"],histos["h_res_theta_Pside"])
histos["h_res_theta"].Draw("hist")
histos["h_res_sed_theta"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_theta"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_theta.SaveAs(fn+"res_theta.pdf")
cnv_res_theta.SaveAs(allpdf)

cnv_res_phi = TCanvas("cnv_res_phi","",500,500)
cnv_res_phi.cd()
hmin,hmax = minmax(histos["h_res_phi"],histos["h_res_phi"])
hmin,hmax = minmax(histos["h_res_phi_Eside"],histos["h_res_phi_Eside"])
hmin,hmax = minmax(histos["h_res_phi_Pside"],histos["h_res_phi_Pside"])
histos["h_res_phi"].Draw("hist")
histos["h_res_sed_phi"].SetLineColor(ROOT.kGreen+2)
histos["h_res_sed_phi"].Draw("hist same")
leg_rec_sed.Draw("same")
cnv_res_phi.SaveAs(fn+"res_phi.pdf")
cnv_res_phi.SaveAs(allpdf)

cnv_rat_E_vs_E = TCanvas("cnv_rat_E_vs_E","",500,500)
cnv_rat_E_vs_E.cd()
histos["h_rat_E_vs_Etru_recvstru"].Draw("col")
cnv_rat_E_vs_E.SaveAs(fn+"rat_E_vs_E.pdf")
cnv_rat_E_vs_E.SaveAs(allpdf)

cnv_res_E_vs_E = TCanvas("cnv_res_E_vs_E","",500,500)
cnv_res_E_vs_E.cd()
histos["h_res_E_vs_Etru_recvstru"].Draw("col")
cnv_res_E_vs_E.SaveAs(fn+"res_E_vs_E.pdf")
cnv_res_E_vs_E.SaveAs(allpdf)

cnv_dE_vs_E = TCanvas("cnv_dE_vs_E","",500,500)
cnv_dE_vs_E.cd()
histos["h_dE_vs_Etru_recvstru"].Draw("col")
cnv_dE_vs_E.SaveAs(fn+"dE_vs_E.pdf")
cnv_dE_vs_E.SaveAs(allpdf)

cnv_res_xy_vs_xy = TCanvas("cnv_res_xy_vs_xy","",1000,500)
cnv_res_xy_vs_xy.Divide(2,1)
cnv_res_xy_vs_xy.cd(1)
histos["h_res_x_vs_xtru_recvstru"].Draw("col")
cnv_res_xy_vs_xy.cd(2)
histos["h_res_y_vs_ytru_recvstru"].Draw("col")
cnv_res_xy_vs_xy.SaveAs(fn+"res_xy_vs_xy.pdf")
cnv_res_xy_vs_xy.SaveAs(allpdf)

cnv_res_xy_vs_xy_log = TCanvas("cnv_res_xy_vs_xy_log_log","",1000,500)
cnv_res_xy_vs_xy_log.Divide(2,1)
cnv_res_xy_vs_xy_log.cd(1)
histos["h_res_x_vs_xtru_recvstru_logbins"].Draw("col")
cnv_res_xy_vs_xy_log.cd(2)
histos["h_res_y_vs_ytru_recvstru_logbins"].Draw("col")
cnv_res_xy_vs_xy_log.SaveAs(fn+"res_xy_vs_xy_log.pdf")
cnv_res_xy_vs_xy_log.SaveAs(allpdf)

cnv_dxy_vs_xy_log = TCanvas("cnv_dxy_vs_xy_log","",1000,1500)
cnv_dxy_vs_xy_log.Divide(2,3)
cnv_dxy_vs_xy_log.cd(1)
histos["h_dx_vs_xtru_clsvstru_logbins"].Draw("col")
cnv_dxy_vs_xy_log.cd(2)
histos["h_dy_vs_ytru_clsvstru_logbins"].Draw("col")
cnv_dxy_vs_xy_log.cd(3)
histos["h_dx_vs_xcls_recvscls_logbins"].Draw("col")
cnv_dxy_vs_xy_log.cd(4)
histos["h_dy_vs_ycls_recvscls_logbins"].Draw("col")
cnv_dxy_vs_xy_log.cd(5)
histos["h_dx_vs_xtru_recvstru_logbins"].Draw("col")
cnv_dxy_vs_xy_log.cd(6)
histos["h_dy_vs_ytru_recvstru_logbins"].Draw("col")
cnv_dxy_vs_xy_log.SaveAs(fn+"dxy_vs_xy_log.pdf")
cnv_dxy_vs_xy_log.SaveAs(allpdf)

cnv_dx_vs_invp_log = TCanvas("cnv_dx_vs_invp_log","",500,500)
cnv_dx_vs_invp_log.cd()
histos["h_dx_vs_invptru_recvstru_logbins"].Draw("col")
cnv_dx_vs_invp_log.SaveAs(fn+"dx_vs_invp_log.pdf")
cnv_dx_vs_invp_log.SaveAs(allpdf+")")


########### write all histos to a root file
allroot = storage+"/output/root/analysis_reco_"+process+".root"
tfileout = TFile(allroot,"RECREATE")
tfileout.cd()
for hname,hist in histos.items(): hist.Write()
tfileout.Write()
tfileout.Close()
