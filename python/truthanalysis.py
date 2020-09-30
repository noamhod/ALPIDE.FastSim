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
storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")

histos = {}
R = 300 ## cm
meGeV = 0.5109989461/1000.
um2cm = 1.e-4
###############################################################

process = proc
basepath = storage+"/data/root/raw/"+path
targetdir_pdf  = basepath.replace("data/root","output/pdf")
targetdir_root = basepath.replace("data/root","output/root")
p = subprocess.Popen("mkdir -p "+targetdir_pdf, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()
p = subprocess.Popen("mkdir -p "+targetdir_root, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()

allpdf  = targetdir_pdf+"/truthtest_"+process+".pdf"
allroot = targetdir_root+"/truthtest_"+process+".root"

###############################################################
### get the resolution and acc*eff
tfaccxeff    = TFile("$STORAGEDIR/output/root/analysis_reco_bppp_accxeff.root","READ")
tfresolution = TFile("$STORAGEDIR/output/root/analysis_reco_bppp_Eresolution.root","READ")
accxeff    = tfaccxeff.Get("turnon")
resolution = tfresolution.Get("fite")


###############################################################

def minmax(h1,h2,f=1.1):
   hmin = h1.GetMinimum() if(h1.GetMinimum()<h2.GetMinimum()) else h2.GetMinimum()
   hmax = h1.GetMaximum() if(h1.GetMaximum()>h2.GetMaximum()) else h2.GetMaximum()
   h1.SetMaximum(hmax*f)
   h2.SetMaximum(hmax*f)
   return hmin,hmax*f

## book
def Book(process):
   ntotmax = 6+1 if(process=="bppp") else 3+1
   histos.update( {"h_ntot":TH1D("h_ntot",";;N", ntotmax,0,ntotmax)} )
   histos["h_ntot"].GetXaxis().SetBinLabel(1,"N_{BXs}")
   histos["h_ntot"].GetXaxis().SetBinLabel(2,"N_{e+} tru_all")
   histos["h_ntot"].GetXaxis().SetBinLabel(3,"N_{e+} tru_acc")
   histos["h_ntot"].GetXaxis().SetBinLabel(4,"N_{e+} rec_emu")
   if(process=="bppp"):
      histos["h_ntot"].GetXaxis().SetBinLabel(5,"N_{e-} tru_all")
      histos["h_ntot"].GetXaxis().SetBinLabel(6,"N_{e-} tru_acc")
      histos["h_ntot"].GetXaxis().SetBinLabel(7,"N_{e-} rec_emu")

   ntrkmin = 0
   ntrkmax = 200
   ntrkbins = 200
   histos.update( {"h_ntrks_positrons":TH1D("h_ntrks_positrons",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_electrons":TH1D("h_ntrks_electrons",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_positrons_rec_emul":TH1D("h_ntrks_positrons_rec_emul",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins,ntrkmin,ntrkmax)} )
   histos.update( {"h_ntrks_electrons_rec_emul":TH1D("h_ntrks_electrons_rec_emul",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins,ntrkmin,ntrkmax)} )

   ntrkmax_small = 50
   ntrkbins_small = 50
   histos.update( {"h_ntrks_positrons_small":TH1D("h_ntrks_positrons_small",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins_small,ntrkmin,ntrkmax_small)} )
   histos.update( {"h_ntrks_electrons_small":TH1D("h_ntrks_electrons_small",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins_small,ntrkmin,ntrkmax_small)} )
   histos.update( {"h_ntrks_positrons_rec_emul_small":TH1D("h_ntrks_positrons_rec_emul_small",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins_small,ntrkmin,ntrkmax_small)} )
   histos.update( {"h_ntrks_electrons_rec_emul_small":TH1D("h_ntrks_electrons_rec_emul_small",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins_small,ntrkmin,ntrkmax_small)} )

   ntrkmax_tiny = 25
   ntrkbins_tiny = 25
   histos.update( {"h_ntrks_positrons_tiny":TH1D("h_ntrks_positrons_tiny",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins_tiny,ntrkmin,ntrkmax_tiny)} )
   histos.update( {"h_ntrks_electrons_tiny":TH1D("h_ntrks_electrons_tiny",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins_tiny,ntrkmin,ntrkmax_tiny)} )
   histos.update( {"h_ntrks_positrons_rec_emul_tiny":TH1D("h_ntrks_positrons_rec_emul_tiny",";Positron multiplicity;N_{e+}/BX/Shot", ntrkbins_tiny,ntrkmin,ntrkmax_tiny)} )
   histos.update( {"h_ntrks_electrons_rec_emul_tiny":TH1D("h_ntrks_electrons_rec_emul_tiny",";Electron multiplicity;N_{e-}/BX/Shot", ntrkbins_tiny,ntrkmin,ntrkmax_tiny)} )
   
   Emin = 0
   Emax = 20 if(proc=="bppp") else 10
   Emax_full = 20
   Ebins = 80 if(proc=="bppp") else 40
   Ebins_fine = 200 if(proc=="bppp") else 100
   Ebins_full = 200
   histos.update( {"h_E_positrons"     :TH1D("h_E_positrons",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins,Emin,Emax)} )
   histos.update( {"h_E_positrons_fine":TH1D("h_E_positrons_fine",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins_fine,Emin,Emax)} )
   histos.update( {"h_E_electrons"     :TH1D("h_E_electrons",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins,Emin,Emax)} )
   histos.update( {"h_E_electrons_fine":TH1D("h_E_electrons_fine",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins_fine,Emin,Emax)} )
   histos.update( {"h_E_positrons_rec_emul"     :TH1D("h_E_positrons_rec_emul",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins,Emin,Emax)} )
   histos.update( {"h_E_positrons_rec_emul_fine":TH1D("h_E_positrons_rec_emul_fine",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins_fine,Emin,Emax)} )
   histos.update( {"h_E_electrons_rec_emul"     :TH1D("h_E_electrons_rec_emul",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins,Emin,Emax)} )
   histos.update( {"h_E_electrons_rec_emul_fine":TH1D("h_E_electrons_rec_emul_fine",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins_fine,Emin,Emax)} )
   histos.update( {"h_E_positrons_full"     :TH1D("h_E_positrons_full",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins_full,Emin,Emax_full)} )
   histos.update( {"h_E_electrons_full"     :TH1D("h_E_electrons_full",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins_full,Emin,Emax_full)} )
   histos.update( {"h_E_positrons_rec_emul_full" :TH1D("h_E_positrons_rec_emul_full",";#it{E} [GeV];N_{e+}/BX/Shot",Ebins_full,Emin,Emax_full)} )
   histos.update( {"h_E_electrons_rec_emul_full" :TH1D("h_E_electrons_rec_emul_full",";#it{E} [GeV];N_{e-}/BX/Shot",Ebins_full,Emin,Emax_full)} )

   pzmin = 0
   pzmax = 20 if(proc=="bppp") else 10
   pzmax_full = 20
   pzbins = 80 if(proc=="bppp") else 40
   pzbins_fine = 200 if(proc=="bppp") else 100
   pzbins_full = 200
   histos.update( {"h_pz_positrons"     :TH1D("h_pz_positrons",";#it{p}_{z} [GeV];N_{e+}/BX/Shot",pzbins,pzmin,pzmax)} )
   histos.update( {"h_pz_positrons_fine":TH1D("h_pz_positrons_fine",";#it{p}_{z} [GeV];N_{e+}/BX/Shot",pzbins_fine,pzmin,pzmax)} )
   histos.update( {"h_pz_positrons_full":TH1D("h_pz_positrons_full",";#it{p}_{z} [GeV];N_{e+}/BX/Shot",pzbins_full,pzmin,pzmax_full)} )
   histos.update( {"h_pz_electrons"     :TH1D("h_pz_electrons",";#it{p}_{z} [GeV];N_{e-}/BX/Shot",pzbins,pzmin,pzmax)} )
   histos.update( {"h_pz_electrons_fine":TH1D("h_pz_electrons_fine",";#it{p}_{z} [GeV];N_{e-}/BX/Shot",pzbins_fine,pzmin,pzmax)} )
   histos.update( {"h_pz_electrons_full":TH1D("h_pz_electrons_full",";#it{p}_{z} [GeV];N_{e-}/BX/Shot",pzbins_full,pzmin,pzmax_full)} )

   pymin = -0.005
   pymax = +0.005
   pybins = 100
   pybins_fine = 200
   histos.update( {"h_py_positrons"     :TH1D("h_py_positrons",";#it{p}_{y} [GeV];N_{e+}/BX/Shot",pybins,pymin,pymax)} )
   histos.update( {"h_py_positrons_fine":TH1D("h_py_positrons_fine",";#it{p}_{y} [GeV];N_{e+}/BX/Shot",pybins_fine,pymin,pymax)} )
   histos.update( {"h_py_electrons"     :TH1D("h_py_electrons",";#it{p}_{y} [GeV];N_{e-}/BX/Shot",pybins,pymin,pymax)} )
   histos.update( {"h_py_electrons_fine":TH1D("h_py_electrons_fine",";#it{p}_{y} [GeV];N_{e-}/BX/Shot",pybins_fine,pymin,pymax)} )

   pxmin = -0.005
   pxmax = +0.005
   pxbins = 100
   pxbins_fine = 200
   histos.update( {"h_px_positrons"     :TH1D("h_px_positrons",";#it{p}_{x} [GeV];N_{e+}/BX/Shot",pxbins,pxmin,pxmax)} )
   histos.update( {"h_px_positrons_fine":TH1D("h_px_positrons_fine",";#it{p}_{x} [GeV];N_{e+}/BX/Shot",pxbins_fine,pxmin,pxmax)} )
   histos.update( {"h_px_electrons"     :TH1D("h_px_electrons",";#it{p}_{x} [GeV];N_{e-}/BX/Shot",pxbins,pxmin,pxmax)} )
   histos.update( {"h_px_electrons_fine":TH1D("h_px_electrons_fine",";#it{p}_{x} [GeV];N_{e-}/BX/Shot",pxbins_fine,pxmin,pxmax)} )
   
   rmax = 0.05*um2cm if(process=="bppp") else 30.*um2cm
   zmax = 30*um2cm   if(process=="bppp") else 200.*um2cm
   histos.update({"h_xyVtx_positrons": TH2D("h_xyVtx_positrons","Positrons vertex x:y;Vertex(x) [cm];Vertex(y) [cm];N_{e+}/BX/Shot",200,-rmax,+rmax, 200,-rmax,+rmax)})
   histos.update({"h_zxVtx_positrons": TH2D("h_zxVtx_positrons","Positrons vertex z:x;Vertex(z) [cm];Vertex(x) [cm];N_{e+}/BX/Shot",200,-zmax,+zmax, 200,-rmax,+rmax)})
   histos.update({"h_zyVtx_positrons": TH2D("h_zyVtx_positrons","Positrons vertex z:y;Vertex(z) [cm];Vertex(y) [cm];N_{e+}/BX/Shot",200,-zmax,+zmax, 200,-rmax,+rmax)})

   ## sumw2
   for hname,hist in histos.items(): hist.Sumw2()

###############################################################

## analysis
def Analyze(n,event):
   npositrons = 0
   nelectrons = 0
   npositrons_acc = 0
   nelectrons_acc = 0
   npositrons_rec_emul = 0
   nelectrons_rec_emul = 0
   ### loop on all tracks
   for j in range(event.E.size()):
      wgt = event.wgt[j]
      p = TLorentzVector()
      p.SetXYZM(event.px[j],event.py[j],event.pz[j],meGeV)
      xVtx = event.vx[j]
      yVtx = event.vy[j]
      zVtx = event.vz[j]
      Etru = event.E[j]
      Erec = Etru*( 1+resolution.GetRandom() )
      prob = accxeff.Eval(Etru)
      
      ### only electrons
      if(event.pdgId[j]==11):
         nelectrons+=wgt
         if(Etru>1.5): nelectrons_acc+=wgt
         nelectrons_rec_emul += prob*wgt
         histos["h_E_electrons"].Fill(Etru,wgt)
         histos["h_E_electrons_fine"].Fill(Etru,wgt)         
         histos["h_E_electrons_full"].Fill(Etru,wgt)         
         histos["h_px_electrons"].Fill(event.px[j],wgt)
         histos["h_px_electrons_fine"].Fill(event.px[j],wgt)
         histos["h_py_electrons"].Fill(event.py[j],wgt)
         histos["h_py_electrons_fine"].Fill(event.py[j],wgt)
         histos["h_pz_electrons"].Fill(event.pz[j],wgt)
         histos["h_pz_electrons_fine"].Fill(event.pz[j],wgt)
         histos["h_pz_electrons_full"].Fill(event.pz[j],wgt)
         histos["h_E_electrons_rec_emul"].Fill(Erec,wgt*prob)
         histos["h_E_electrons_rec_emul_fine"].Fill(Erec,wgt*prob)
         histos["h_E_electrons_rec_emul_full"].Fill(Erec,wgt*prob)
         
      ### only positrons
      if(event.pdgId[j]==-11):
         npositrons+=wgt
         if(Etru>1.5): npositrons_acc+=wgt
         npositrons_rec_emul += prob*wgt
         histos["h_E_positrons"].Fill(Etru,wgt)
         histos["h_E_positrons_fine"].Fill(Etru,wgt)
         histos["h_E_positrons_full"].Fill(Etru,wgt)
         histos["h_px_positrons"].Fill(event.px[j],wgt)
         histos["h_px_positrons_fine"].Fill(event.px[j],wgt)
         histos["h_py_positrons"].Fill(event.py[j],wgt)
         histos["h_py_positrons_fine"].Fill(event.py[j],wgt)
         histos["h_pz_positrons"].Fill(event.pz[j],wgt)
         histos["h_pz_positrons_fine"].Fill(event.pz[j],wgt)
         histos["h_pz_positrons_full"].Fill(event.pz[j],wgt)
         histos["h_xyVtx_positrons"].Fill(xVtx,yVtx)
         histos["h_zxVtx_positrons"].Fill(zVtx,xVtx)
         histos["h_zyVtx_positrons"].Fill(zVtx,yVtx)
         histos["h_E_positrons_rec_emul"].Fill(Erec,wgt*prob)
         histos["h_E_positrons_rec_emul_fine"].Fill(Erec,wgt*prob)
         histos["h_E_positrons_rec_emul_full"].Fill(Erec,wgt*prob)
   
   histos["h_ntrks_positrons"].Fill(npositrons)
   histos["h_ntrks_positrons_small"].Fill(npositrons)
   histos["h_ntrks_positrons_tiny"].Fill(npositrons)
   
   histos["h_ntrks_positrons_rec_emul"].Fill(npositrons_rec_emul)
   histos["h_ntrks_positrons_rec_emul_small"].Fill(npositrons_rec_emul)
   histos["h_ntrks_positrons_rec_emul_tiny"].Fill(npositrons_rec_emul)
   
   histos["h_ntrks_electrons"].Fill(nelectrons)
   histos["h_ntrks_electrons_small"].Fill(nelectrons)
   histos["h_ntrks_electrons_tiny"].Fill(nelectrons)
   
   histos["h_ntrks_electrons_rec_emul"].Fill(nelectrons_rec_emul)
   histos["h_ntrks_electrons_rec_emul_small"].Fill(nelectrons_rec_emul)
   histos["h_ntrks_electrons_rec_emul_tiny"].Fill(nelectrons_rec_emul)
   
   histos["h_ntot"].AddBinContent(2,npositrons)
   histos["h_ntot"].AddBinContent(3,npositrons_acc)
   histos["h_ntot"].AddBinContent(4,npositrons_rec_emul)
   if(process=="bppp"):
      histos["h_ntot"].AddBinContent(5,nelectrons)
      histos["h_ntot"].AddBinContent(6,nelectrons_acc)
      histos["h_ntot"].AddBinContent(7,nelectrons_rec_emul)

###############################################################

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

###############################################################

## scale to BX/Shot
for hname,hist in histos.items():
   if(hname=="h_ntot"): hist.AddBinContent(nevents) ## keep track on number of events
   else:                hist.Scale(1./nevents)      ## otherwise, normalise to number of events

###############################################################

##### summarize pdf plots
fn = allpdf.replace(".pdf","_")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
hmin,hmax = minmax(histos["h_ntrks_positrons"],histos["h_ntrks_positrons_rec_emul"])
histos["h_ntrks_positrons"].Draw("hist")
histos["h_ntrks_positrons_rec_emul"].SetLineColor(ROOT.kRed)
histos["h_ntrks_positrons_rec_emul"].Draw("hist same")
cnv.SaveAs(fn+"ntrks.pdf")
cnv.SaveAs(allpdf+"(")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
hmin,hmax = minmax(histos["h_ntrks_positrons_small"],histos["h_ntrks_positrons_rec_emul_small"])
histos["h_ntrks_positrons_small"].Draw("hist")
histos["h_ntrks_positrons_rec_emul_small"].SetLineColor(ROOT.kRed)
histos["h_ntrks_positrons_rec_emul_small"].Draw("hist same")
cnv.SaveAs(fn+"ntrks.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",1000,1000)
cnv.Divide(2,2)
cnv.cd(1)
hmin,hmax = minmax(histos["h_E_positrons"],histos["h_E_positrons_rec_emul"])
histos["h_E_positrons"].Draw("hist")
histos["h_E_positrons_rec_emul"].SetLineColor(ROOT.kRed)
histos["h_E_positrons_rec_emul"].Draw("hist same")
cnv.cd(2)
hmin,hmax = minmax(histos["h_E_positrons_fine"],histos["h_E_positrons_rec_emul_fine"])
histos["h_E_positrons_fine"].Draw("hist")
histos["h_E_positrons_rec_emul_fine"].SetLineColor(ROOT.kRed)
histos["h_E_positrons_rec_emul_fine"].Draw("hist same")
cnv.cd(3)
histos["h_E_electrons_full"].Draw("hist")
cnv.cd(4)
histos["h_E_electrons_full"].Draw("hist")
cnv.SaveAs(fn+"Energy.pdf")
cnv.SaveAs(allpdf)


cnv = TCanvas("cnv","",1500,1000)
cnv.Divide(3,2)
cnv.cd(1)
histos["h_px_positrons"].Draw("hist")
cnv.cd(2)
histos["h_py_positrons"].Draw("hist")
cnv.cd(3)
histos["h_pz_positrons"].Draw("hist")
cnv.cd(4)
histos["h_px_electrons"].Draw("hist")
cnv.cd(5)
histos["h_py_electrons"].Draw("hist")
cnv.cd(6)
histos["h_pz_electrons"].Draw("hist")
cnv.SaveAs(fn+"Momenta.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",1500,1000)
cnv.Divide(3,2)
cnv.cd(1)
histos["h_px_positrons_fine"].Draw("hist")
cnv.cd(2)
histos["h_py_positrons_fine"].Draw("hist")
cnv.cd(3)
histos["h_pz_positrons_fine"].Draw("hist")
cnv.cd(4)
histos["h_px_electrons_fine"].Draw("hist")
cnv.cd(5)
histos["h_py_electrons_fine"].Draw("hist")
cnv.cd(6)
histos["h_pz_electrons_fine"].Draw("hist")
cnv.SaveAs(fn+"Momenta_fine.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
histos["h_xyVtx_positrons"].Draw("colz")
cnv.cd(2)
histos["h_zxVtx_positrons"].Draw("colz")
cnv.cd(3)
histos["h_zyVtx_positrons"].Draw("colz")
cnv.SaveAs(fn+"vtx.pdf")
cnv.SaveAs(allpdf+")")


######################################
tfileout = TFile(allroot,"RECREATE")
tfileout.cd()
for hname,hist in histos.items(): hist.Write()
tfileout.Write()
tfileout.Close()