#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TLorentzVector, TVector3, TCanvas
import argparse
parser = argparse.ArgumentParser(description='truthphotons.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-e', metavar='energy',  required=False,  help='beam energy')
parser.add_argument('-g', metavar='photons?', required=False,help='photons only? [default=0/1]')
argus  = parser.parse_args()
proc   = argus.p
ebeam  = argus.e
photon = (argus.g=="1")

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage = os.path.expandvars("$STORAGEDIR")

histos = {}
R1 = 300 ## cm
R2 = 900 ## cm
um2cm = 1.e-4
###############################################################


## book
def Book(process):
   histos.update( {"h_E":TH1D("h_E",";E [GeV];Photons", 200,0,20)} )
   rmax = 0.05*um2cm if(process=="bppp") else 30.*um2cm
   zmax = 30*um2cm   if(process=="bppp") else 200.*um2cm
   histos.update({"h_xVtx":TH1D("h_xVtx",";Vertex(x) [cm];Photons",200,-rmax,+rmax)})
   histos.update({"h_yVtx":TH1D("h_yVtx",";Vertex(y) [cm];Photons",200,-rmax,+rmax)})
   histos.update({"h_zVtx":TH1D("h_zVtx",";Vertex(z) [cm];Photons",200,-zmax,+zmax)})
   histos.update({"h_xyVtx":TH2D("h_xyVtx",";Vertex(x) [cm];Vertex(y) [cm];Photons",200,-rmax,+rmax, 200,-rmax,+rmax)})
   histos.update({"h_R1_x":TH1D("h_R1_x",  "Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];Photons",500,-4,+4)})
   histos.update({"h_R1_y":TH1D("h_R1_y",  "Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};y [cm];Photons",500,-4,+4)})
   histos.update({"h_R1_xy":TH2D("h_R1_xy","Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];y [cm];Photons",500,-4,+4, 500,-4,+4)})
   histos.update({"h_R1_xE":TH2D("h_R1_xE","Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];E [GeV];Photons",500,-4,+4, 200,0,20)})
   histos.update({"h_R1_yE":TH2D("h_R1_yE","Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};y [cm];E [GeV];Photons",500,-4,+4, 200,0,20)})
   histos.update({"h_R2_x":TH1D("h_R2_x",  "Photons at r="+str(int(R2))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];Photons",500,-4,+4)})
   histos.update({"h_R2_y":TH1D("h_R2_y",  "Photons at r="+str(int(R2))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};y [cm];Photons",500,-4,+4)})
   histos.update({"h_R2_xy":TH2D("h_R2_xy","Photons at r="+str(int(R2))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];y [cm];Photons",500,-4,+4, 500,-4,+4)})
   histos.update({"h_R2_xE":TH2D("h_R2_xE","Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];E [GeV];Photons",500,-4,+4, 200,0,20)})
   histos.update({"h_R2_yE":TH2D("h_R2_yE","Photons at r="+str(int(R1))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};y [cm];E [GeV];Photons",500,-4,+4, 200,0,20)})
   ntrkmin  = 4780000 if(process=="bppp") else 175000000000
   ntrkmax  = 4810000 if(process=="bppp") else 270000000000
   ntrkbins = 300 if(process=="bppp")     else 200
   histos.update( {"h_photons":TH1D("h_photons",";Photons multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )

## analysis
def Analyze(n,event):
   nphot = 0
   for j in range(event.E.size()):
      wgt = event.wgt[j]
      nphot += wgt
      ### all photons
      histos["h_E"].Fill(event.E[j],wgt)
      p = TLorentzVector()
      p.SetPxPyPzE(event.px[j],event.py[j],event.pz[j],event.E[j])
      v = p.Vect()
      vxHat = (v.X()/v.Mag())
      vyHat = (v.Y()/v.Mag())
      vzHat = (v.Z()/v.Mag())

      r1 = TVector3(R1*vxHat,R1*vyHat,R1*vzHat)
      histos["h_R1_x"].Fill(r1.X(),wgt)
      histos["h_R1_y"].Fill(r1.Y(),wgt)
      histos["h_R1_xy"].Fill(r1.X(),r1.Y(),wgt)
      histos["h_R1_xE"].Fill(r1.X(),event.E[j],wgt)
      histos["h_R1_yE"].Fill(r1.Y(),event.E[j],wgt)
      
      r2 = TVector3(R2*vxHat,R2*vyHat,R2*vzHat)
      histos["h_R2_x"].Fill(r2.X(),wgt)
      histos["h_R2_y"].Fill(r2.Y(),wgt)
      histos["h_R2_xy"].Fill(r2.X(),r2.Y(),wgt)
      histos["h_R2_xE"].Fill(r2.X(),event.E[j],wgt)
      histos["h_R2_yE"].Fill(r2.Y(),event.E[j],wgt)

      xVtx = event.vx[j]
      yVtx = event.vy[j]
      zVtx = event.vz[j]
      histos["h_xVtx"].Fill(xVtx,wgt)
      histos["h_yVtx"].Fill(yVtx,wgt)
      histos["h_zVtx"].Fill(zVtx,wgt)
      histos["h_xyVtx"].Fill(xVtx,yVtx,wgt)
   histos["h_photons"].Fill(nphot)
   return nphot
   
## event loop
def EventLoop(tree):
   nevents = tree.GetEntries()
   print("with %d events" % nevents)
   nphot_min = 1e20
   nphot_max = -1
   n=0 ### init n
   for event in tree:
      if(n%100==0 and n>0): print("  processed %d events" % n)
      nphot = Analyze(n,event)
      nphot_min = nphot if(nphot<nphot_min) else nphot_min
      nphot_max = nphot if(nphot>nphot_max) else nphot_max
      n+=1
   print("Total events processed: ",n)
   print("nphot_min=%g, nphot_max=%g" % (nphot_min,nphot_max))
   return nevents

## file and tree
def Run(tfilename,ttreename):
   print("getting tree from ",tfilename)
   tfile = TFile(tfilename,"READ")
   tree = tfile.Get(ttreename)
   nevents = EventLoop(tree)
   return nevents

###############################################################

## actually run
process = proc
beamenergy = ebeam
tfilename = storage+"/data/root/raw_photons_"+process+".root"
Book(process)
nevents = Run(tfilename,"tt")
print("nevents=",nevents)

###############################################################

allpdf = storage+"/output/pdf/truthphotons_"+process+".pdf"
fn = allpdf.replace(".pdf","_")

## summarize
cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_E"].Draw("hist")
cnv.SaveAs(fn+"Energy.pdf")
cnv.SaveAs(allpdf+"(")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_photons"].Draw("hist")
cnv.SaveAs(fn+"nphotons.pdf")
cnv.SaveAs(allpdf)


cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_xVtx"].Draw("hist")
cnv.SaveAs(fn+"xVtx.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_yVtx"].Draw("hist")
cnv.SaveAs(fn+"yVtx.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_zVtx"].Draw("hist")
cnv.SaveAs(fn+"zVtx.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_xyVtx"].Draw("colz")
cnv.SaveAs(fn+"xyVtx.pdf")
cnv.SaveAs(allpdf)


cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_R1_x"].Draw("hist")
cnv.SaveAs(fn+"xAt3m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_R1_y"].Draw("hist")
cnv.SaveAs(fn+"yAt3m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R1_xy"].Draw("colz")
cnv.SaveAs(fn+"xyAt3m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R1_xE"].Draw("colz")
cnv.SaveAs(fn+"xeAt3m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R1_yE"].Draw("colz")
cnv.SaveAs(fn+"yeAt3m.pdf")
cnv.SaveAs(allpdf)


cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_R2_x"].Draw("hist")
cnv.SaveAs(fn+"xAt9m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogy()
histos["h_R2_y"].Draw("hist")
cnv.SaveAs(fn+"yAt9m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R2_xy"].Draw("colz")
cnv.SaveAs(fn+"xyAt9m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R2_xE"].Draw("colz")
cnv.SaveAs(fn+"xeAt9m.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
cnv.SetLogz()
histos["h_R2_yE"].Draw("colz")
cnv.SaveAs(fn+"yeAt9m.pdf")
cnv.SaveAs(allpdf+")")