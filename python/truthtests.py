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

histos = {}
R = 300 ## cm
meGeV = 0.5109989461/1000.
um2cm = 1.e-4
###############################################################


## book
def Book(process):
   Eeemax   = 20  if(process=="bppp") else 25
   nEeebins = 200 if(process=="bppp") else 250
   histos.update( {"h_Eee":TH1D("h_Eee",";E_{ee} [GeV];Tracks", nEeebins,0,Eeemax)} )
   histos.update( {"h_Mee":TH1D("h_Mee",";m_{ee} [GeV];Tracks", 100,0,+0.01)} )
   dEeemax = 15 if(process=="bppp") else 20
   histos.update( {"h_dEee":TH1D("h_dEee",";|E_{1}-E_{2}| [GeV];Tracks", 100,0,dEeemax)} )
   histos.update( {"h_E":TH1D("h_E",";E [GeV];Tracks", 200,0,20)} )
   histos.update( {"h_Ee+":TH1D("h_Ee+","Positrons only;E [GeV];Tracks", 200,0,20)} )
   ntrkmax = 200  if(process=="bppp") else 1250
   ntrkmin = 0    if(process=="bppp") else 1000
   ntrkbins = 100 if(process=="bppp") else 50
   histos.update( {"h_ntrks":TH1D("h_ntrks",";Track multiplicity;Events", ntrkbins,ntrkmin,ntrkmax)} )
   rmax = 0.05*um2cm if(process=="bppp") else 30.*um2cm
   zmax = 30*um2cm   if(process=="bppp") else 200.*um2cm
   histos.update({"h_xVtx":TH1D("h_xVtx","Positrons only;Vertex(x) [cm];Tracks",200,-rmax,+rmax)})
   histos.update({"h_yVtx":TH1D("h_yVtx","Positrons only;Vertex(y) [cm];Tracks",200,-rmax,+rmax)})
   histos.update({"h_zVtx":TH1D("h_zVtx","Positrons only;Vertex(z) [cm];Tracks",200,-zmax,+zmax)})
   histos.update({"h_xyVtx":TH2D("h_xyVtx","Positrons only;Vertex(x) [cm];Vertex(y) [cm];Tracks",200,-rmax,+rmax, 200,-rmax,+rmax)})
   histos.update({"h_x":TH1D("h_x","Positrons only, at r="+str(int(R))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];Tracks",200,-1,+1)})
   histos.update({"h_y":TH1D("h_y","Positrons only, at r="+str(int(R))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};y [cm];Tracks",200,-1,+1)})
   histos.update({"h_xy":TH2D("h_xy","Positrons only, at r="+str(int(R))+" cm, with #hat{r}=(p_{x},p_{y},p_{z})/|#vec{p}| and #vec{r}=r#upoint#hat{r};x [cm];y [cm];Tracks",200,-1,+1, 200,-1,+1)})


## analysis
def Analyze(n,event):
   histos["h_ntrks"].Fill(event.E.size())
   npairs = 0
   npositrons = 0
   for j in range(event.E.size()):
      wgt = event.wgt[j]
      ### all e+/e-
      histos["h_E"].Fill(event.E[j],wgt)
      if(event.pdgId[j]<0): histos["h_Ee+"].Fill(event.E[j],wgt)
      ### e+e- quantities
      if(j<event.E.size()-1):
         isel = (event.pdgId[j]==11) ## avoid double counting
         is0q = (event.pdgId[j]+event.pdgId[j+1]==0)
         is1w = (int(event.wgt[j])==1 and int(event.wgt[j+1])==1)
         if(isel and is0q and is1w):
            npairs+=1
            p1 = TLorentzVector()
            p2 = TLorentzVector()
            # p1.SetPxPyPzE(event.px[j],event.py[j],event.pz[j],event.E[j])
            # p2.SetPxPyPzE(event.px[j+1],event.py[j+1],event.pz[j+1],event.E[j+1])
            p1.SetXYZM(event.px[j],event.py[j],event.pz[j],meGeV)
            p2.SetXYZM(event.px[j+1],event.py[j+1],event.pz[j+1],meGeV)
            p = p1+p2
            wgt12 = wgt*event.wgt[j+1]
            histos["h_Eee"].Fill(p.E(),wgt12)
            histos["h_Mee"].Fill(p.M(),wgt12)
            histos["h_dEee"].Fill(abs(event.E[j]-event.E[j+1]),wgt12)
      ### only positrons
      if(event.pdgId[j]==-11):
         npositrons+=1
         p = TLorentzVector()
         # p.SetPxPyPzE(event.px[j],event.py[j],event.pz[j],event.E[j])
         p.SetXYZM(event.px[j],event.py[j],event.pz[j],meGeV)
         v = p.Vect()
         vxHat = (v.X()/v.Mag())
         vyHat = (v.Y()/v.Mag())
         vzHat = (v.Z()/v.Mag())
         r = TVector3(R*vxHat,R*vyHat,R*vzHat)
         histos["h_x"].Fill(r.X())
         histos["h_y"].Fill(r.Y())
         histos["h_xy"].Fill(r.X(),r.Y())
         xVtx = event.vx[j]
         yVtx = event.vy[j]
         zVtx = event.vz[j]
         histos["h_xVtx"].Fill(xVtx)
         histos["h_yVtx"].Fill(yVtx)
         histos["h_zVtx"].Fill(zVtx)
         # print(zVtx)
         histos["h_xyVtx"].Fill(xVtx,yVtx)

   
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

###############################################################

## actually run
process = proc
beamenergy = ebeam
tfilename = "../raw_"+process+".root"
Book(process)
nevents = Run(tfilename,"tt")
print("nevents=",nevents)

###############################################################

allpdf = "../output/pdf/truthtest_"+process+".pdf"
fn = allpdf.replace(".pdf","_")

## summarize
cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_E"].Draw("hist")
cnv.SaveAs(fn+"Energy.pdf")
cnv.SaveAs(allpdf+"(")

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_Eee"].Draw("hist")
cnv.SaveAs(fn+"Eee.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_Mee"].Draw("hist")
cnv.SaveAs(fn+"Mee.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_dEee"].Draw("hist")
cnv.SaveAs(fn+"dEee.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_ntrks"].Draw("hist")
cnv.SaveAs(fn+"ntrks.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_Ee+"].Draw("hist")
cnv.SaveAs(fn+"Epositrons.pdf")
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
histos["h_x"].Draw("hist")
cnv.SaveAs(fn+"xAtTracker.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_y"].Draw("hist")
cnv.SaveAs(fn+"yAtTracker.pdf")
cnv.SaveAs(allpdf)

cnv = TCanvas("cnv","",500,500)
cnv.cd()
histos["h_xy"].Draw("colz")
cnv.SaveAs(fn+"xyAtTracker.pdf")
cnv.SaveAs(allpdf+")")