#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TPad, TView, TLatex, TLegend, TGaxis
import argparse
parser = argparse.ArgumentParser(description='analysis.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-b', metavar='Bfield',  required=True,  help='Magnetic filed in kG [0-20]')
parser.add_argument('-n', metavar='nevents', required=True,  help='maximum number of events')
argus = parser.parse_args()
proc  = argus.p
bfield = argus.b+"kG"
btitle = "%.1f [T]" % (float(argus.b)/10)
ptitle = "BPPP" if(proc=="bppp") else "Trident"
Nmax = int(argus.n)
if(Nmax==0): Nmax = 1e10

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
# ROOT.gStyle.SetPadBottomMargin(0.15)
# ROOT.gStyle.SetPadLeftMargin(0.16)
   
#############################################

process = proc

### stave geometry
np=-1
Hstave = 1.5  # cm
Lstave = 50 #27.12
Rbeampipe = 2.413 # cm
RoffsetBfield = 5.7 if(process=="bppp") else 14 # cm
x1L = -RoffsetBfield-Lstave 
x1R = -RoffsetBfield        
x2L = +RoffsetBfield        
x2R = +RoffsetBfield+Lstave 
yUp = +Hstave/2.        
yDn = -Hstave/2.     
xStavesOverlap = 4
zStavesOffset = 1.2
xInnerLlist = [x1L,x1R]
xInnerRlist = [x2R,x2L]
xOuterLlist = [x1L+xStavesOverlap-Lstave,x1R+xStavesOverlap-Lstave]
xOuterRlist = [x2R-xStavesOverlap+Lstave,x2L-xStavesOverlap+Lstave]
ylist  = [yDn,yUp]
xInL = array.array('d', xInnerLlist)
xInR = array.array('d', xInnerRlist)
xOutL = array.array('d', xOuterLlist)
xOutR = array.array('d', xOuterRlist)
zIn1 = array.array('d', [300,300])
zIn2 = array.array('d', [310,310])
zIn3 = array.array('d', [320,320])
zIn4 = array.array('d', [330,330])
zOut1 = array.array('d', [300+zStavesOffset,300+zStavesOffset])
zOut2 = array.array('d', [310+zStavesOffset,310+zStavesOffset])
zOut3 = array.array('d', [320+zStavesOffset,320+zStavesOffset])
zOut4 = array.array('d', [330+zStavesOffset,330+zStavesOffset])
yS  = array.array('d', ylist)
np = len(xInnerLlist)
staveInnerLXZ = [TPolyLine(np,xInL,zIn1),TPolyLine(np,xInL,zIn2),TPolyLine(np,xInL,zIn3),TPolyLine(np,xInL,zIn4)]
staveInnerRXZ = [TPolyLine(np,xInR,zIn1),TPolyLine(np,xInR,zIn2),TPolyLine(np,xInR,zIn3),TPolyLine(np,xInR,zIn4)]
staveOuterLXZ = [TPolyLine(np,xOutL,zOut1),TPolyLine(np,xOutL,zOut2),TPolyLine(np,xOutL,zOut3),TPolyLine(np,xOutL,zOut4)]
staveOuterRXZ = [TPolyLine(np,xOutR,zOut1),TPolyLine(np,xOutR,zOut2),TPolyLine(np,xOutR,zOut3),TPolyLine(np,xOutR,zOut4)]
staveYZ  = [TPolyLine(np,yS,zIn1),TPolyLine(np,yS,zIn2),TPolyLine(np,yS,zIn3),TPolyLine(np,yS,zIn4)]
for i in range(len(staveInnerLXZ)):
   staveInnerLXZ[i].SetLineColor(ROOT.kCyan)
   staveInnerRXZ[i].SetLineColor(ROOT.kCyan)
   staveOuterLXZ[i].SetLineColor(ROOT.kCyan)
   staveOuterRXZ[i].SetLineColor(ROOT.kCyan)
   staveYZ[i].SetLineColor(ROOT.kCyan)


wDipoleInX = 330/10
hDipoleInY = 108/10
lDipoleInZ = 1029/10
zDipoleInOffset = 100
xDipoleIn = array.array('d', [-wDipoleInX/2,-wDipoleInX/2,+wDipoleInX/2,+wDipoleInX/2,-wDipoleInX/2])
yDipoleIn = array.array('d', [-hDipoleInY/2,-hDipoleInY/2,+hDipoleInY/2,+hDipoleInY/2,-hDipoleInY/2])
zDipoleIn = array.array('d', [zDipoleInOffset,zDipoleInOffset+lDipoleInZ,zDipoleInOffset+lDipoleInZ,zDipoleInOffset,zDipoleInOffset])
np = len(xDipoleIn)
dipoleInXZ = TPolyLine(np,xDipoleIn,zDipoleIn)
dipoleInYZ = TPolyLine(np,yDipoleIn,zDipoleIn)
dipoleInXZ.SetLineColor(ROOT.kRed)
dipoleInYZ.SetLineColor(ROOT.kRed)

wDipoleOutX = 1196/10
hDipoleOutY = 672/10
lDipoleOutZ = 1396/10
lDipoleInZDummy = 1238/10
deltax = ((lDipoleOutZ-lDipoleInZDummy)/2)*math.tan(30*math.pi/180)
deltaz = ((lDipoleOutZ-lDipoleInZDummy)/2)
zDipoleOutOffset = zDipoleInOffset - (lDipoleOutZ-lDipoleInZ)/2
zDipoleInDummyOffset = zDipoleOutOffset + (lDipoleOutZ-lDipoleInZDummy)/2
xDipoleOutL = array.array('d', [-wDipoleOutX/2,   -wDipoleOutX/2,               -wDipoleInX/2-deltax,         -wDipoleInX/2-deltax/2,                    -wDipoleInX/2-deltax/2,              -wDipoleInX/2,                        -wDipoleInX/2,        -wDipoleInX/2,        -wDipoleInX/2-deltax, -wDipoleOutX/2, ])
zDipoleOutL = array.array('d', [zDipoleOutOffset, zDipoleOutOffset+lDipoleOutZ, zDipoleOutOffset+lDipoleOutZ, zDipoleOutOffset+lDipoleOutZ-deltaz*(2/3), zDipoleOutOffset+lDipoleOutZ-deltaz, zDipoleInDummyOffset+lDipoleInZDummy, zDipoleInDummyOffset, zDipoleInDummyOffset, zDipoleOutOffset,     zDipoleOutOffset])
xDipoleOutR = array.array('d', [wDipoleOutX/2,    wDipoleOutX/2,                wDipoleInX/2+deltax,          wDipoleInX/2+deltax/2,                     wDipoleInX/2+deltax/2,               wDipoleInX/2,                         wDipoleInX/2,         wDipoleInX/2,         wDipoleInX/2+deltax,  wDipoleOutX/2, ])
zDipoleOutR = array.array('d', [zDipoleOutOffset, zDipoleOutOffset+lDipoleOutZ, zDipoleOutOffset+lDipoleOutZ, zDipoleOutOffset+lDipoleOutZ-deltaz*(2/3), zDipoleOutOffset+lDipoleOutZ-deltaz, zDipoleInDummyOffset+lDipoleInZDummy, zDipoleInDummyOffset, zDipoleInDummyOffset, zDipoleOutOffset,     zDipoleOutOffset])

yDipoleOut = array.array('d', [-hDipoleOutY/2,-hDipoleOutY/2,+hDipoleOutY/2,+hDipoleOutY/2,-hDipoleOutY/2])
zDipoleOut = array.array('d', [zDipoleOutOffset,zDipoleOutOffset+lDipoleOutZ,zDipoleOutOffset+lDipoleOutZ,zDipoleOutOffset,zDipoleOutOffset])
np = len(zDipoleOut)
dipoleOutYZ = TPolyLine(np,yDipoleOut,zDipoleOut)
dipoleOutYZ.SetLineColor(ROOT.kBlack)

np = len(zDipoleOutL)
dipoleOutXZL = TPolyLine(np,xDipoleOutL,zDipoleOutL)
dipoleOutXZL.SetLineColor(ROOT.kBlack)
dipoleOutXZR = TPolyLine(np,xDipoleOutR,zDipoleOutR)
dipoleOutXZR.SetLineColor(ROOT.kBlack)

dipoleExit = TPolyLine(2,array.array('d',[-100,+100]),array.array('d',[zDipoleOutOffset+lDipoleOutZ,zDipoleOutOffset+lDipoleOutZ]))

##################################################################

tfile = TFile("../data/root/Bfields/dig_"+process+"_"+bfield+".root","READ")
fn = "../output/pdf/BfieldStudy_"+process+"_"+bfield+"_"
fr = fn
fr = fr.replace("pdf","root")

# the inner magnet volume is: x*y*z = 326*108*1029 mm^3
# the outer magnet dimensions are: x*y*z= 1196*672*1520 mm^3

hxE = TH2D("h2_E_vs_x",";x [cm];E [GeV];Tracks", 1000,-100,+100, 200, 0,  20)
hxz = TH2D("h2_z_vs_x",";x [cm];z [cm];Tracks",  1000,-100,+100, 2000,0,+400)
hyz = TH2D("h2_z_vs_y",";y [cm];z [cm];Tracks",  1000,-100,+100, 2000,0,+400)

ttree = tfile.Get("dig")
Ntotal = ttree.GetEntries()
Nprint = 2
Nevents = 0
for event in ttree:
   for i in range(event.trkpts_fullrange.size()):
      for n in range(event.trkpts_fullrange[i].GetN()):
         x = ROOT.Double()
         y = ROOT.Double()
         z = ROOT.Double()
         event.trkpts_fullrange[i].GetPoint(n,x,y,z)
         E = event.trkp4[i].E()
         hxz.Fill(x,z)
         hyz.Fill(y,z)
         if(abs(z-(zDipoleOutOffset+lDipoleOutZ))<=0.1): hxE.Fill(x,E)
   if(Nevents%Nprint==0): print("Done %g out of %g" % (Nevents,Ntotal))
   if(Nevents==Nmax): break
   Nevents += 1

##################################################################

fOut = TFile(fr+"all.root","RECREATE")
fOut.cd()


cnv = TCanvas("cnv_xyzE","",1000,1000)
cnv.Divide(2,1)
cnv.cd(1)
pad1 = TPad("pad1","",0,0,1,1)
pad2 = TPad("pad2","",0,0,1,1)
pad2.SetFillStyle(4000) # will be transparent
pad1.Draw()
pad1.cd()
pad1.SetTicks(1,1)
pad1.SetLogz()
pad1.SetGridy()
pad1.SetGridx()
hxz.SetTitle(ptitle+" for #it{B} = "+btitle)
hxz.Draw("col")
dipoleExit.Draw("same")
dipoleOutXZL.Draw("same")
dipoleOutXZR.Draw("same")
dipoleInXZ.Draw("same")
for i in range(len(staveInnerLXZ)):
   staveInnerLXZ[i].Draw("same")
   staveInnerRXZ[i].Draw("same")
   staveOuterLXZ[i].Draw("same")
   staveOuterRXZ[i].Draw("same")
dipoleExit.SetLineColor(ROOT.kViolet)
dipoleExit.SetLineStyle(2)
pad1.Update()
cnv.cd(1)
# compute the pad range with suitable margins
ymin = hxE.GetYaxis().GetXmin()
ymax = hxE.GetYaxis().GetXmax()
dy = (ymax-ymin)/0.8 #10 margins top and bottom
xmin = hxE.GetXaxis().GetXmin()
xmax = hxE.GetXaxis().GetXmax()
dx = (xmax-xmin)/0.8; #10 per cent margins left and right
pad2.Range(xmin-0.1*dx,ymin-0.1*dy,xmax+0.1*dx,ymax+0.1*dy)
pad2.Draw()
pad2.cd()
hxE.SetMarkerColor(ROOT.kViolet)
hxE.Draw("][scat same")
pad2.Update()
# draw axis on the right side of the pad
axis = TGaxis(xmax,ymin,xmax,ymax,ymin,ymax,50510,"+L")
axis.SetLineColor(ROOT.kViolet)
axis.SetLabelColor(ROOT.kViolet)
axis.SetTitle("E [GeV]")
axis.SetTitleColor(ROOT.kViolet)
axis.SetLabelSize(0.035)
axis.SetTitleSize(0.035)
axis.Draw()
#### yz
cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogz()
ROOT.gPad.SetGridy()
ROOT.gPad.SetGridx()
# hyz = tfile.Get("h2_z_vs_y")
hyz.SetTitle(ptitle+" for #it{B} = "+btitle)
hyz.Draw("col")
dipoleOutYZ.Draw("same")
dipoleInYZ.Draw("same")
for i in range(len(staveYZ)):
   staveYZ[i].Draw("same")
cnv.SaveAs(fn+"xyz.pdf(")
cnv.SaveAs(fn+"xz_and_yz_and_xE.pdf")
cnv.SaveAs(fn+"xz_and_yz_and_xE.png")
cnv.Write()


cnv = TCanvas("cnv_xyz","",1000,1000)
cnv.Divide(2,1)
cnv.cd(1)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogz()
ROOT.gPad.SetGridy()
ROOT.gPad.SetGridx()
hxz.SetTitle(ptitle+" for #it{B} = "+btitle)
hxz.Draw("col")
dipoleOutXZL.Draw("same")
dipoleOutXZR.Draw("same")
dipoleInXZ.Draw("same")
for i in range(len(staveInnerLXZ)):
   staveInnerLXZ[i].Draw("same")
   staveInnerRXZ[i].Draw("same")
   staveOuterLXZ[i].Draw("same")
   staveOuterRXZ[i].Draw("same")
cnv.cd(2)
ROOT.gPad.SetTicks(1,1)
ROOT.gPad.SetLogz()
ROOT.gPad.SetGridy()
ROOT.gPad.SetGridx()
# hyz = tfile.Get("h2_z_vs_y")
hyz.SetTitle(ptitle+" for #it{B} = "+btitle)
hyz.Draw("col")
dipoleOutYZ.Draw("same")
dipoleInYZ.Draw("same")
for i in range(len(staveYZ)):
   staveYZ[i].Draw("same")
cnv.SaveAs(fn+"xyz.pdf")
cnv.SaveAs(fn+"xz_and_yz.pdf")
cnv.SaveAs(fn+"xz_and_yz.png")
cnv.Write()

cnv = TCanvas("cnv_xE","",1000,1000)
hxE.SetMarkerColor(ROOT.kViolet)
hxE.Draw("scat")
cnv.SaveAs(fn+"xyz.pdf)")
cnv.SaveAs(fn+"xE.pdf")
cnv.SaveAs(fn+"xE.png")
cnv.Write()

hxz.Write()
hyz.Write()
hxE.Write()

fOut.Write()
fOut.Close()