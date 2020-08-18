#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TF1, TF2, TF3, TCanvas, TPolyLine

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)

zActiveDipoleLengh = 102.9 #cm
zActiveDipoleEntrance = 100 #cm
zActiveDipoleExit = 202.9 #cm
xActiveDipoleWidth = 33 #cm
yActiveDipoleHeight = 10.8 #cm
BfieldGamLaser = 1.6 #T
BfieldEleLaser = 1 #T

strength = BfieldGamLaser

x_xy = array.array('d', [-xActiveDipoleWidth/2,-xActiveDipoleWidth/2,+xActiveDipoleWidth/2,+xActiveDipoleWidth/2,-xActiveDipoleWidth/2])
y_xy = array.array('d', [-yActiveDipoleHeight/2,+yActiveDipoleHeight/2,+yActiveDipoleHeight/2,-yActiveDipoleHeight/2,-yActiveDipoleHeight/2])
n = len(x_xy)
xy = TPolyLine(n,x_xy,y_xy)
xy.SetLineColor(ROOT.kBlack)

x_xz = array.array('d', [-xActiveDipoleWidth/2,-xActiveDipoleWidth/2,+xActiveDipoleWidth/2,+xActiveDipoleWidth/2,-xActiveDipoleWidth/2])
y_xz = array.array('d', [-zActiveDipoleLengh/2,+zActiveDipoleLengh/2,+zActiveDipoleLengh/2,-zActiveDipoleLengh/2,-zActiveDipoleLengh/2])
n = len(x_xz)
xz = TPolyLine(n,x_xz,y_xz)
xz.SetLineColor(ROOT.kBlack)

x_yz = array.array('d', [-yActiveDipoleHeight/2,-yActiveDipoleHeight/2,+yActiveDipoleHeight/2,+yActiveDipoleHeight/2,-yActiveDipoleHeight/2])
y_yz = array.array('d', [-zActiveDipoleLengh/2,+zActiveDipoleLengh/2,+zActiveDipoleLengh/2,-zActiveDipoleLengh/2,-zActiveDipoleLengh/2])
n = len(x_yz)
yz = TPolyLine(n,x_yz,y_yz)
yz.SetLineColor(ROOT.kBlack)


x = array.array('d', [-xActiveDipoleWidth/2,-xActiveDipoleWidth/2])
y = array.array('d', [0,strength])
n = len(x)
Rx = TPolyLine(n,x,y)
Rx.SetLineColor(ROOT.kBlack)
x = array.array('d', [+xActiveDipoleWidth/2,+xActiveDipoleWidth/2])
y = array.array('d', [0,strength])
n = len(x)
Lx = TPolyLine(n,x,y)
Lx.SetLineColor(ROOT.kBlack)

x = array.array('d', [-yActiveDipoleHeight/2,-yActiveDipoleHeight/2])
y = array.array('d', [0,strength])
n = len(x)
Ry = TPolyLine(n,x,y)
Ry.SetLineColor(ROOT.kBlack)
x = array.array('d', [+yActiveDipoleHeight/2,+yActiveDipoleHeight/2])
y = array.array('d', [0,strength])
n = len(x)
Ly = TPolyLine(n,x,y)
Ly.SetLineColor(ROOT.kBlack)

x = array.array('d', [-zActiveDipoleLengh/2,-zActiveDipoleLengh/2])
y = array.array('d', [0,strength])
n = len(x)
Rz = TPolyLine(n,x,y)
Rz.SetLineColor(ROOT.kBlack)
x = array.array('d', [+zActiveDipoleLengh/2,+zActiveDipoleLengh/2])
y = array.array('d', [0,strength])
n = len(x)
Lz = TPolyLine(n,x,y)
Lz.SetLineColor(ROOT.kBlack)



s3func = '[0] * ((2-ROOT::Math::erfc(([1]/2+x)/[2]))/2)*((2-ROOT::Math::erfc(([1]/2-x)/[2]))/2)\
              * ((2-ROOT::Math::erfc(([3]/2+y)/[4]))/2)*((2-ROOT::Math::erfc(([3]/2-y)/[4]))/2)\
              * ((2-ROOT::Math::erfc(([5]/2+z)/[6]))/2)*((2-ROOT::Math::erfc(([5]/2-z)/[6]))/2)'
fieldxyz = TF3("fieldxyz", s3func, -25,+25, -10,+10, -80,+80)


s2func = '[0] * ((2-ROOT::Math::erfc(([1]/2+x)/[2]))/2) * ((2-ROOT::Math::erfc(([1]/2-x)/[2]))/2)\
              * ((2-ROOT::Math::erfc(([3]/2+y)/[4]))/2) * ((2-ROOT::Math::erfc(([3]/2-y)/[4]))/2)'
fieldxy = TF2("fieldxy", s2func, -25,+25, -10,+10)
fieldxz = TF2("fieldxz", s2func, -25,+25, -80,+80)
fieldyz = TF2("fieldyz", s2func, -10,+10, -80,+80)

s1func = '[0] * ((2-ROOT::Math::erfc(([1]/2+x)/[2]))/2) * ((2-ROOT::Math::erfc(([1]/2-x)/[2]))/2)'
fieldx = TF1("fieldx", s1func, -25,+25)
fieldy = TF1("fieldy", s1func, -10,+10)
fieldz = TF1("fieldz", s1func, -80,+80)

fieldxy.SetTitle("B-field x:y")
fieldxz.SetTitle("B-field x:z")
fieldyz.SetTitle("B-field y:z")

fieldx.SetTitle("B-field x")
fieldy.SetTitle("B-field y")
fieldz.SetTitle("B-field z")


fieldxy.SetParameter('p0', strength)
fieldxy.SetParameter('p1', xActiveDipoleWidth)
fieldxy.SetParameter('p2', 4)
fieldxy.SetParameter('p3', yActiveDipoleHeight)
fieldxy.SetParameter('p4', 2)

fieldxz.SetParameter('p0', strength)
fieldxz.SetParameter('p1', xActiveDipoleWidth)
fieldxz.SetParameter('p2', 4)
fieldxz.SetParameter('p3', zActiveDipoleLengh)
fieldxz.SetParameter('p4', 16)

fieldyz.SetParameter('p0', strength)
fieldyz.SetParameter('p1', yActiveDipoleHeight)
fieldyz.SetParameter('p2', 2)
fieldyz.SetParameter('p3', zActiveDipoleLengh)
fieldyz.SetParameter('p4', 16)



fieldxyz.SetParameter('p0', strength)
fieldxyz.SetParameter('p1', xActiveDipoleWidth)
fieldxyz.SetParameter('p2', 4)
fieldxyz.SetParameter('p3', yActiveDipoleHeight)
fieldxyz.SetParameter('p4', 2)
fieldxyz.SetParameter('p5', zActiveDipoleLengh)
fieldxyz.SetParameter('p6', 16)



fieldxy.SetNpy(1000)
fieldxz.SetNpy(1000)
fieldyz.SetNpy(1000)


fieldx.SetParameter('p0', strength)
fieldx.SetParameter('p1', xActiveDipoleWidth)
fieldx.SetParameter('p2', 4)

fieldy.SetParameter('p0', strength)
fieldy.SetParameter('p1', yActiveDipoleHeight)
fieldy.SetParameter('p2', 2)

fieldz.SetParameter('p0', strength)
fieldz.SetParameter('p1', zActiveDipoleLengh)
fieldz.SetParameter('p2', 16)


cnv = TCanvas("cnv","",1200,400)
cnv.Divide(3,1)
cnv.cd(1)
fieldx.SetMaximum(strength*1.1)
fieldx.SetMinimum(0)
fieldx.Draw()
Rx.Draw("same")
Lx.Draw("same")
cnv.cd(2)
fieldy.SetMaximum(strength*1.1)
fieldy.SetMinimum(0)
fieldy.Draw()
Ry.Draw("same")
Ly.Draw("same")
cnv.cd(3)
fieldz.SetMaximum(strength*1.1)
fieldz.SetMinimum(0)
fieldz.Draw()
Rz.Draw("same")
Lz.Draw("same")
cnv.SaveAs("BfieldMap_1d.png")
cnv.SaveAs("BfieldMap.pdf(")



cnv = TCanvas("cnv","",1200,400)
cnv.Divide(3,1)
cnv.cd(1)
# fieldxy.Draw('surf1')
fieldxy.SetMaximum(strength)
fieldxy.SetMinimum(0)
fieldxy.Draw("colz")
xy.Draw("same")
cnv.cd(2)
# fieldxz.Draw('surf1')
fieldxz.SetMaximum(strength)
fieldxz.SetMinimum(0)
fieldxz.Draw("colz")
xz.Draw("same")
cnv.cd(3)
# fieldyz.Draw('surf1')
fieldyz.SetMaximum(strength)
fieldyz.SetMinimum(0)
fieldyz.Draw("colz")
yz.Draw("same")
cnv.SaveAs("BfieldMap_2d.png")
cnv.SaveAs("BfieldMap.pdf")


cnv = TCanvas("cnv","",500,500)
fieldxyz.SetClippingBoxOn(0,0,0)
fieldxyz.SetFillColor(30)
fieldxyz.SetLineColor(15)
fieldxyz.Draw('GLS')
cnv.SaveAs("BfieldMap_3d.png")
cnv.SaveAs("BfieldMap.pdf)")

print("field at r=(0.0,0.0,50.0): %.3f [T]" % fieldxyz.Eval(0,0,50))
print("field at r=(0.5,8.0,60.0): %.3f [T]" % fieldxyz.Eval(0.5,8,60))
print("field at r=(0.0,0.0,0.0) : %.3f [T]" % fieldxyz.Eval(0.0,0,0))

