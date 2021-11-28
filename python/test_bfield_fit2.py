#!/usr/bin/python
import os
import math
import subprocess
import array
from array import array
import numpy as np
import ctypes
import ROOT
from ROOT import TFile, TH1D, TH2D, TF1, TCanvas, TPad, TGraph, TGraphAsymmErrors, TLegend, TLatex

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)

def label(text,x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1)
   s.SetTextAlign(13)
   s.SetTextColor(col)
   s.SetTextSize(0.044)
   if(boldit): text = "#font[72]{"+text+"}"
   s.DrawLatex(x,y,text)


def graph(name,h2,hq50,hq01,hq99):
   x,y     = array('d'), array('d')
   eyh,eyl = array('d'), array('d')
   exh,exl = array('d'), array('d')
   for b in range(1,hq50.GetNbinsX()+1):
      xc  = hq50.GetBinCenter(b)
      q50 = hq50.GetBinContent(b) 
      q99 = hq99.GetBinContent(b)
      q01 = hq01.GetBinContent(b)
      N = 0
      for by in range(1,h2.GetNbinsY()+1): N += h2.GetBinContent(b,by)
      if(N<2): continue
      x.append( xc )
      y.append( q50 )
      eyh.append( q99-q50 )
      eyl.append( q50-q01 )
      exh.append( 0 )
      exl.append( 0 )
   if(len(x)==0): return None
   gr = TGraphAsymmErrors(len(x),x,y,exl,exh,eyl,eyh)
   gr.SetName(name)
   gr.SetLineColor( ROOT.kBlack )
   gr.SetLineWidth( 1 )
   gr.SetMarkerColor( ROOT.kBlack )
   gr.SetMarkerStyle( 24 )
   gr.SetMarkerSize( 0.1 )
   return gr


def fit(message,name,graph,sfunc,xmin,xmax):
   print("\n"+message+":")
   f = TF1(name, sfunc, xmin,xmax)
   f.SetLineColor(ROOT.kRed)
   f.SetLineWidth(1)
   res = graph.Fit(f,"MERS")
   chi2dof = f.GetChisquare()/f.GetNDF() if(f.GetNDF()>0) else -1
   print(name+": chi2/Ndof=",chi2dof)
   return f


def band(name,h2,tf1,xmin,xmax,width=0.1):
   xup,yup = array('d'), array('d')
   xdn,ydn = array('d'), array('d')
   for b in range(1,h2.GetNbinsX()+1):
      xc = h2.GetXaxis().GetBinCenter(b)
      yc = tf1.Eval(xc)
      if(yc==0 or (xc<xmin or xc>xmax)): continue
      xup.append( xc )
      xdn.append( xc )
      # yup.append( yc*(1+width/(yc)) )
      # ydn.append( yc*(1-width/(yc)) )
      # yup.append( yc*(1+width) )
      # ydn.append( yc*(1-width) )
      yup.append( yc+width )
      ydn.append( yc-width )
   grUp = TGraph(len(xup),xup,yup)
   grDn = TGraph(len(xdn),xdn,ydn)
   grUp.SetName(name)
   grUp.SetLineColor( ROOT.kBlack )
   grUp.SetLineWidth( 1 )
   grDn.SetName(name)
   grDn.SetLineColor( ROOT.kBlack )
   grDn.SetLineWidth( 1 )
   return grUp,grDn
   

process = "glaser"


tfile = TFile("../data/root/dig/dig_"+process+"_.root","READ")

# hnames = ["h2_y_vs_x_exit",
#           "h2_y_vs_x_L1I","h2_y_vs_x_L1O","h2_y_vs_x_L4I","h2_y_vs_x_L4O",
#           "h2_E_vs_x_L1I","h2_E_vs_x_L1O","h2_E_vs_x_L4I","h2_E_vs_x_L4O",
#           "h2_dx14_vs_x_L4I","h2_dx14_vs_x_L4O" ]
hnames = ["h2_y_vs_x_exit",
          "h2_y_vs_x_L1I","h2_y_vs_x_L4I",
          "h2_E_vs_x_L1I","h2_E_vs_x_L4I",
          "h2_dx14_vs_x_L4I"]

fits = {"E_vs_x":"[0]+[1]/([2]+x)", "dx14_vs_x":"pol1"}
xmins = {"Eside":-50,"Pside":+5}
xmaxs = {"Eside":-5,"Pside":+50}
qMed = 50
qUp = 67
qDn = 33

histos    = {}
quantiles = {}
graphs    = {}
tf1s      = {}

### get the histos
for hname in hnames:
   name = hname
   h2 = tfile.Get(hname+"_Eside").Clone(hname)
   h2.Add(tfile.Get(hname+"_Pside"))
   histos.update( {name:h2} )
   histos[name].SetDirectory(0)
   ## get median quantile
   quantiles.update( {name+"_q"+str(qMed):histos[name].QuantilesX(qMed/100)} )
   quantiles[name+"_q"+str(qMed)].SetLineColor(ROOT.kBlack)
   quantiles[name+"_q"+str(qMed)].SetMarkerColor(ROOT.kBlack)
   quantiles[name+"_q"+str(qMed)].SetMarkerStyle(24)
   quantiles[name+"_q"+str(qMed)].SetMarkerSize(0.8)
   ## get 99% quantile
   quantiles.update({name+"_q"+str(qUp):histos[name].QuantilesX(qUp/100)})
   quantiles[name+"_q"+str(qUp)].SetLineColor(ROOT.kBlack)
   ## get 1% quantile
   quantiles.update({name+"_q"+str(qDn):histos[name].QuantilesX(qDn/100)})
   quantiles[name+"_q"+str(qDn)].SetLineColor(ROOT.kBlack)
   ## get the graph:
   gr = graph(name+"_gr", histos[name], quantiles[name+"_q"+str(qMed)], quantiles[name+"_q"+str(qDn)], quantiles[name+"_q"+str(qUp)])
   if(gr): graphs.update({name+"_gr":gr})


outname = "inputs_for_reco_"+process

### plot
cnv = TCanvas("c","",1200,500)
cnv.SaveAs(outname+".pdf(")
for name,h in histos.items():
   cnv = TCanvas("c_"+name,"",1200,500)
   cnv.Divide(2,1)
   cnv.cd(1)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetGrid()
   h.Draw("col")
   cnv.cd(2)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetGrid()
   h.Draw("col")
   if(name+"_gr" in graphs):
      graphs[name+"_gr"].Draw("ps")
      # if("y_vs_x" in name): continue
      sfunc = fits["E_vs_x"] if("E_vs_x" in name) else fits["dx14_vs_x"]
      tf1s.update({name+"_Eside":fit(name+" Eside",name+"_Eside",graphs[name+"_gr"],sfunc,xmins["Eside"],xmaxs["Eside"])})
      tf1s.update({name+"_Pside":fit(name+" Pside",name+"_Pside",graphs[name+"_gr"],sfunc,xmins["Pside"],xmaxs["Pside"])})
      tf1s[name+"_Eside"].Draw("same")
      tf1s[name+"_Pside"].Draw("same")
      bUp_Eside,bDn_Eside = band(name+"_band_Eside",h,tf1s[name+"_Eside"],xmins["Eside"],xmaxs["Eside"])
      bUp_Pside,bDn_Pside = band(name+"_band_Pside",h,tf1s[name+"_Pside"],xmins["Pside"],xmaxs["Pside"])
      bUp_Eside.Draw("l")
      bUp_Pside.Draw("l")
      bDn_Eside.Draw("l")
      bDn_Pside.Draw("l")
      
   cnv.SaveAs(outname+".pdf")
cnv = TCanvas("c","",1200,500)
cnv.SaveAs(outname+".pdf)")


fout = TFile(outname+".root","RECREATE")
fout.cd()
for name,h in histos.items():    h.Write()
for name,q in quantiles.items(): q.Write()
for name,g in graphs.items():    g.Write()
for name,f in tf1s.items():      f.Write()
fout.Write()
fout.Close()


# print("\nPositron side:")
# fEvsX_pos = TF1("fEvsX_pos", "[0]+[1]/([2]+x)", 1,20)
# fEvsX_pos.SetLineColor(ROOT.kRed)
# fEvsX_pos.SetLineWidth(1)
# res_pos = gr.Fit(fEvsX_pos,"MERS")
# chi2dof = fEvsX_pos.GetChisquare()/fEvsX_pos.GetNDF() if(fEvsX_pos.GetNDF()>0) else -1
# print("fEvsX_pos: chi2/Ndof=",chi2dof)
#
# print("\nElectron side:")
# fEvsX_ele = TF1("fEvsX_ele", "[0]+[1]/([2]+x)", -20,1)
# fEvsX_ele.SetLineColor(ROOT.kBlue)
# fEvsX_ele.SetLineWidth(1)
# res_ele = gr.Fit(fEvsX_ele,"MERS")
# chi2dof = fEvsX_ele.GetChisquare()/fEvsX_ele.GetNDF() if(fEvsX_ele.GetNDF()>0) else -1
# print("fEvsX_ele: chi2/Ndof=",chi2dof)
#
#
# leg1 = TLegend(0.55,0.8,0.87,0.87)
# # leg1.SetFillStyle(4000) # will be transparent
# leg1.SetFillColor(0)
# leg1.SetTextFont(42)
# leg1.SetBorderSize(0)
# hxE.SetFillColor(hxE.GetMarkerColor())
# hxE.SetLineColor(hxE.GetMarkerColor())
# leg1.AddEntry(hxE,"Simulation","F")
#
#
#
# leg2 = TLegend(0.55,0.6,0.87,0.87)
# # leg2.SetFillStyle(4000) # will be transparent
# leg2.SetFillColor(0)
# leg2.SetTextFont(42)
# leg2.SetBorderSize(0)
# leg2.AddEntry(gr,"Median^{+Q99%}_{-Q1%}","f")
# leg2.AddEntry(fEvsX_pos,"Fit (e^{+} side)","f")
# leg2.AddEntry(fEvsX_ele,"Fit (e^{-} side)","f")
#
# cnv = TCanvas("c1","",1200,500)
# cnv.Divide(2,1)
# cnv.cd(1)
# ROOT.gPad.SetTicks(1,1)
# ROOT.gPad.SetGrid()
# hxE.Draw("scat")
# leg1.Draw("same")
# label("B_{y}^{0} = -0.95 T",0.2,0.83)
# label("(non-uniform)",0.2,0.75)
# cnv.cd(2)
# ROOT.gPad.SetTicks(1,1)
# ROOT.gPad.SetGrid()
# hxE.Draw("scat")
# gr.Draw("ps")
# fEvsX_pos.Draw("same")
# fEvsX_ele.Draw("same")
# leg2.Draw("same")
# label("B_{y}^{0} = -0.95 T",0.2,0.83)
# label("(non-uniform)",0.2,0.75)
# cnv.SaveAs("test_bfield_fit.pdf")