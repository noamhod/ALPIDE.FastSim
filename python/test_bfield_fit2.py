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
ROOT.gStyle.SetPadTopMargin(0.05)
ROOT.gStyle.SetPadBottomMargin(0.11)
ROOT.gStyle.SetPadLeftMargin(0.12)
ROOT.gStyle.SetPadRightMargin(0.12)

def label(text,x,y,col=ROOT.kBlack,boldit=False):
   s = TLatex()
   s.SetNDC(1)
   s.SetTextFont(132)
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
   
   

process = "elaser"
isflat  = "_flat" ## "_flat" or ""

tfname = "../data/root/dig/"+process+"/phase0/gpc/5.0/dig_"+process+"_"+isflat+".root"
if("flat" in isflat): tfname = "../data/root/dig/"+process+"/flat/dig_"+process+"_"+isflat+".root"
tfile = TFile(tfname,"READ")

hnames = ["h2_y_vs_x_exit",
          "h2_y_vs_x_L1I","h2_y_vs_x_L4I","h2_y_vs_x_L1O","h2_y_vs_x_L4O",
          "h2_E_vs_x_L1I","h2_E_vs_x_L4I","h2_E_vs_x_L1O","h2_E_vs_x_L4O",
          "h2_dx14_vs_x_L4I","h2_dx14_vs_x_L4O","h2_dx14_vs_x_L4X",
          "h2_dy14_vs_y_L4I","h2_dy14_vs_y_L4O","h2_dy14_vs_y_L4X"]
   

# fits = {"E_vs_x":"[0]+[1]/([2]+x)", "dx14_vs_x":"pol1"}
fits = {"E_vs_x":"[0]/([1]+x)", "dx14_vs_x":"pol1"}
xmins = {"Eside":-60,"Pside":+4}
xmaxs = {"Eside":-4,"Pside":+60}
ymin = -0.688128
ymax = +0.811872
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
   h2 = tfile.Get(hname+"_Pside").Clone(hname)
   if(process=="glaser"): h2.Add(tfile.Get(hname+"_Eside"))
   histos.update( {name:h2} )
   histos[name].SetDirectory(0)
   # histos[name].SetTitle(hname.replace("h2_","").replace("_"," "))
   histos[name].SetTitle("")
   
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


outname = "inputs_for_reco_"+process+isflat

storage = os.getenv("$STORAGEDIR","../")
print("storage=",storage)
pdfsdir = storage+"/output/pdf/inputs_for_reco/"+process+"/"
if(isflat!=""): pdfsdir = pdfsdir+"flat/"
os.makedirs(pdfsdir, exist_ok=True)

### plot
cnv = TCanvas("c","",1200,500)
cnv.SaveAs(outname+".pdf(")
for name,h in histos.items():
   if(h.GetEntries()==0): continue
   text1 = "#gamma+laser" if(process=="glaser") else "e+laser"
   if(isflat!=""): text1 = text1+" (flat E)"
   
   if("y_vs_x" in name):
      cnv = TCanvas("c_"+name,"",700,500)
      ROOT.gPad.SetTicks(1,1)
      ROOT.gPad.SetGrid()
      h.Draw("col")
      text2 = name.replace("h2_","").replace("_"," ").replace("dx","#Deltax")
      label(text1,0.2,0.84)
      label(text2,0.2,0.79)
      cnv.RedrawAxis()
      cnv.SaveAs(outname+".pdf")
      cnv.SaveAs(pdfsdir+name.replace("h2_","")+".pdf")
   else:
      cnv = TCanvas("c_"+name,"",1200,500)
      cnv.Divide(2,1)
      cnv.cd(1)
      ROOT.gPad.SetTicks(1,1)
      ROOT.gPad.SetGrid()
      h.Draw("col")
      text2 = name.replace("h2_","").replace("_"," ").replace("dx","#Deltax").replace("dy","#Deltay").replace("14","")
      label(text1,0.2,0.84)
      label(text2,0.2,0.79)
      cnv.cd(2)
      ROOT.gPad.SetTicks(1,1)
      ROOT.gPad.SetGrid()
      h.Draw("col")
      if(name+"_gr" in graphs):
         graphs[name+"_gr"].Draw("ps")
         sfunc = fits["E_vs_x"] if("E_vs_x" in name) else fits["dx14_vs_x"]
         isdy = ("dy" in name)
         if(process=="glaser"): 
            if(not isdy): tf1s.update({name+"_Eside":fit(name+" Eside",name+"_Eside",graphs[name+"_gr"],sfunc,xmins["Eside"],xmaxs["Eside"])})
            else:         tf1s.update({name+"_Eside":fit(name+" Eside",name+"_Eside",graphs[name+"_gr"],sfunc,ymin,ymax)})
         if(not isdy): tf1s.update({name+"_Pside":fit(name+" Pside",name+"_Pside",graphs[name+"_gr"],sfunc,xmins["Pside"],xmaxs["Pside"])})
         else:         tf1s.update({name+"_Pside":fit(name+" Pside",name+"_Pside",graphs[name+"_gr"],sfunc,ymin,ymax)})
         if(process=="glaser"):
            if(not isdy): tf1s[name+"_Eside"].Draw("same")
            else:         tf1s[name+"_Eside"].Draw("same")
         if(not isdy): tf1s[name+"_Pside"].Draw("same")
         else:         tf1s[name+"_Pside"].Draw("same")
         # if("dx14_vs_x" in name):
         #    bUp_Eside,bDn_Eside = band(name+"_band_Eside",h,tf1s[name+"_Eside"],xmins["Eside"],xmaxs["Eside"]) if(process=="glaser") else 0,0
         #    bUp_Pside,bDn_Pside = band(name+"_band_Pside",h,tf1s[name+"_Pside"],xmins["Pside"],xmaxs["Pside"])
         #    if(process=="glaser"): bUp_Eside.Draw("l")
         #    bUp_Pside.Draw("l")
         #    if(process=="glaser"): bDn_Eside.Draw("l")
         #    bDn_Pside.Draw("l")
      leg = TLegend(0.60,0.73,0.87,0.87)
      leg.SetFillStyle(4000) # will be transparent
      leg.SetFillColor(0)
      leg.SetTextFont(42)
      leg.SetBorderSize(0)
      leg.AddEntry(graphs[name+"_gr"],"Median","pl")
      leg.AddEntry(tf1s[name+"_Pside"],"Fit","l")
      leg.Draw("same")
      cnv.RedrawAxis()
      cnv.SaveAs(outname+".pdf")
      cnv.SaveAs(pdfsdir+name.replace("h2_","")+".pdf")

cnv = TCanvas("c","",1200,500)
cnv.SaveAs(outname+".pdf)")


fout = TFile("../output/root/"+outname+".root","RECREATE")
fout.cd()
for name,h in histos.items():    h.Write()
for name,q in quantiles.items(): q.Write()
for name,g in graphs.items():    g.Write()
for name,f in tf1s.items():      f.Write()
fout.Write()
fout.Close()