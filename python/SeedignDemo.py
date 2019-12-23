#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)

#############################################
### electron mass:
me = 0.51099895/1000. ### GeV
me2 = me*me
cm2m = 1.e-2
cm2um = 1.e4
um2cm = 1.e-4

### magnetic field
B  = 1.4 # Tesla
LB = 1   # meters

### geometry:
zDipoleExit = 202.9
xDipoleExitMin = +1  ## cm
xDipoleExitMax = +30 ## cm
xDipoleExit = (xDipoleExitMax-xDipoleExitMin)/2.
yDipoleExitMin = -0.05 ## cm --> TODO: need tuning
yDipoleExitMax = +0.05 ## cm --> TODO: need tuning
xAbsMargins = 0.50 # cm --> TODO: need tuning
yAbsMargins = 0.05 # cm --> TODO: need tuning

### background stuff
NnoiseClusters = 1000

### histos
histos = {
   "h_dxrel":TH1D("h_dxrel",";(x_{cls}-x_{tru})/x_{tru};Tracks",200,-0.01,+0.01),
   "h_dyrel":TH1D("h_dyrel",";(y_{cls}-y_{tru})/y_{tru};Tracks",200,-0.25,+0.25),
}


#############################################
def GetFits():
   finname = "../output/root/fits_E_vs_x.root"
   tf = TFile(finname,"READ")
   fitsEx = {
              "Dipole":{"electrons":tf.Get("Ele_EofX0"), "positrons":tf.Get("Pos_EofX0")},
              "Layer1":{"electrons":tf.Get("Ele_EofX1"), "positrons":tf.Get("Pos_EofX1")},
              "Layer2":{"electrons":tf.Get("Ele_EofX2"), "positrons":tf.Get("Pos_EofX2")},
              "Layer3":{"electrons":tf.Get("Ele_EofX3"), "positrons":tf.Get("Pos_EofX3")},
              "Layer4":{"electrons":tf.Get("Ele_EofX4"), "positrons":tf.Get("Pos_EofX4")},
            }
   return fitsEx

def getgeometry(dipole=False):
   tfile = TFile("../data/root/geometry.root","READ")
   geometry = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
              tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
              tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
              tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2")]
   if(dipole): geometry.append(tfile.Get("TPolyLine3D;1"))
   return geometry

def PlaneStr(i):
   if  (i==0): return "Dipole"
   elif(i==1): return "Layer1"
   elif(i==2): return "Layer2"
   elif(i==3): return "Layer3"
   elif(i==4): return "Layer4"
   print("ERROR: unsupported plane index:",i)
   quit()
   return ""

def GetE(x,iplane,particles):
   if(iplane==0):
      if(abs(x)<1):
         print("ERROR: |x|<1 for iplane=0:",x)
         quit()
      if(abs(x)>30):
         print("ERROR: |x|>30 for iplane=0:",x)
         quit()
   else:   
      if(abs(x)<4):
         print("ERROR: |x|<4 for iplane>0:",x)
         quit()
      if(abs(x)>75):
         print("ERROR: |x|>75 for iplane>0:",x)
         quit()
   if(particles!="electrons" and particles!="positrons"):
      print("ERROR: particles string unsuppoerted:",particles)
      quit()
   if(particles=="electrons" and x<0):
      print("ERROR: electrons muxt come at positive x")
      quit()
   if(particles=="positrons" and x>0):
      print("ERROR: positrons muxt come at negative x")
      quit()
   E = fitsEx[PlaneStr(iplane)][particles].Eval(x)
   return E

def rUnit2(r1,r2):
   r = (r2-r1).Unit()
   return r ## TVector2

def rUnit3(r1,r2):
   r = (r2-r1).Unit()
   return r ## TVector3

def p3(r1,r2,E):
   r = rUnit3(r1,r2)
   p = TVector3()
   p.SetXYZ(E*r.X(),E*r.Y(),E*r.Z())
   return p ## TVector3

def p4(r1,r2,E):
   p = p3(r1,r2,E)
   tlv = TLorentzVector()
   tlv.SetXYZM(p.Px(),p.Py(),p.Pz(),me)
   return tlv ## TLorentzVector

def setnoiseclusters(noiseclusters,N):
   rnd = TRandom()
   rnd.SetSeed()
   for i in range(N):
      x = rnd.Uniform(4,31)
      y = rnd.Uniform(-0.75,+0.75)
      z = rnd.Integer(4)
      if(z==0): z = 300
      if(z==1): z = 310
      if(z==2): z = 320
      if(z==3): z = 330
      noiseclusters.SetNextPoint(x,y,z)

def xofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   a = dx/dz
   b = r1[0]-a*r1[2]
   x = a*z+b
   return x

def yofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dy = r2[1]-r1[1]
   a = dy/dz
   b = r1[1]-a*r1[2]
   y = a*z+b
   return y
   
def zofx(r1,r2,x):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   a = dz/dx
   b = r1[2]-a*r1[0]
   z = a*x+b
   return z

def countpoints(points3d):
   npoints = 0
   for i in range(points3d.GetN()):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      points3d.GetPoint(i,r[0],r[1],r[2])
      if(abs(r[0])<xDipoleExitMin): continue
      if(r[2]<300 or r[2]>330):     continue
      npoints += 1
   return npoints

def getpoint(pm, i):
   r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   pm.GetPoint(i,r[0],r[1],r[2])
   return r[0],r[1],r[2]

def windowline(points,isclosed=True):
   line = TPolyLine3D()
   arr = []
   for i in range(points.GetN()):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      points.GetPoint(i,r[0],r[1],r[2])
      if(r[0]==0): continue
      arr.append([r[0],r[1],r[2]])      
   if(isclosed):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      points.GetPoint(0,r[0],r[1],r[2])
      arr.append([r[0],r[1],r[2]])
   for point in arr:
      line.SetNextPoint(point[0],point[1],point[2])
   return line

def setyzwindow(window,cluster):
   x1,y1,z1 = getpoint(cluster,0)
   y1min = y1-yAbsMargins if((y1-yAbsMargins)>-0.75) else -0.75 ## must be within stave volume
   y1max = y1+yAbsMargins if((y1+yAbsMargins)<+0.75) else +0.75 ## must be within stave volume
   window.SetNextPoint(x1,y1min,z1)
   window.SetNextPoint(x1,y1max,z1)
   window.SetNextPoint(xDipoleExit,yDipoleExitMax,zDipoleExit)
   window.SetNextPoint(xDipoleExit,yDipoleExitMin,zDipoleExit)
   wondow_line = windowline(window)
   wondow_line.SetLineWidth(1)
   wondow_line.SetLineColor(ROOT.kBlack)
   return wondow_line

def setwidewindow(window,cluster):
   x1,y1,z1 = getpoint(cluster,0)
   x1min = x1-xAbsMargins if((x1-xAbsMargins)>4)  else 4  ## must be within stave volume
   x1max = x1+xAbsMargins if((x1+xAbsMargins)<31) else 31 ## must be within stave volume
   window.SetNextPoint(x1min,y1,z1)                   ## top left
   window.SetNextPoint(x1max,y1,z1)                   ## top right
   window.SetNextPoint(xDipoleExitMax,y1,zDipoleExit) ## bottom right
   window.SetNextPoint(xDipoleExitMin,y1,zDipoleExit) ## bottom left
   wondow_line = windowline(window)
   wondow_line.SetLineWidth(1)
   wondow_line.SetLineColor(ROOT.kBlack)
   return wondow_line

def setnarrwindow(window,cluster1,cluster2):
   x1,y1,z1 = getpoint(cluster1,0)   
   x2,y2,z2 = getpoint(cluster2,0)   
   window.SetNextPoint(x1-xAbsMargins,y1,z1) ## top left
   window.SetNextPoint(x1+xAbsMargins,y1,z1) ## top right
   window.SetNextPoint(x2+xAbsMargins,y2,z2) ## bottom right
   window.SetNextPoint(x2-xAbsMargins,y2,z2) ## bottom left
   wondow_line = windowline(window)
   wondow_line.SetLineWidth(1)
   wondow_line.SetLineColor(ROOT.kBlue)
   return wondow_line

def getxminmax(window,zTest):
   windowarr = []
   for i in range(window.GetN()):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      window.GetPoint(i,r[0],r[1],r[2])
      if(r[0]==0): continue
      windowarr.append([r[0],r[1],r[2]])
   ### get the equation of the 
   xmin = xofz(windowarr[0],windowarr[3],zTest)
   xmax = xofz(windowarr[1],windowarr[2],zTest)
   return xmin,xmax
   
def getyminmax(window,zTest):
   windowarr = []
   for i in range(window.GetN()):
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      window.GetPoint(i,r[0],r[1],r[2])
      if(r[0]==0): continue
      windowarr.append([r[0],r[1],r[2]])
   ### get the equation of the 
   ymin = yofz(windowarr[0],windowarr[3],zTest)
   ymax = yofz(windowarr[1],windowarr[2],zTest)
   return ymin,ymax
   
def isinwindowx(points,i,window):
   ### get the test point
   rTest = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   points.GetPoint(i,rTest[0],rTest[1],rTest[2])
   if(rTest[0]<1): return False
   ### get the window x extremes at zTest
   xmin,xmax = getxminmax(window,rTest[2])
   if(rTest[0]>=xmin and rTest[0]<=xmax): return True
   return False
   
def isinwindowyz(points,i,window):
   ### get the test point
   rTest = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   points.GetPoint(i,rTest[0],rTest[1],rTest[2])
   if(rTest[0]<1): return False
   ### get the window x extremes at zTest
   ymin,ymax = getyminmax(window,rTest[2])
   # print("z=%g: y=%g --> ymin=%g, ymax=%g" % (rTest[2],rTest[1],ymin,ymax))
   if(rTest[1]>=ymin and rTest[1]<=ymax): return True
   return False
   
def trimwide(clusters,clusters_wide,windowx,windowyz):
   for i in range(clusters.GetN()):
      acceptx  = isinwindowx(clusters,i,windowx)
      acceptyz = isinwindowyz(clusters,i,windowyz)
      if(not acceptx):  continue
      if(not acceptyz): continue
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      clusters.GetPoint(i,r[0],r[1],r[2])
      clusters_wide.SetNextPoint(r[0],r[1],r[2])
      
def trimsigwide(clusters_sig,clusters_sig_wide,windowx,windowyz):
   for i in range(clusters_sig.GetN()):
      acceptx  = isinwindowx(clusters_sig,i,windowx)
      acceptyz = isinwindowyz(clusters_sig,i,windowyz)
      if(not acceptx):  continue
      if(not acceptyz): continue
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      clusters_sig.GetPoint(i,r[0],r[1],r[2])
      clusters_sig_wide.SetNextPoint(r[0],r[1],r[2])

def trimnarr(clusters_wide,clusters_narr,window):
   for i in range(clusters_wide.GetN()):
      accept = isinwindowx(clusters_wide,i,window)
      if(not accept): continue
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      clusters_wide.GetPoint(i,r[0],r[1],r[2])
      clusters_narr.SetNextPoint(r[0],r[1],r[2])
      
def trimsignarr(clusters_sig_wide,clusters_sig_narr,window):
   for i in range(clusters_sig_wide.GetN()):
      accept = isinwindowx(clusters_sig_wide,i,window)
      if(not accept): continue
      r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
      clusters_sig_wide.GetPoint(i,r[0],r[1],r[2])
      clusters_sig_narr.SetNextPoint(r[0],r[1],r[2])

def makeseed(cluster1,cluster2,particles):
   p = TLorentzVector()
   r1 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   r2 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   cluster1.GetPoint(0,r1[0],r1[1],r1[2])
   cluster2.GetPoint(0,r2[0],r2[1],r2[2])

   # E1 = GetE(r1[0],4,"electrons") ## energy of the first seed cluster
   # E2 = GetE(r2[0],1,"electrons") ## energy of the second seed cluster
   # E = (E1+E2)/2.
   
   x0 = 0
   z0 = zofx(r1,r2,x0)
   xExit = xofz(r1,r2,zDipoleExit)*cm2m
   q = 1 if(particles=="electrons") else -1
   H = (zDipoleExit-z0)*cm2m
   R = H*(LB)/xExit + xExit ## look this up in my sides
   P = 0.3*B*R
   
   v1 = TVector2(r1[2],r1[1])
   v2 = TVector2(r2[2],r2[1])
   u = rUnit2(v2,v1)
   uz = u.X()
   uy = u.Y()
   px = 0
   py = P*uy
   pz = P*uz
   p.SetPxPyPzE(px,py,pz,math.sqrt(P*P+me2))
   print("Seed: E=%g, pT=%g, eta=%g, phi=%g, theta=%g" % (p.E(),p.Pt(),p.Eta(),p.Phi(),p.Theta()) )
   
   return p


#############################################
### get tracks and clusters
def Read(event,points):
   for i in range(event.polm_clusters.size()): ## clusters' vectors are always written out (even if empty) for all gen tracks!
      wgt  = event.wgtgen[i]
      pgen = event.pgen[i]
      ### cut on acceptance
      if(event.acctrkgen[i]!=1): continue 
      if(i==0):
         print("Genr: P=%g, E=%g, pT=%g, eta=%g, phi=%g, theta=%g, rapidity=%g" % (pgen.P(),pgen.E(),pgen.Pt(),pgen.Eta(),pgen.Phi(),pgen.Theta(), pgen.Rapidity() ) )
         print("Genr: P=%g, E=%g, px=%g, py=%g, pz=%g" % (pgen.P(),pgen.E(),pgen.Px(),pgen.Py(),pgen.Pz() ) )
      ### loop over clusters along track (simulated)
      for jxy in range(event.polm_clusters[i].GetN()):
         rgen = [ ROOT.Double(), ROOT.Double(), ROOT.Double() ]
         rcls = [ ROOT.Double(), ROOT.Double(), ROOT.Double() ]
         event.polm_clusters[i].GetPoint(jxy,rcls[0],rcls[1],rcls[2]) ### the clusters
         if(rcls[2]<300): continue
         for kxy in range(event.polm_gen[i].GetN()):
            event.polm_gen[i].GetPoint(kxy,rgen[0],rgen[1],rgen[2]) ### the truth track
            if(rgen[2]==rcls[2]): break
         # print("rgen=",rgen)
         # print("rcls=",rcls)
         if(rcls[2]==330): points["cluster_seed1"].SetNextPoint(rcls[0],rcls[1],rcls[2])
         if(rcls[2]==320): points["cluster_layr3"].SetNextPoint(rcls[0],rcls[1],rcls[2])
         if(rcls[2]==310): points["cluster_layr2"].SetNextPoint(rcls[0],rcls[1],rcls[2])
         if(rcls[2]==300): points["cluster_seed2"].SetNextPoint(rcls[0],rcls[1],rcls[2])
         dxrel = (rcls[0]-rgen[0])/rgen[0]
         dyrel = (rcls[1]-rgen[1])/rgen[1]
         # print("dxrel=%g, dyrel=%g" % (dxrel,dyrel) )
         histos["h_dxrel"].Fill(dxrel)
         histos["h_dyrel"].Fill(dyrel)

## event loop
def EventLoop(tree,nmax=0):
   nevents = tree.GetEntries()
   print("with %d events" % nevents)
   n=0 ### init n
   for event in tree:
      if(n>nmax): break
      if(n%10==0 and n>0): print("  processed %d events" % n)
      Read(event,points)
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

## drawing
def draw(name,title,points,suffix):
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.ShowAxis()
   # view.SetRange(-80,-50,0, +80,+50,350)
   view.SetRange(0,-10,190, +30,+10,340)
   geom = getgeometry() ### get the geometry from the root file
   for g in geom: g.Draw()
   window_wide.Draw("same")
   window_narr.Draw("same")
   window_yz.Draw("same")
   points["cluster_seed1"+suffix].Draw("same")
   points["cluster_layr3"+suffix].Draw("same")
   points["cluster_layr2"+suffix].Draw("same")
   points["cluster_seed2"+suffix].Draw("same")
   print(title+' --> N clusters:',countpoints(points["cluster_noise"+suffix]))
   points["cluster_noise"+suffix].Draw("same")
   cnv.SaveAs(name)
   rootname = name.replace(".pdf","."+title+".root").replace("(","").replace(")","")
   cnv.SaveAs(rootname)

def draw2(name,h1,h2):
   cnv = TCanvas("","",1000,500)
   cnv.Divide(2,1)
   cnv.cd(1)
   h1.Draw()
   cnv.cd(2)
   h2.Draw()
   cnv.SaveAs(name)

##########################################################################################
##########################################################################################
##########################################################################################
### define the data structures to hold clusters and windows
points = { "cluster_seed1" : TPolyMarker3D(),
           "cluster_seed2" : TPolyMarker3D(),
           "cluster_layr2" : TPolyMarker3D(),
           "cluster_layr3" : TPolyMarker3D(),
           
           "cluster_seed1_wide" : TPolyMarker3D(),
           "cluster_seed2_wide" : TPolyMarker3D(),
           "cluster_layr2_wide" : TPolyMarker3D(),
           "cluster_layr3_wide" : TPolyMarker3D(),
           
           "cluster_seed1_narr" : TPolyMarker3D(),
           "cluster_seed2_narr" : TPolyMarker3D(),
           "cluster_layr2_narr" : TPolyMarker3D(),
           "cluster_layr3_narr" : TPolyMarker3D(),
           
           "cluster_noise"      : TPolyMarker3D(), ## about 100 tracks per event per side
           "cluster_noise_wide" : TPolyMarker3D(), ## after removing the points not in the wide window
           "cluster_noise_narr" : TPolyMarker3D(), ## after removing the points not in the narrow window
           
           "window_wide"   : TPolyMarker3D(),
           "window_narr"   : TPolyMarker3D(),
           "window_yz"     : TPolyMarker3D(),
}

### read fits to root file
fitsEx = GetFits()

### set the signal clusters
nevents = Run("../data/root/rec_bppp.root","res")
for name,marker3d in points.items():
   if("cluster_seed" in name or "cluster_layr" in name): marker3d.SetMarkerColor(ROOT.kRed)

### initialise the noise clusters with 100 clusters
setnoiseclusters(points["cluster_noise"],NnoiseClusters)

### set the yz window
window_yz = setyzwindow(points["window_yz"],points["cluster_seed1"])

### set the wide window starting from cluster_seed1 (window corners must be added clockwise!)
window_wide = setwidewindow(points["window_wide"],points["cluster_seed1"])

### discard all clusters which are not in the wide window
trimwide(points["cluster_noise"],points["cluster_noise_wide"],points["window_wide"],points["window_yz"])
trimsigwide(points["cluster_seed1"],points["cluster_seed1_wide"],points["window_wide"],points["window_yz"])
trimsigwide(points["cluster_layr3"],points["cluster_layr3_wide"],points["window_wide"],points["window_yz"])
trimsigwide(points["cluster_layr2"],points["cluster_layr2_wide"],points["window_wide"],points["window_yz"])
trimsigwide(points["cluster_seed2"],points["cluster_seed2_wide"],points["window_wide"],points["window_yz"])

### TODO: choose one cluster in layer 1 as the second seed
##

### set the narrow window  (window corners must be added clockwise!)
window_narr = setnarrwindow(points["window_narr"],points["cluster_seed1"],points["cluster_seed2"])

### discard all clusters which are not in the narrow window
trimnarr(points["cluster_noise_wide"],points["cluster_noise_narr"],points["window_narr"])
trimsignarr(points["cluster_seed1_wide"],points["cluster_seed1_narr"],points["window_narr"])
trimsignarr(points["cluster_layr3_wide"],points["cluster_layr3_narr"],points["window_narr"])
trimsignarr(points["cluster_layr2_wide"],points["cluster_layr2_narr"],points["window_narr"])
trimsignarr(points["cluster_seed2_wide"],points["cluster_seed2_narr"],points["window_narr"])

### TODO: check if there are at least 1 cluster in both layer 2 and layer 3 within the narrow window
##

### get the seed 4-momentum
pseed = makeseed(points["cluster_seed1"],points["cluster_seed2"],"electrons")

### TODO: propagate the seed to the IP via the dipole
##

### TODO: give the seed kinematics and the 4 clusters to the KF and perform the fit 
## 

### TODO: redo last step for all other clusters in layer 2 and layer 4
##

### draw
pdfname = "../output/pdf/seedingdemo.pdf"
draw(pdfname+"(", "AllNoiseClusters",points,"")
draw(pdfname,     "WideWinNoiseClusters",points,"_wide")
draw(pdfname+")", "NarrWinNoiseClusters",points,"_narr")
pdfname = "../output/pdf/seedingdemo_clusters.pdf"
draw2(pdfname,histos["h_dxrel"],histos["h_dyrel"])