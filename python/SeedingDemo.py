#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import config as cfg
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend
import argparse
parser = argparse.ArgumentParser(description='SeedingDemo.py...')
parser.add_argument('-p', metavar='process', required=True,  help='physics process [trident or bppp]')
parser.add_argument('-s', metavar='sides',   required=False, help='detector side [e+, e-, e+e-]')
parser.add_argument('-e', metavar='energy',  required=False, help='beam energy')
argus = parser.parse_args()
proc  = argus.p
sides = "e+" if(proc=="trident") else "e+e-" ## jsut the default
if(argus.s is not None): sides = argus.s
if(proc=="trident" and "e-" in sides):
   print("ERROR: do not run tracking in the electron side for trident")
   quit()
print("Running with proc=%s and sides=%s" % (proc,sides))

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError

#############################################
### read configuration from file
cfgmap = cfg.set(proc,sides,True)

#############################################
NeventsToDraw = 10

### electron mass:
me = 0.51099895/1000. ### GeV
me2 = me*me
cm2m = 1.e-2
cm2um = 1.e4
um2cm = 1.e-4

### magnetic field
B  = 1.4 if(proc=="trident") else 2.0 # Tesla
LB = 1   # meters

### possible energies 
Emax = 17.5 if(proc=="trident") else 16 # GeV
Emin = 1.00 if(proc=="trident") else 2 # GeV

### geometry:
zDipoleExit = 202.9
xDipoleExitMinAbs = 1.5 if(proc=="bppp") else 4   ## cm --> TODO: need tuning
xDipoleExitMaxAbs = 25  if(proc=="bppp") else 30  ## cm --> TODO: need tuning
yDipoleExitMin = -0.05 ## cm --> TODO: need tuning
yDipoleExitMax = +0.05 ## cm --> TODO: need tuning
xAbsMargins = 0.025 # cm --> TODO: need tuning
yAbsMargins = 0.025 if(proc=="bppp") else 0.1 # cm --> TODO: need tuning

### stave geometry
Hstave    = 1.5  # cm
Lstave    = 27 if(proc=="bppp") else 50 # cm
Rbeampipe = 2.413 # 4 # cm
RoffsetBfield22BPPP = 7.0  # cm for BPPP in B=2.2T
RoffsetBfield20BPPP = 5.7  # cm for BPPP in B=2.0T
RoffsetBfield14BPPP = 4.0  # cm for BPPP in B=1.4T
RoffsetBfield = RoffsetBfield20BPPP if(proc=="bppp") else 14 # cm
xPsideL = -RoffsetBfield-Lstave
xPsideR = -RoffsetBfield       
xEsideL = +RoffsetBfield       
xEsideR = +RoffsetBfield+Lstave
yUp = +Hstave/2.
yDn = -Hstave/2.

### for the histogram
detXmin = xPsideL
detXmax = xEsideR
if(proc=="trident"): detXmax = xPsideR
if(proc=="bppp" and sides=="e+"): detXmax = xPsideR
if(proc=="bppp" and sides=="e-"): detXmin = xEsideL



### background stuff
NnoiseClusters = 50  ## uniformly distributed clusters in x:y for each layer
NbkgTracks     = 200 ## uniformly distributed tracks coming from the inner part of the beampipe at the dipole-exit plane

layers = [1,2,3,4]

### histos
histos = {
   "h_dxrel":TH1D("h_dxrel",";(x_{cls}-x_{tru})/x_{tru};Tracks",200,-0.01,+0.01),
   "h_dyrel":TH1D("h_dyrel",";(y_{cls}-y_{tru})/y_{tru};Tracks",200,-0.25,+0.25),
}


#############################################
def getgeometry(dipole=False):
   tfile = TFile("../data/root/"+proc+"_geometry.root","READ")
   geometry = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
              tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
              tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
              tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2")]
   if(dipole): geometry.append(tfile.Get("TPolyLine3D;1"))
   return geometry

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

def xofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   if(dx==0):
      print("ERROR in xofz: dx=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
   a = dx/dz
   b = r1[0]-a*r1[2]
   x = a*z+b
   return x

def yofz(r1,r2,z):
   dz = r2[2]-r1[2]
   dy = r2[1]-r1[1]
   if(dz==0):
      print("ERROR in yofz: dz=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
   a = dy/dz
   b = r1[1]-a*r1[2]
   y = a*z+b
   return y
   
def zofx(r1,r2,x):
   dz = r2[2]-r1[2]
   dx = r2[0]-r1[0]
   if(dx==0):
      print("ERROR in zofx: dx=0 --> r1[0]=%g,r2[0]=%g, r1[1]=%g,r2[1]=%g, r1[2]=%g,r2[2]=%g" % (r1[0],r2[0],r1[1],r2[1],r1[2],r2[2]))
      quit()
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

def getyzwindow(cluster,i):
   x1,y1,z1 = getpoint(cluster,i)
   particles = "electrons" if(x1>0) else "positrons"
   y1min = y1-yAbsMargins if((y1-yAbsMargins)>-0.75) else -0.75 ## must be within stave volume
   y1max = y1+yAbsMargins if((y1+yAbsMargins)<+0.75) else +0.75 ## must be within stave volume
   window_pts = TPolyMarker3D()
   window_pts.SetNextPoint(x1,y1min,z1)
   window_pts.SetNextPoint(x1,y1max,z1)
   xmiddle = +xDipoleExitMinAbs+(abs(x1)*0.99-xDipoleExitMinAbs)/2
   if(particles=="electrons"):
      window_pts.SetNextPoint(+xmiddle,yDipoleExitMax,zDipoleExit)
      window_pts.SetNextPoint(+xmiddle,yDipoleExitMin,zDipoleExit)
   else:
      window_pts.SetNextPoint(-xmiddle,yDipoleExitMax,zDipoleExit)
      window_pts.SetNextPoint(-xmiddle,yDipoleExitMin,zDipoleExit)
   window_lin = windowline(window_pts)
   window_lin.SetLineWidth(1)
   window_lin.SetLineColor(ROOT.kBlack)
   return window_pts,window_lin

def getwidewindow(cluster,i):
   x1,y1,z1 = getpoint(cluster,i)
   particles = "electrons" if(x1>0) else "positrons"
   x1min = 0
   x1max = 0
   window_pts = TPolyMarker3D()
   if(particles=="electrons"): 
      x1min = x1-xAbsMargins if((x1-xAbsMargins)>xEsideL) else xEsideL ## must be within stave volume
      x1max = x1+xAbsMargins if((x1+xAbsMargins)<xEsideR) else xEsideR ## must be within stave volume
      window_pts.SetNextPoint(x1min,y1,z1) ## top left
      window_pts.SetNextPoint(x1max,y1,z1) ## top right
      # window_pts.SetNextPoint(xDipoleExitMaxAbs,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(x1max*0.99,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(xDipoleExitMinAbs,y1,zDipoleExit) ## bottom left
   else:
      x1min = x1-xAbsMargins if((x1-xAbsMargins)>xPsideL) else xPsideL ## must be within stave volume
      x1max = x1+xAbsMargins if((x1+xAbsMargins)<xPsideR) else xPsideR ## must be within stave volume
      window_pts.SetNextPoint(x1min,y1,z1) ## top left
      window_pts.SetNextPoint(x1max,y1,z1) ## top right
      # window_pts.SetNextPoint(-xDipoleExitMinAbs,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(-xDipoleExitMinAbs,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(x1min*0.99,y1,zDipoleExit) ## bottom left
  
   window_lin = windowline(window_pts)
   window_lin.SetLineWidth(1)
   window_lin.SetLineColor(ROOT.kBlack)
   return window_pts,window_lin

def getnarrwindow(cluster1,cluster2,i1,i2):
   x1,y1,z1 = getpoint(cluster1,i1)   
   x2,y2,z2 = getpoint(cluster2,i2)
   window_pts = TPolyMarker3D()
   window_pts.SetNextPoint(x1-xAbsMargins,y1,z1) ## top left
   window_pts.SetNextPoint(x1+xAbsMargins,y1,z1) ## top right
   window_pts.SetNextPoint(x2+xAbsMargins,y2,z2) ## bottom right
   window_pts.SetNextPoint(x2-xAbsMargins,y2,z2) ## bottom left
   window_lin = windowline(window_pts)
   window_lin.SetLineWidth(1)
   window_lin.SetLineColor(ROOT.kBlue)
   return window_pts,window_lin

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
   
def isinwindowxz(points,i,window):
   ### get the test point
   rTest = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   points.GetPoint(i,rTest[0],rTest[1],rTest[2])
   if(abs(rTest[0])<1): return False
   ### get the window x extremes at zTest
   xmin,xmax = getxminmax(window,rTest[2])
   if(rTest[0]>=xmin and rTest[0]<=xmax): return True
   return False
   
def isinwindowyz(points,i,window):
   ### get the test point
   rTest = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   points.GetPoint(i,rTest[0],rTest[1],rTest[2])
   if(abs(rTest[0])<1): return False
   ### get the window x extremes at zTest
   ymin,ymax = getyminmax(window,rTest[2])
   if(rTest[1]>=ymin and rTest[1]<=ymax): return True
   return False
   
def trimwide(allpoints,points_wide,windowxz,windowyz,xpivot):
   for layer in layers:
      for i in range(allpoints["Cls"][layer].GetN()):
         acceptxz = isinwindowxz(allpoints["Cls"][layer],i,windowxz)
         acceptyz = isinwindowyz(allpoints["Cls"][layer],i,windowyz)
         if(not acceptxz): continue
         if(not acceptyz): continue
         r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
         allpoints["Cls"][layer].GetPoint(i,r[0],r[1],r[2])
         # if(isel(r[0]) and r[0]>1.001*xpivot):     continue
         if(isel(r[0]) and r[0]>=xpivot):     continue
         # if(not isel(r[0]) and r[0]<1.001*xpivot): continue
         if(not isel(r[0]) and r[0]<=xpivot): continue
         points_wide["Cls"][layer].SetNextPoint(r[0],r[1],r[2])
         points_wide["IsSig"][layer].append(allpoints["IsSig"][layer][i])
         points_wide["TrkId"][layer].append(allpoints["TrkId"][layer][i])

def trimnarr(points_wide,points_narr,window):
   for layer in layers:
      for i in range(points_wide["Cls"][layer].GetN()):
         accept = isinwindowxz(points_wide["Cls"][layer],i,window)
         if(not accept): continue
         r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
         points_wide["Cls"][layer].GetPoint(i,r[0],r[1],r[2])
         points_narr["Cls"][layer].SetNextPoint(r[0],r[1],r[2])
         points_narr["IsSig"][layer].append(points_wide["IsSig"][layer][i])
         points_narr["TrkId"][layer].append(points_wide["TrkId"][layer][i])

def makeseed(cluster1,cluster2,i1,i2,particles):
   p = TLorentzVector()
   r1 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   r2 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   cluster1.GetPoint(i1,r1[0],r1[1],r1[2])
   cluster2.GetPoint(i2,r2[0],r2[1],r2[2])
   
   rnd = TRandom()
   rnd.SetSeed()
   posneg = rnd.Uniform(-1,+1)
   pxgaus = rnd.Gaus(7.2e-4,5.e-4)
   
   x0 = 0
   z0 = zofx(r1,r2,x0)
   xExit = math.abs(xofz(r1,r2,zDipoleExit))*cm2m
   H = math.abs((zDipoleExit-z0))*cm2m
   R = H*(LB)/xExit + xExit ## look this up in my slides
   P = 0.3*B*R
   # P = 0.3*B*R/1.001
   
   v1 = TVector2(r1[2],r1[1])
   v2 = TVector2(r2[2],r2[1])
   u = rUnit2(v2,v1)
   uz = u.X()
   uy = u.Y()
   px = pxgaus if(posneg>=0) else -pxgaus
   py = P*uy
   pz = P*uz
   # p.SetPxPyPzE(px,py,pz,math.sqrt(P*P+me2))
   p.SetPxPyPzE(px,py,pz,math.sqrt(px*px + py*py + pz*pz + me2))
   
   return p

def getparticlename(cluster,i):
   r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   cluster.GetPoint(i,r[0],r[1],r[2])
   if  (r[0]>+1): return "electrons"
   elif(r[0]<-1): return "positrons"
   else:
      print("ERROR: x=",r[0])
      print("  cannot determine particle name")
      quit()
   return "unknown"


#############################################

def isel(x):
   if(x<0): return False
   return True


### define the data structures to hold clusters and windows
def initpoints():
   points = { 
              "IsSig" : {1:[],2:[],3:[],4:[]},
              "TrkId" : {1:[],2:[],3:[],4:[]},
              "Cls"   : {1:TPolyMarker3D(),2:TPolyMarker3D(),3:TPolyMarker3D(),4:TPolyMarker3D()}
   }   
   return points
   
def AddPoint(points,rcls,sig=False,trkid=-1):
   layer = -1
   if  (rcls[2]==330): layer = 4
   elif(rcls[2]==320): layer = 3
   elif(rcls[2]==310): layer = 2
   elif(rcls[2]==300): layer = 1
   else:
      print("ERROR: cannot add cluster at layer with z=",rcls[2])
      quit()
   points["Cls"][layer].SetNextPoint(rcls[0],rcls[1],rcls[2])
   # if(sig): points["Cls"][layer].SetMarkerColor(ROOT.kRed)
   points["IsSig"][layer].append(sig)
   points["TrkId"][layer].append(trkid)


def getNnon0(clusters):
   N = 0
   for j1 in range(clusters.GetN()):
      r = [ ROOT.Double(), ROOT.Double(), ROOT.Double() ]
      clusters.GetPoint(j1,r[0],r[1],r[2]) ### the clusters
      if(r[2]==0): break
      N += 1
   return N
   
def getLogicSidesArr():
   sidesarr = ["Eside","Pside"]
   if(proc=="trident"): sidesarr = ["Pside"]
   if(proc=="bppp" and sides=="e+"): sidesarr = ["Pside"]
   if(proc=="bppp" and sides=="e-"): sidesarr = ["Eside"]
   return sidesarr
   
def drawall(name,pointsEside,pointsPside,dodraw):
   if(not dodraw): return
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.ShowAxis()
   # view.SetRange(-80,-50,0, +80,+50,350)
   if  (sides=="e-"): view.SetRange(0,-10,190, xEsideR,+10,340)
   elif(sides=="e+"): view.SetRange(xPsideL,-10,190, 0,+10,340)
   else:              view.SetRange(xPsideL,-10,190, xEsideR,+10,340)
   geom = getgeometry() ### get the geometry from the root file
   for g in geom: g.Draw()
   for layer in layers:
      pointsEside["Cls"][layer].Draw("same")
      pointsPside["Cls"][layer].Draw("same")
   cnv.SaveAs(name)
   
   
def draw(name,points,dodraw,particles="",window_yz=None,window_xz=None):
   if(not dodraw): return
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.ShowAxis()
   # view.SetRange(-80,-50,0, +80,+50,350)
   if  (particles=="electrons"): view.SetRange(0,-10,190, xEsideR,+10,340)
   elif(particles=="positrons"): view.SetRange(xPsideL,-10,190, 0,+10,340)
   else:                         view.SetRange(xPsideL,-10,190, xEsideR,+10,340)
   geom = getgeometry() ### get the geometry from the root file
   for g in geom: g.Draw()
   if(window_yz is not None): window_yz.Draw("same")
   if(window_xz is not None): window_xz.Draw("same")
   for layer in layers: points["Cls"][layer].Draw("same")
   cnv.SaveAs(name)
   

def line3d(z,m_xz,c_xz,m_yz,c_yz):
   # the intersection of those two planes and
   # the function for the line would be:
   # z = m_yz * y + c_yz
   # z = m_xz * x + c_xz
   # or:
   x = (z - c_xz)/m_xz
   y = (z - c_yz)/m_yz
   return x,y

def seed3dfit(name,r1,r2,r3,r4,dodraw):
   g = TGraph2D()
   g.SetMarkerSize(3)
   g.SetMarkerStyle(20)
   g.SetMarkerColor(ROOT.kRed)
   g.SetLineColor(ROOT.kRed)
   g.SetPoint(0,r1[2],r1[1],r1[0])
   g.SetPoint(1,r2[2],r2[1],r2[0])
   g.SetPoint(2,r3[2],r3[1],r3[0])
   g.SetPoint(3,r4[2],r4[1],r4[0])
   g.GetXaxis().SetRangeUser(290,340)
   g.GetYaxis().SetRangeUser(-0.8,+0.8)
   g.GetZaxis().SetRangeUser(0 if(isel(r1[0])) else xPsideL, xEsideR if(isel(r1[0])) else 0 )
   
   x = np.array([r1[0],r2[0],r3[0],r4[0]])
   y = np.array([r1[1],r2[1],r3[1],r4[1]])
   z = np.array([r1[2],r2[2],r3[2],r4[2]])
   
   # this will find the slope and x-intercept of a plane
   # parallel to the y-axis that best fits the data
   A_xz = np.vstack((x, np.ones(len(x)))).T
   # m_xz, c_xz = np.linalg.lstsq(A_xz, z,rcond=None)[0]
   result_xz = np.linalg.lstsq(A_xz, z,rcond=None)
   m_xz, c_xz = result_xz[0]
   residuals_xz = result_xz[1]

   # again for a plane parallel to the x-axis
   A_yz = np.vstack((y, np.ones(len(y)))).T
   # m_yz, c_yz = np.linalg.lstsq(A_yz, z,rcond=None)[0]
   result_yz = np.linalg.lstsq(A_yz, z,rcond=None)
   m_yz, c_yz = result_yz[0]
   residuals_yz = result_yz[1]

   if(dodraw):
      zz = np.array([300,310,320,330])
      xx,yy = line3d(zz, m_xz,c_xz,m_yz,c_yz)
      lfit = TPolyLine3D()
      for i in range(4):
         lfit.SetNextPoint(zz[i],yy[i],xx[i])
      lfit.SetLineColor(ROOT.kBlue)
      cnv = TCanvas("","",2000,2000)
      view = TView.CreateView(1)
      xviewmin = 0 if(isel(r1[0])) else xPsideL
      xviewmax = xEsideR if(isel(r1[0])) else 0
      view.SetRange(290,-0.8, xviewmin , 340,+0.8,xviewmax)
      view.ShowAxis()
      g.Draw("p0")
      lfit.Draw("smae")
      cnv.SaveAs(name)
   
   return residuals_xz,residuals_yz

def seed3dfitSVD(name,r1,r2,r3,r4,dodraw):
   g = TGraph2D()
   g.SetMarkerSize(3)
   g.SetMarkerStyle(20)
   g.SetMarkerColor(ROOT.kRed)
   g.SetLineColor(ROOT.kRed)
   g.SetPoint(0,r1[2],r1[1],r1[0])
   g.SetPoint(1,r2[2],r2[1],r2[0])
   g.SetPoint(2,r3[2],r3[1],r3[0])
   g.SetPoint(3,r4[2],r4[1],r4[0])
   g.GetXaxis().SetRangeUser(290,340)
   g.GetYaxis().SetRangeUser(-0.8,+0.8)
   g.GetZaxis().SetRangeUser( 0 if(isel(r1[0])) else xPsideL, xEsideR if(isel(r1[0])) else 0 )
   
   x = np.array([r1[0],r2[0],r3[0],r4[0]])
   y = np.array([r1[1],r2[1],r3[1],r4[1]])
   z = np.array([r1[2],r2[2],r3[2],r4[2]])

   data = np.concatenate((x[:, np.newaxis], 
                          y[:, np.newaxis], 
                          z[:, np.newaxis]), 
                         axis=1)

   # Calculate the mean of the points, i.e. the 'center' of the cloud
   datamean = data.mean(axis=0)

   # Do an SVD on the mean-centered data (Singular Value Decomposition)
   uu, dd, vv = np.linalg.svd(data - datamean) 

   # Now vv[0] contains the first principal component, i.e. the direction
   # vector of the 'best fit' line in the least squares sense.

   # Now generate some points along this best fit line, for plotting.

   # I use -7, 7 since the spread of the data is roughly 14
   # and we want it to have mean 0 (like the points we did
   # the svd on). Also, it's a straight line, so we only need 2 points.
   # linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
   linepts = vv[0] * np.mgrid[-50:50:2j][:, np.newaxis]

   # shift by the mean to get the line in the right place
   linepts += datamean
   
   if(dodraw):
      lfit = TPolyLine3D()
      for point in linepts:
         lfit.SetNextPoint(point[2],point[1],point[0])
      lfit.SetLineColor(ROOT.kBlue)
      cnv = TCanvas("","",2000,2000)
      view = TView.CreateView(1)
      xviewmin = 0 if(isel(r1[0])) else xPsideL
      xviewmax = xEsideR if(isel(r1[0])) else 0
      view.SetRange(290,-0.8, xviewmin , 340,+0.8,xviewmax)
      view.ShowAxis()
      g.Draw("p0")
      lfit.Draw("smae")
      cnv.SaveAs(name)
   
   return linepts, dd ## dd is a 1D array of the data singular values
   
   
def seed2dfit(name,r1,r2,r3,r4,dodraw):
   gxyz = TGraph2D()
   gxyz.SetMarkerSize(1)
   gxyz.SetMarkerStyle(24)
   gxyz.SetMarkerColor(ROOT.kBlack)
   gxyz.SetPoint(0,r1[2],r1[1],r1[0])
   gxyz.SetPoint(1,r2[2],r2[1],r2[0])
   gxyz.SetPoint(2,r3[2],r3[1],r3[0])
   gxyz.SetPoint(3,r4[2],r4[1],r4[0])
   gxyz.GetXaxis().SetRangeUser(290,340)
   gxyz.GetYaxis().SetRangeUser(-0.8,+0.8)
   gxyz.GetZaxis().SetRangeUser( 0 if(isel(r1[0])) else xPsideL, xEsideR if(isel(r1[0])) else 0 )
      
   gxz = TGraph()
   gxz.SetMarkerSize(1)
   gxz.SetMarkerStyle(24)
   gxz.SetMarkerColor(ROOT.kBlack)
   gxz.SetPoint(0,r1[2],r1[0])
   gxz.SetPoint(1,r2[2],r2[0])
   gxz.SetPoint(2,r3[2],r3[0])
   gxz.SetPoint(3,r4[2],r4[0])
   gxz.GetXaxis().SetRangeUser(290,340)
   if(isel(r1[0])): gxz.GetYaxis().SetRangeUser(0,xEsideR)
   else:            gxz.GetYaxis().SetRangeUser(xPsideL,0)
   
   gyz = TGraph()
   gyz.SetMarkerSize(1)
   gyz.SetMarkerStyle(24)
   gyz.SetMarkerColor(ROOT.kBlack)
   gyz.SetLineColor(ROOT.kRed)
   gyz.SetPoint(0,r1[2],r1[1])
   gyz.SetPoint(1,r2[2],r2[1])
   gyz.SetPoint(2,r3[2],r3[1])
   gyz.SetPoint(3,r4[2],r4[1])
   gyz.GetXaxis().SetRangeUser(290,340)
   gyz.GetYaxis().SetRangeUser(-0.8,+0.8)
   
   fxz = TF1("fxz","pol1",300,330)
   fitr_xz = gxz.Fit(fxz,"Q")
   fitf_xz = gxz.FindObject("fxz")
   
   fyz = TF1("fyz","pol1",300,330)
   fitr_yz = gyz.Fit(fyz,"Q")
   fitf_yz = gyz.FindObject("fyz")
   
   if(dodraw):
      lfit = TPolyLine3D()
      for z in [300,310,320,330]:
         x = fitf_xz.Eval(z)
         y = fitf_yz.Eval(z)
         lfit.SetNextPoint(z,y,x)
      lfit.SetLineColor(ROOT.kBlue)
      cnv = TCanvas("","",1500,500)
      cnv.Divide(3,1)
      cnv.cd(1)
      view = TView.CreateView(1)
      xviewmin = 0 if(isel(r1[0])) else xPsideL
      xviewmax = xEsideR if(isel(r1[0])) else 0
      view.SetRange(290,-0.8, xviewmin , 340,+0.8,xviewmax)
      view.ShowAxis()
      gxyz.Draw("p0")
      lfit.Draw("smae")
      cnv.cd(2)
      gxz.Draw()
      cnv.cd(3)
      gyz.Draw()
      cnv.SaveAs(name)
   
   chi2_xz = fitf_xz.GetChisquare()/fitf_xz.GetNDF()
   chi2_yz = fitf_yz.GetChisquare()/fitf_yz.GetNDF()
   prob_xz = fitf_xz.GetProb()
   prob_yz = fitf_yz.GetProb()
   return  chi2_xz,prob_xz,chi2_yz,prob_yz


##########################################################################################
##########################################################################################
##########################################################################################
### read fits to root file
# fitsEx = GetFits()

svd0Seed = ROOT.std.vector( float )()
svd1Seed = ROOT.std.vector( float )()
svd2Seed = ROOT.std.vector( float )()
chi2xzSeed = ROOT.std.vector( float )()
chi2yzSeed = ROOT.std.vector( float )()
residxzSeed = ROOT.std.vector( float )()
residyzSeed = ROOT.std.vector( float )()
issigSeed = ROOT.std.vector( int )()
iGenMatch = ROOT.std.vector( int )()
x1Seed = ROOT.std.vector( float )()
y1Seed = ROOT.std.vector( float )()
z1Seed = ROOT.std.vector( float )()
x2Seed = ROOT.std.vector( float )()
y2Seed = ROOT.std.vector( float )()
z2Seed = ROOT.std.vector( float )()
x3Seed = ROOT.std.vector( float )()
y3Seed = ROOT.std.vector( float )()
z3Seed = ROOT.std.vector( float )()
x4Seed = ROOT.std.vector( float )()
y4Seed = ROOT.std.vector( float )()
z4Seed = ROOT.std.vector( float )()

# x1Cluster_intrksys = ROOT.std.vector( float )()
# y1Cluster_intrksys = ROOT.std.vector( float )()
# z1Cluster_intrksys = ROOT.std.vector( float )()
# x2Cluster_intrksys = ROOT.std.vector( float )()
# y2Cluster_intrksys = ROOT.std.vector( float )()
# z2Cluster_intrksys = ROOT.std.vector( float )()
# x3Cluster_intrksys = ROOT.std.vector( float )()
# y3Cluster_intrksys = ROOT.std.vector( float )()
# z3Cluster_intrksys = ROOT.std.vector( float )()
# x4Cluster_intrksys = ROOT.std.vector( float )()
# y4Cluster_intrksys = ROOT.std.vector( float )()
# z4Cluster_intrksys = ROOT.std.vector( float )()

pxSeed = ROOT.std.vector( float )()
pySeed = ROOT.std.vector( float )()
pzSeed = ROOT.std.vector( float )()
eSeed  = ROOT.std.vector( float )()
pxGen  = ROOT.std.vector( float )()
pyGen  = ROOT.std.vector( float )()
pzGen  = ROOT.std.vector( float )()
eGen   = ROOT.std.vector( float )()
qGen   = ROOT.std.vector( float )()
iGen   = ROOT.std.vector( int )()
tF = TFile("../data/root/seeds_"+proc+".root","RECREATE")
tF.cd()
tT = TTree("seeds","seeds")

tT.Branch('svd0Seed',svd0Seed)
tT.Branch('svd1Seed',svd1Seed)
tT.Branch('svd2Seed',svd2Seed)
tT.Branch('chi2xzSeed',chi2xzSeed)
tT.Branch('chi2yzSeed',chi2yzSeed)
tT.Branch('residxzSeed',residxzSeed)
tT.Branch('residyzSeed',residyzSeed)
tT.Branch('issigSeed',issigSeed)
tT.Branch('iGenMatch',iGenMatch)
tT.Branch('x1Seed',x1Seed)
tT.Branch('y1Seed',y1Seed)
tT.Branch('z1Seed',z1Seed)
tT.Branch('x2Seed',x2Seed)
tT.Branch('y2Seed',y2Seed)
tT.Branch('z2Seed',z2Seed)
tT.Branch('x3Seed',x3Seed)
tT.Branch('y3Seed',y3Seed)
tT.Branch('z3Seed',z3Seed)
tT.Branch('x4Seed',x4Seed)
tT.Branch('y4Seed',y4Seed)
tT.Branch('z4Seed',z4Seed)
tT.Branch('pxSeed',pxSeed)
tT.Branch('pySeed',pySeed)
tT.Branch('pzSeed',pzSeed)
tT.Branch('eSeed',eSeed)
tT.Branch('pxGen',pxGen)
tT.Branch('pyGen',pyGen)
tT.Branch('pzGen',pzGen)
tT.Branch('eGen',eGen)
tT.Branch('qGen',qGen)
tT.Branch('iGen',iGen)


histos = { "h_residuals_xz_sig": TH1D("residuals_xz_sig",";residuals_{xz};Tracks", 500,0,0.5),
           "h_residuals_yz_sig": TH1D("residuals_yz_sig",";residuals_{yz};Tracks", 500,0,500), 
           "h_residuals_xz_bkg": TH1D("residuals_xz_bkg",";residuals_{xz};Tracks", 500,0,0.5),
           "h_residuals_yz_bkg": TH1D("residuals_yz_bkg",";residuals_{yz};Tracks", 500,0,500),
           
           "h_svd_dd0_sig": TH1D("svd_dd0_sig",";svd_{dd0};Tracks", 500,21,24),
           "h_svd_dd0_bkg": TH1D("svd_dd0_bkg",";svd_{dd0};Tracks", 500,21,24), 
           "h_svd_dd1_sig": TH1D("svd_dd1_sig",";svd_{dd1};Tracks", 500,0,0.1),
           "h_svd_dd1_bkg": TH1D("svd_dd1_bkg",";svd_{dd1};Tracks", 500,0,0.1),
           "h_svd_dd2_sig": TH1D("svd_dd2_sig",";svd_{dd2};Tracks", 500,0,0.05),
           "h_svd_dd2_bkg": TH1D("svd_dd2_bkg",";svd_{dd2};Tracks", 500,0,0.05),
           
           "h_prob_xz_sig": TH1D("prob_xz_sig",";prob_{xz};Tracks", 500,0,1.0),
           "h_prob_yz_sig": TH1D("prob_yz_sig",";prob_{yz};Tracks", 500,0,1.0), 
           "h_prob_xz_bkg": TH1D("prob_xz_bkg",";prob_{xz};Tracks", 500,0,1.0),
           "h_prob_yz_bkg": TH1D("prob_yz_bkg",";prob_{yz};Tracks", 500,0,1.0),
           
           "h_chi2ndf_xz_sig": TH1D("chi2ndf_xz_sig",";chi2ndf_{xz};Tracks", 500,0,0.001),
           "h_chi2ndf_yz_sig": TH1D("chi2ndf_yz_sig",";chi2ndf_{yz};Tracks", 500,0,0.001), 
           "h_chi2ndf_xz_bkg": TH1D("chi2ndf_xz_bkg",";chi2ndf_{xz};Tracks", 500,0,0.001),
           "h_chi2ndf_yz_bkg": TH1D("chi2ndf_yz_bkg",";chi2ndf_{yz};Tracks", 500,0,0.001),
           
           "h_seed_resE" : TH1D("seed_resE", ";(E_{seed}-E_{gen})/E_{gen};Tracks",    100,-0.05,+0.05),
           "h_seed_resPz": TH1D("seed_resPz",";(Pz_{seed}-Pz_{gen})/Pz_{gen};Tracks", 100,-0.05,+0.05), 
           "h_seed_resPy": TH1D("seed_resPy",";(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,-10,+10),
           
           "h_seed_resE_vs_x"  : TH2D("seed_resE_vs_x",  ";x;(E_{seed}-E_{gen})/E_{gen};Tracks",    100,detXmin,detXmax, 100,-0.05,+0.05),
           "h_seed_resPy_vs_x" : TH2D("seed_resPy_vs_x", ";x;(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,detXmin,detXmax, 100,-10,+10),
           
           "h_N_sigacc":        TH1D("N_sigacc",        ";Track multiplicity;Events", 100,30,330),
           "h_N_all_seeds":     TH1D("N_all_seeds",     ";Track multiplicity;Events", 100,30,330),
           "h_N_matched_seeds": TH1D("N_matched_seeds", ";Track multiplicity;Events", 100,30,330),
           "h_N_good_seeds":    TH1D("N_good_seeds",    ";Track multiplicity;Events", 100,30,330),
           
           "h_seeding_score": TH1D("h_seeding_score", ";N_{seeds}^{matched}/N_{signa}^{in.acc} [%];Events", 20,91,101),
           "h_seeding_pool":  TH1D("h_seeding_pool",  ";N_{seeds}^{all}/N_{signa}^{in.acc} [%];Events", 50,90,590),
}
sidesarr = getLogicSidesArr()
pdfname = "../output/pdf/seedingdemo_"+proc+".pdf"
intfile = TFile("../data/root/rec_"+proc+".root","READ")
intree = intfile.Get("res")
nevents = intree.GetEntries()
print("with %d events" % nevents)
nmax = 1000000
n=0 ### init n
for event in intree:
   Nsigall = 0
   Nsigacc = 0
   Nseeds = 0
   Nmatched = 0
   Ngood = 0
   
   ## clear the output vectors
   svd0Seed.clear()
   svd1Seed.clear()
   svd2Seed.clear()
   chi2xzSeed.clear()
   chi2yzSeed.clear()
   residxzSeed.clear()
   residyzSeed.clear()
   issigSeed.clear()
   iGenMatch.clear()
   x1Seed.clear()
   y1Seed.clear()
   z1Seed.clear()
   x2Seed.clear()
   y2Seed.clear()
   z2Seed.clear()
   x3Seed.clear()
   y3Seed.clear()
   z3Seed.clear()
   x4Seed.clear()
   y4Seed.clear()
   z4Seed.clear()
   
   # x1Cluster_intrksys.clear()
   # y1Cluster_intrksys.clear()
   # z1Cluster_intrksys.clear()
   # x2Cluster_intrksys.clear()
   # y2Cluster_intrksys.clear()
   # z2Cluster_intrksys.clear()
   # x3Cluster_intrksys.clear()
   # y3Cluster_intrksys.clear()
   # z3Cluster_intrksys.clear()
   # x4Cluster_intrksys.clear()
   # y4Cluster_intrksys.clear()
   # z4Cluster_intrksys.clear()
   
   pxSeed.clear()
   pySeed.clear()
   pzSeed.clear()
   eSeed.clear()
   
   pxGen.clear()
   pyGen.clear()
   pzGen.clear()
   eGen.clear()
   qGen.clear()
   iGen.clear()
   
   
   
   ### start the loop
   if(n>nmax): break
   
   ### draw?
   dodraw = (n<=NeventsToDraw)
   
   ### container for all clusters
   allpointsEside = initpoints()
   allpointsPside = initpoints()
   
   ## clusters' vectors are always written out (even if empty) for all gen tracks!
   ## each generated track in the vector always has 4 clusters accessed via TPolyMarker3D::GetPoint()
   for i in range(event.polm_clusters.size()):
      
      ###############################################################
      if(proc=="bppp"):
         if(sides=="e+" and event.qgen[i]<0): continue ## only positrons
         if(sides=="e-" and event.qgen[i]>0): continue ## only electrons
      if(proc=="trident" and event.qgen[i]<0): continue ## only positrons
      ###############################################################
      
      Nsigall += 1
      
      wgt  = event.wgtgen[i]
      pgen = event.pgen[i]
      ### cut on acceptance
      if(event.acctrkgen[i]!=1): continue
      
      Nsigacc += 1
      
      ### write the truth track momentum and its index
      pxGen.push_back(pgen.Px())
      pyGen.push_back(pgen.Py())
      pzGen.push_back(pgen.Pz())
      eGen.push_back(pgen.E())
      qGen.push_back(event.qgen[i])
      iGen.push_back(i)
      
      ### loop over all clusters of the track and put in the allpoints classified by the layer
      for jxy in range(event.polm_clusters[i].GetN()):
         rcls = [ ROOT.Double(), ROOT.Double(), ROOT.Double() ]
         event.polm_clusters[i].GetPoint(jxy,rcls[0],rcls[1],rcls[2]) ### the clusters
         if(rcls[0]>0): AddPoint(allpointsEside,rcls,True,i)
         if(rcls[0]<0): AddPoint(allpointsPside,rcls,True,i)
   Nsig4 = getNnon0(allpointsEside["Cls"][4])+getNnon0(allpointsPside["Cls"][4])
   
   
   ### embed some noise ***clusters***
   rnd = TRandom()
   rnd.SetSeed()
   for kN in range(NnoiseClusters):
      for layer in layers:
         for side in sidesarr:
            x = 0
            if(side=="Pside"): x = rnd.Uniform(xPsideL,xPsideR)
            if(side=="Eside"): x = rnd.Uniform(xEsideL,xEsideR)
            y = rnd.Uniform(-0.75,+0.75)
            if(layer==1): z = 300
            if(layer==2): z = 310
            if(layer==3): z = 320
            if(layer==4): z = 330
            rnoise = [x,y,z]
            if(side=="Pside"): AddPoint(allpointsPside,rnoise)
            if(side=="Eside"): AddPoint(allpointsEside,rnoise)
   
   ### embed background ***tracks***
   resolution = 0.001 ## cm (10 um)
   rnd = TRandom()
   rnd.SetSeed()
   nbkgtracks = 0
   while (nbkgtracks<NbkgTracks):
      # production vertex is uniform on the inner perimeter of the beampipe
      R = cfgmap["Rbeampipe"]-cfgmap["Wbeampipe"]
      phi = rnd.Uniform(0,2*ROOT.TMath.Pi())
      x0 = R*ROOT.TMath.Cos(phi)
      y0 = R*ROOT.TMath.Sin(phi)
      z0 = cfgmap["zDipoleExit"]
      # chose a point in the exit of the tracker 
      z4 = cfgmap["zLayer4"]
      x4 = rnd.Uniform(1.1*cfgmap["xPsideL"],   1.1*cfgmap["xEsideR"])
      y4 = rnd.Uniform(-1.1*cfgmap["Hstave"]/2, 1.1*cfgmap["Hstave"]/2)
      # require that the track is not pointing in the reange between the 2 trackr sides
      if(x4>cfgmap["xPsideR"] and x4<cfgmap["xEsideL"]): continue
      nbkgtracks += 1 
      
      # place 3 points on layers 1,2,3 along the line between r4 and r0, with some error (resolution)
      r0 = [x0,y0,z0]
      r4 = [x4,y4,z4]
      z1 = cfgmap["zLayer1"]
      x1 = xofz(r4,r0,z1)+rnd.Gaus(0,resolution)
      y1 = yofz(r4,r0,z1)+rnd.Gaus(0,resolution)
      z2 = cfgmap["zLayer2"]
      x2 = xofz(r4,r0,z2)+rnd.Gaus(0,resolution)
      y2 = yofz(r4,r0,z2)+rnd.Gaus(0,resolution)
      z3 = cfgmap["zLayer3"]
      x3 = xofz(r4,r0,z3)+rnd.Gaus(0,resolution)
      y3 = yofz(r4,r0,z3)+rnd.Gaus(0,resolution)
      r1 = [x1,y1,z1]
      r2 = [x2,y2,z2]
      r3 = [x3,y3,z3]
      if(side=="Pside"):
         AddPoint(allpointsPside,r1)
         AddPoint(allpointsPside,r2)
         AddPoint(allpointsPside,r3)
      if(side=="Eside"):
         AddPoint(allpointsEside,r1)
         AddPoint(allpointsEside,r2)
         AddPoint(allpointsEside,r3)

   ### all points (signal, background and noise)   
   Nbkgsig4 = getNnon0(allpointsEside["Cls"][4])+getNnon0(allpointsPside["Cls"][4])
   
   
   ### just draw the full event
   drawall(pdfname+"(",allpointsEside,allpointsPside,dodraw)
   
   
   ### loop on the 2 sides
   for side in sidesarr:
      allpoints = allpointsEside if(side=="Eside") else allpointsPside
      
      ### the initial pool for pivot clusters
      Nall4 = getNnon0(allpoints["Cls"][4])
      
      ### loop over the clusters and start the seeding
      for j4 in range(Nall4):
         r4 = getpoint(allpoints["Cls"][4],j4)
         xpivot = r4[0]
         ### electron / positron?
         particlename = getparticlename(allpoints["Cls"][4],j4)
         ### get the yz window
         winpts_yz,winlin_yz = getyzwindow(allpoints["Cls"][4],j4)
         ### set the wide window starting from cluster_seed1 (window corners must be added clockwise!)
         winpts_xz_wide,winlin_xz_wide = getwidewindow(allpoints["Cls"][4],j4)
         ### discard all clusters which are not in the wide window
         widepoints = initpoints()
         trimwide(allpoints,widepoints,winpts_xz_wide,winpts_yz,xpivot)
         Nwide1 = getNnon0(widepoints["Cls"][1])
         # if(Nwide1<1): print("Failed Nwide1")
         ### draw the wide window
         draw(pdfname,widepoints,dodraw,particlename,winlin_yz,winlin_xz_wide)
         
         ### choose one cluster in layer 1 as the second seed
         for j1 in range(Nwide1):
            ### get the narrow window  (window corners must be added clockwise!)
            winpts_xz_narr,winlin_xz_narr = getnarrwindow(allpoints["Cls"][4],widepoints["Cls"][1],j4,j1)
            ### discard all clusters which are not in the narrow window
            narrpoints = initpoints()
            trimnarr(widepoints,narrpoints,winpts_xz_narr)
            Nnarr2 = getNnon0(narrpoints["Cls"][2])
            Nnarr3 = getNnon0(narrpoints["Cls"][3])
            
            ### check if there are at least 1 cluster in both layer 2 and layer 3 within the narrow window
            if(Nnarr2<1 or Nnarr3<1): continue
            
            ### draw the narrow window
            draw(pdfname,narrpoints,dodraw,particlename,None,winlin_xz_narr)
            
            ### get the seed - note: could be that there are several combinations but the seed momentum would be identical
            pseed = makeseed(allpoints["Cls"][4],widepoints["Cls"][1],j4,j1,particlename)
            # if(pseed.E()>Emax or pseed.E()<Emin): print("pseed.E()>Emax or pseed.E()<Emin")
            if(pseed.E()>Emax or pseed.E()<Emin): continue
            
            ### set the cluster in layer 1
            r1 = getpoint(widepoints["Cls"][1],j1)
            
            ### loop on the clusters in layer 2 and 3:
            for j2 in range(Nnarr2):
            # for j2 in range(narrpoints["Cls"][2].GetN()):
               r2 = getpoint(narrpoints["Cls"][2],j2)
               for j3 in range(Nnarr3):
               # for j3 in range(narrpoints["Cls"][3].GetN()):
                  r3 = getpoint(narrpoints["Cls"][3],j3)
                  Nseeds += 1
                  
                  issig1 = widepoints["IsSig"][1][j1]
                  issig2 = narrpoints["IsSig"][2][j2]
                  issig3 = narrpoints["IsSig"][3][j3]
                  issig4 = allpoints["IsSig"][4][j4]
                  
                  trkid4 = allpoints["TrkId"][4][j4]
                  trkid3 = narrpoints["TrkId"][3][j3]
                  trkid2 = narrpoints["TrkId"][2][j2]
                  trkid1 = widepoints["TrkId"][1][j1]
                  issig = (issig4 and issig1 and issig2 and issig3)
                  trkid = str(trkid4) if(trkid4==trkid1 and trkid1==trkid2 and trkid2==trkid3) else "mult"
                  issiguniq = (issig and trkid!="mult")
                  if(issiguniq): Nmatched += 1
      
                  ### two independent 2d fits
                  chi2_xz,prob_xz,chi2_yz,prob_yz = seed2dfit(pdfname,r1,r2,r3,r4,dodraw)
                  if(issiguniq):
                     histos["h_chi2ndf_xz_sig"].Fill(chi2_xz)
                     histos["h_chi2ndf_yz_sig"].Fill(chi2_yz)
                     histos["h_prob_xz_sig"].Fill(prob_xz)
                     histos["h_prob_yz_sig"].Fill(prob_yz)
                  else:
                     histos["h_chi2ndf_xz_bkg"].Fill(chi2_xz)
                     histos["h_chi2ndf_yz_bkg"].Fill(chi2_yz)
                     histos["h_prob_xz_bkg"].Fill(prob_xz)
                     histos["h_prob_yz_bkg"].Fill(prob_yz)
      
                  ### a single 3d fit
                  res_xz, res_yz = seed3dfit(pdfname,r1,r2,r3,r4,dodraw)
                  if(issiguniq):
                     histos["h_residuals_xz_sig"].Fill(res_xz)
                     histos["h_residuals_yz_sig"].Fill(res_yz)
                  else:
                     histos["h_residuals_xz_bkg"].Fill(res_xz)
                     histos["h_residuals_yz_bkg"].Fill(res_yz)
                     
                  ### a single 3d fit SVD
                  lfitpts, dd = seed3dfitSVD(pdfname,r1,r2,r3,r4,dodraw)
                  
                  ### set again the pseed according to lfit
                  cluster1 = TPolyMarker3D()
                  cluster2 = TPolyMarker3D()
                  cluster1.SetNextPoint(lfitpts[0][0],lfitpts[0][1],lfitpts[0][2])
                  cluster2.SetNextPoint(lfitpts[1][0],lfitpts[1][1],lfitpts[1][2])
                  if(isel(lfitpts[0][0])): pseed = makeseed(cluster2,cluster1,0,0,particlename)
                  else:                    pseed = makeseed(cluster1,cluster2,0,0,particlename)
                  if(pseed.E()>Emax or pseed.E()<Emin): continue
                  
                  ### the SVD alg
                  if(issiguniq):
                     histos["h_svd_dd0_sig"].Fill(dd[0])
                     histos["h_svd_dd1_sig"].Fill(dd[1])
                     histos["h_svd_dd2_sig"].Fill(dd[2])
                  else:
                     histos["h_svd_dd0_bkg"].Fill(dd[0])
                     histos["h_svd_dd1_bkg"].Fill(dd[1])
                     histos["h_svd_dd2_bkg"].Fill(dd[2])
                  
                  ### get the generated matched track momentum
                  pgen = TLorentzVector()
                  igen = -1
                  for k in range(iGen.size()):
                     if(iGen[k]==trkid4):
                        igen = k
                        break
                  
                  ### write out the good seeds
                  svd0Seed.push_back(dd[0])
                  svd1Seed.push_back(dd[1])
                  svd2Seed.push_back(dd[2])
                  chi2xzSeed.push_back(chi2_xz)
                  chi2yzSeed.push_back(chi2_yz)
                  residxzSeed.push_back(res_xz)
                  residyzSeed.push_back(res_yz)
                  issigSeed.push_back(issiguniq)
                  iGenMatch.push_back(igen)
                  x1Seed.push_back(r1[0])
                  y1Seed.push_back(r1[1])
                  z1Seed.push_back(r1[2])
                  x2Seed.push_back(r2[0])
                  y2Seed.push_back(r2[1])
                  z2Seed.push_back(r2[2])
                  x3Seed.push_back(r3[0])
                  y3Seed.push_back(r3[1])
                  z3Seed.push_back(r3[2])
                  x4Seed.push_back(r4[0])
                  y4Seed.push_back(r4[1])
                  z4Seed.push_back(r4[2])
                  
                  pxSeed.push_back(pseed.Px())
                  pySeed.push_back(pseed.Py())
                  pzSeed.push_back(pseed.Pz())
                  eSeed.push_back(pseed.E())
                  
                  ### cut on some quality
                  isgood = (dd[1]<0.005 and dd[2]<0.0025)
                  if(not isgood): continue
                  Ngood += 1
                  
                  ### check perforrmance of seeding
                  pgen.SetPxPyPzE(pxGen[igen],pyGen[igen],pzGen[igen],eGen[igen])
                  resE = (pseed.E()-pgen.E())/pgen.E()
                  resPz = (pseed.Pz()-pgen.Pz())/pgen.Pz()
                  resPy = (pseed.Py()-pgen.Py())/pgen.Py()
                  histos["h_seed_resE"].Fill(resE)
                  histos["h_seed_resPz"].Fill(resPz)
                  histos["h_seed_resPy"].Fill(resPy)
                  histos["h_seed_resE_vs_x"].Fill(r4[0],resE)
                  histos["h_seed_resPy_vs_x"].Fill(r4[0],resPy)
                  
   histos["h_N_sigacc"].Fill(Nsigacc)       
   histos["h_N_all_seeds"].Fill(Nseeds)       
   histos["h_N_matched_seeds"].Fill(Nmatched)       
   histos["h_N_good_seeds"].Fill(Ngood)
   histos["h_seeding_score"].Fill(Nmatched/Nsigacc*100)
   histos["h_seeding_pool"].Fill(Nseeds/Nsigacc*100)
      
   if(dodraw):
      cnv = TCanvas("","",2000,2000)
      cnv.SaveAs(pdfname+")")
   print("Event: %g --> Nsigall=%g, Nsigacc=%g, Nseeds=%g, Nmatched=%g, Ngood=%g --> Seeds matching performance: Nmatched/Nsigacc=%5.1f%%" % (n,Nsigall,Nsigacc,Nseeds,Nmatched,Ngood,Nmatched/Nsigacc*100))

   tT.Fill()
   if(n%10==0 and n>0): print("  processed %d events" % n)
   n+=1
print("Total events processed: ",n)










cnv = TCanvas("","",1000,1000)
cnv.Divide(2,2)
cnv.cd(1)
histos["h_chi2ndf_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_xz_sig"].Draw()
histos["h_chi2ndf_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_xz_bkg"].Draw("same")
cnv.cd(2)
histos["h_chi2ndf_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_yz_sig"].Draw()
histos["h_chi2ndf_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_yz_bkg"].Draw("same")
cnv.cd(3)
histos["h_prob_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_prob_xz_sig"].Draw()
histos["h_prob_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_prob_xz_bkg"].Draw("same")
cnv.cd(4)
histos["h_prob_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_prob_yz_sig"].Draw()
histos["h_prob_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_prob_yz_bkg"].Draw("same")
cnv.SaveAs("../output/pdf/chi2ndf_"+proc+".pdf")

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_residuals_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_xz_sig"].Draw()
histos["h_residuals_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_xz_bkg"].Draw("same")
cnv.cd(2)
histos["h_residuals_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_yz_sig"].Draw()
histos["h_residuals_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_yz_bkg"].Draw("same")
cnv.SaveAs("../output/pdf/resid3dfit_"+proc+".pdf")

cnv = TCanvas("","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
histos["h_svd_dd0_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd0_sig"].Draw()
histos["h_svd_dd0_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd0_bkg"].Draw("same")
cnv.cd(2)
histos["h_svd_dd1_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd1_sig"].Draw()
histos["h_svd_dd1_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd1_bkg"].Draw("same")
cnv.cd(3)
histos["h_svd_dd2_sig"].SetLineColor(ROOT.kRed);   histos["h_svd_dd2_sig"].Draw()
histos["h_svd_dd2_bkg"].SetLineColor(ROOT.kBlack); histos["h_svd_dd2_bkg"].Draw("same")
cnv.SaveAs("../output/pdf/svd3dfit_"+proc+".pdf")

cnv = TCanvas("","",1000,1000)
cnv.Divide(2,2)
cnv.cd(1); histos["h_seed_resE"].Draw("hist")
cnv.cd(2); histos["h_seed_resPy"].Draw("hist")
cnv.cd(3); histos["h_seed_resE_vs_x"].Draw("col")
cnv.cd(4); histos["h_seed_resPy_vs_x"].Draw("col")
cnv.SaveAs("../output/pdf/seedsres_"+proc+".pdf")

cnv = TCanvas("","",1500,500)
cnv.Divide(3,1)
cnv.cd(1)
histos["h_N_sigacc"].SetLineColor(ROOT.kBlack); histos["h_N_sigacc"].Draw()
histos["h_N_all_seeds"].SetLineColor(ROOT.kBlue); histos["h_N_all_seeds"].Draw("same")
histos["h_N_matched_seeds"].SetLineColor(ROOT.kGreen); histos["h_N_matched_seeds"].Draw("same")
histos["h_N_good_seeds"].SetLineColor(ROOT.kRed); histos["h_N_good_seeds"].Draw("same")
cnv.cd(2)
histos["h_seeding_pool"].Draw("hist")
cnv.cd(3)
histos["h_seeding_score"].Draw("hist")
cnv.SaveAs("../output/pdf/seedsmult_"+proc+".pdf")

tF.cd()
tT.Write()
tF.Write()
tF.Close()





