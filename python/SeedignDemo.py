#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TH1D, TH2D, TH3D, TF1, TF2, TGraph, TGraph2D, TRandom, TVector2, TVector3, TLorentzVector, TPolyMarker3D, TPolyLine3D, TPolyLine, TCanvas, TView, TLatex, TLegend

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# ROOT.gErrorIgnoreLevel = ROOT.kWarning
ROOT.gErrorIgnoreLevel = ROOT.kError

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

###
Emax = 17.5 # GeV

### geometry:
zDipoleExit = 202.9
xDipoleExitMinAbs = 1  ## cm
xDipoleExitMaxAbs = 30 ## cm
xDipoleExit = (xDipoleExitMaxAbs-xDipoleExitMinAbs)/2.
yDipoleExitMin = -0.05 ## cm --> TODO: need tuning
yDipoleExitMax = +0.05 ## cm --> TODO: need tuning
xAbsMargins = 0.025 # cm --> TODO: need tuning
yAbsMargins = 0.025 # cm --> TODO: need tuning

### background stuff
NnoiseClusters = 100

layers = [1,2,3,4]

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

def getyzwindow(cluster,i):
   x1,y1,z1 = getpoint(cluster,i)
   particles = "electrons" if(x1>0) else "positrons"
   y1min = y1-yAbsMargins if((y1-yAbsMargins)>-0.75) else -0.75 ## must be within stave volume
   y1max = y1+yAbsMargins if((y1+yAbsMargins)<+0.75) else +0.75 ## must be within stave volume
   window_pts = TPolyMarker3D()
   window_pts.SetNextPoint(x1,y1min,z1)
   window_pts.SetNextPoint(x1,y1max,z1)
   if(particles=="electrons"): 
      window_pts.SetNextPoint(xDipoleExit,yDipoleExitMax,zDipoleExit)
      window_pts.SetNextPoint(xDipoleExit,yDipoleExitMin,zDipoleExit)
   else:
      window_pts.SetNextPoint(-xDipoleExit,yDipoleExitMax,zDipoleExit)
      window_pts.SetNextPoint(-xDipoleExit,yDipoleExitMin,zDipoleExit)
   window_lin = windowline(window_pts)
   window_lin.SetLineWidth(1)
   window_lin.SetLineColor(ROOT.kBlack)
   return window_pts,window_lin

def getwidewindow(cluster,i):
   x1,y1,z1 = getpoint(cluster,i)
   particles = "electrons" if(x1>0) else "positrons"
   x1min = 0
   x1max = 0
   if(particles=="electrons"): 
      x1min = x1-xAbsMargins if((x1-xAbsMargins)>4)  else 4  ## must be within stave volume
      x1max = x1+xAbsMargins if((x1+xAbsMargins)<31) else 31 ## must be within stave volume
   else:
      x1min = x1-xAbsMargins if((x1-xAbsMargins)>-31) else -31 ## must be within stave volume
      x1max = x1+xAbsMargins if((x1+xAbsMargins)<-4)  else -4  ## must be within stave volume
   window_pts = TPolyMarker3D()
   window_pts.SetNextPoint(x1min,y1,z1)                   ## top left
   window_pts.SetNextPoint(x1max,y1,z1)                   ## top right
   if(particles=="electrons"): 
      window_pts.SetNextPoint(xDipoleExitMaxAbs,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(xDipoleExitMinAbs,y1,zDipoleExit) ## bottom left
   else:
      window_pts.SetNextPoint(-xDipoleExitMinAbs,y1,zDipoleExit) ## bottom right
      window_pts.SetNextPoint(-xDipoleExitMinAbs,y1,zDipoleExit) ## bottom left
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
   
def trimwide(allpoints,points_wide,windowxz,windowyz,xpivot):
   for layer in layers:
      for i in range(allpoints["Cls"][layer].GetN()):
         acceptxz = isinwindowxz(allpoints["Cls"][layer],i,windowxz)
         acceptyz = isinwindowyz(allpoints["Cls"][layer],i,windowyz)
         if(not acceptxz): continue
         if(not acceptyz): continue
         r = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
         allpoints["Cls"][layer].GetPoint(i,r[0],r[1],r[2])
         if(isel(r[0]) and r[0]>1.001*xpivot):     continue
         if(not isel(r[0]) and r[0]<1.001*xpivot): continue
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
         # if(isel(r[0]) and r[0]>1.001*xpivot):     continue
         # if(not isel(r[0]) and r[0]<1.001*xpivot): continue
         points_narr["Cls"][layer].SetNextPoint(r[0],r[1],r[2])
         points_narr["IsSig"][layer].append(points_wide["IsSig"][layer][i])
         points_narr["TrkId"][layer].append(points_wide["TrkId"][layer][i])

def makeseed(cluster1,cluster2,i1,i2,particles):
   p = TLorentzVector()
   r1 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   r2 = [ROOT.Double(), ROOT.Double(), ROOT.Double()]
   cluster1.GetPoint(i1,r1[0],r1[1],r1[2])
   cluster2.GetPoint(i2,r2[0],r2[1],r2[2])

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
   
def draw(name,points,window_yz=None,window_xz=None):
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.ShowAxis()
   # view.SetRange(-80,-50,0, +80,+50,350)
   view.SetRange(0,-10,190, +30,+10,340)
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

def seed3dfit(name,r1,r2,r3,r4):
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
   g.GetZaxis().SetRangeUser(0,30)
   
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
   # print("residuals {xz,yz}={%g,%g}" % (residuals_xz,residuals_yz))

   # zz = np.linspace(290,340)
   zz = np.array([300,310,320,330])
   # xx,yy = lin(zz)
   xx,yy = line3d(zz, m_xz,c_xz,m_yz,c_yz)
   lfit = TPolyLine3D()
   for i in range(4):
      # print("[%g] x=%g, y=%g --> z=%g" % (i,xx[i],yy[i],zz[i]))
      lfit.SetNextPoint(zz[i],yy[i],xx[i])
   lfit.SetLineColor(ROOT.kBlue)
   
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.SetRange(290,-0.8,0, 340,+0.8,30)
   view.ShowAxis()
   g.Draw("p0")
   lfit.Draw("smae")
   cnv.SaveAs(name)   
   return residuals_xz,residuals_yz

def seed3dfitSVD(name,r1,r2,r3,r4):
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
   g.GetZaxis().SetRangeUser(0,30)
   
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
   linepts = vv[0] * np.mgrid[-20:20:2j][:, np.newaxis]

   # shift by the mean to get the line in the right place
   linepts += datamean
   
   lfit = TPolyLine3D()
   for point in linepts:
      lfit.SetNextPoint(point[2],point[1],point[0])
   lfit.SetLineColor(ROOT.kBlue)
   
   cnv = TCanvas("","",2000,2000)
   view = TView.CreateView(1)
   view.SetRange(290,-0.8,0, 340,+0.8,30)
   view.ShowAxis()
   g.Draw("p0")
   lfit.Draw("smae")
   cnv.SaveAs(name)   
   return dd ## a 1D array of the data singular values
   
   
def seed2dfit(name,r1,r2,r3,r4):
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
   gxyz.GetZaxis().SetRangeUser(0,30)
      
   gxz = TGraph()
   gxz.SetMarkerSize(1)
   gxz.SetMarkerStyle(24)
   gxz.SetMarkerColor(ROOT.kBlack)
   gxz.SetPoint(0,r1[2],r1[0])
   gxz.SetPoint(1,r2[2],r2[0])
   gxz.SetPoint(2,r3[2],r3[0])
   gxz.SetPoint(3,r4[2],r4[0])
   gxz.GetXaxis().SetRangeUser(290,340)
   gxz.GetYaxis().SetRangeUser(0,+30)
   
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
   view.SetRange(290,-0.8,0, 340,+0.8,30)
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
fitsEx = GetFits()

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
pxSeed = ROOT.std.vector( float )()
pySeed = ROOT.std.vector( float )()
pzSeed = ROOT.std.vector( float )()
eSeed  = ROOT.std.vector( float )()
pxGen  = ROOT.std.vector( float )()
pyGen  = ROOT.std.vector( float )()
pzGen  = ROOT.std.vector( float )()
eGen   = ROOT.std.vector( float )()
iGen   = ROOT.std.vector( int )()
tF = TFile("../data/root/seeds.root","RECREATE")
tF.cd()
tT = TTree("seeds","seeds")
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
           
           "h_seed_resE":  TH1D("seed_resE", ";(E_{seed}-E_{gen})/E_{gen};Tracks",    100,-3,+3),
           "h_seed_resPz": TH1D("seed_resPz",";(Pz_{seed}-Pz_{gen})/Pz_{gen};Tracks", 100,-3,+3), 
           "h_seed_resPy": TH1D("seed_resPy",";(Py_{seed}-Py_{gen})/Py_{gen};Tracks", 100,-5,+5),
           
           "h_N_signal":             TH1D("N_signal",             ";Track multiplicity;Events", 75,0,150),
           "h_N_all_seeds":          TH1D("N_all_seeds",          ";Track multiplicity;Events", 75,0,150),
           "h_N_matched_seeds":      TH1D("N_matched_seeds",      ";Track multiplicity;Events", 75,0,150),
           "h_N_good_seeds":         TH1D("N_good_seeds",         ";Track multiplicity;Events", 75,0,150),
}

pdfname = "../output/pdf/seedingdemo.pdf"
intfile = TFile("../data/root/rec_bppp.root","READ")
intree = intfile.Get("res")
nevents = intree.GetEntries()
print("with %d events" % nevents)
nmax = 10000
n=0 ### init n
for event in intree:
   Nsignal = 0
   Nseeds = 0
   Nmatched = 0
   Ngood = 0
   
   ## clear the output vectors
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
   pxSeed.clear()
   pySeed.clear()
   pzSeed.clear()
   eSeed.clear()
   pxGen.clear()
   pyGen.clear()
   pzGen.clear()
   eGen.clear()
   iGen.clear()
   
   ### start the loop
   if(n>nmax): break
   
   ### container for all clusters
   allpoints = initpoints()
   
   ## clusters' vectors are always written out (even if empty) for all gen tracks!
   ## each generated track in the vector always has 4 clusters accessed via TPolyMarker3D::GetPoint()
   # print(" NAllSigTrks=",event.polm_clusters.size())
   
   for i in range(event.polm_clusters.size()):
      
      #############################
      if(event.qgen[i]>0): continue
      #############################
      
      Nsignal += 1
      
      wgt  = event.wgtgen[i]
      pgen = event.pgen[i]
      ### cut on acceptance
      if(event.acctrkgen[i]!=1): continue 
      
      ### write the truth track momentum and its index
      pxGen.push_back(pgen.Px())
      pyGen.push_back(pgen.Py())
      pzGen.push_back(pgen.Pz())
      eGen.push_back(pgen.E())
      iGen.push_back(i)
      # print("Genr[%g]: E=%g, pT=%g, eta=%g, phi=%g, theta=%g" % (i,pgen.E(),pgen.Pt(),pgen.Eta(),pgen.Phi(),pgen.Theta()) )
      
      ### loop over all clusters of the track and put in the allpoints classified by the layer
      for jxy in range(event.polm_clusters[i].GetN()):
         rcls = [ ROOT.Double(), ROOT.Double(), ROOT.Double() ]
         event.polm_clusters[i].GetPoint(jxy,rcls[0],rcls[1],rcls[2]) ### the clusters
         AddPoint(allpoints,rcls,True,i)
   Nsig4 = getNnon0(allpoints["Cls"][4])
      
   ### embed some noise clusters
   rnd = TRandom()
   rnd.SetSeed()
   for kN in range(NnoiseClusters):
      for layer in layers:
         # for side in [0,1]:
         for side in [1]:
            x = 0
            if(side==0): x = rnd.Uniform(-31,-4)
            if(side==1): x = rnd.Uniform(4,31)
            y = rnd.Uniform(-0.75,+0.75)
            if(layer==1): z = 300
            if(layer==2): z = 310
            if(layer==3): z = 320
            if(layer==4): z = 330
            rnoise = [x,y,z]
            AddPoint(allpoints,rnoise)
   Nbkgsig4 = getNnon0(allpoints["Cls"][4])
   # print(" Nbkgsig4=",Nbkgsig4)
   
   
   draw(pdfname+"(",allpoints)
   
   
   Nall4 = Nbkgsig4   
   ### loop over the clusters and start the seeding
   # print(" Nall4=",Nall4)
   # for j4 in range(allpoints["Cls"][4].GetN()):
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
      # print("  Nwide1=",Nwide1)
      
      draw(pdfname,widepoints,winlin_yz,winlin_xz_wide)
      
      # print("  Nwide1=",widepoints["Cls"][1].GetN())
      ### choose one cluster in layer 1 as the second seed
      # for j1 in range(widepoints["Cls"][1].GetN()):
      # if(Nwide1<1): print("failed [sig?="+str(allpoints["IsSig"][4][j4])+"]: with Nwide1<1")
      for j1 in range(Nwide1):
         ### get the narrow window  (window corners must be added clockwise!)
         winpts_xz_narr,winlin_xz_narr = getnarrwindow(allpoints["Cls"][4],widepoints["Cls"][1],j4,j1)
         ### discard all clusters which are not in the narrow window
         narrpoints = initpoints()
         trimnarr(widepoints,narrpoints,winpts_xz_narr)
         Nnarr2 = getNnon0(narrpoints["Cls"][2])
         Nnarr3 = getNnon0(narrpoints["Cls"][3])
         
         # print("   Nnarr1=",Nnarr1)
         # print("   Nnarr2=",Nnarr2)
         # print("   Nnarr3=",Nnarr3)
         
         # if(Nnarr2<1): print("failed: with Nnarr2<1")
         # if(Nnarr3<1): print("failed: with Nnarr3<1")
         
         ### check if there are at least 1 cluster in both layer 2 and layer 3 within the narrow window
         if(Nnarr2<1): continue
         if(Nnarr3<1): continue
         
         draw(pdfname,narrpoints,None,winlin_xz_narr)
         
         ### get the seed - note: could be that there are several combinations but the seed momentum would be identical
         pseed = makeseed(allpoints["Cls"][4],widepoints["Cls"][1],j4,j1,particlename)
         if(pseed.E()>Emax): continue
         
         r1 = getpoint(widepoints["Cls"][1],j1)
         
         NseedsPerWindow = Nnarr2*Nnarr3
         # print("NseedsPerWindow=",NseedsPerWindow)
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
               
               # print("r1=",r1)
               # print("r2=",r2)
               # print("r3=",r3)
               # print("r4=",r4)
               # print("Seed: IsSig=%g, TrkIds={%g,%g,%g,%g} --> E=%g, pT=%g, eta=%g, phi=%g, theta=%g" % (issig,t1,t2,t3,t4,pseed.E(),pseed.Pt(),pseed.Eta(),pseed.Phi(),pseed.Theta()) )

               ### two independent 2d fits
               chi2_xz,prob_xz,chi2_yz,prob_yz = seed2dfit(pdfname,r1,r2,r3,r4)
               # print("chi2_xz=%g, prob_xz=%g, chi2_yz=%g, prob_yz=%g" % (chi2_xz,prob_xz,chi2_yz,prob_yz))
               if(issiguniq):
                  histos["h_chi2ndf_xz_sig"].Fill(chi2_xz)
                  histos["h_chi2ndf_yz_sig"].Fill(chi2_yz)
               else:
                  histos["h_chi2ndf_xz_bkg"].Fill(chi2_xz)
                  histos["h_chi2ndf_yz_bkg"].Fill(chi2_yz)

               ### a single 3d fit
               res_xz, res_yz = seed3dfit(pdfname,r1,r2,r3,r4)
               if(issiguniq):
                  histos["h_residuals_xz_sig"].Fill(res_xz)
                  histos["h_residuals_yz_sig"].Fill(res_yz)
               else:
                  histos["h_residuals_xz_bkg"].Fill(res_xz)
                  histos["h_residuals_yz_bkg"].Fill(res_yz)
                  
               ### a single 3d fit SVD
               dd = seed3dfitSVD(pdfname,r1,r2,r3,r4)
               # print("dd=",dd)
               if(issiguniq):
                  histos["h_svd_dd0_sig"].Fill(dd[0])
                  histos["h_svd_dd1_sig"].Fill(dd[1])
                  histos["h_svd_dd2_sig"].Fill(dd[2])
               else:
                  histos["h_svd_dd0_bkg"].Fill(dd[0])
                  histos["h_svd_dd1_bkg"].Fill(dd[1])
                  histos["h_svd_dd2_bkg"].Fill(dd[2])
               
               isgood = (dd[1]<0.005 and dd[2]<0.0025)
               if(not isgood): continue
               Ngood += 1
               
               pgen = TLorentzVector()
               igen = -1
               for k in range(iGen.size()):
                  if(iGen[k]==j4):
                     igen = k
                     break
               pgen.SetPxPyPzE(pxGen[igen],pyGen[igen],pzGen[igen],eGen[igen])
               resE = (pseed.E()-pgen.E())/pgen.E()
               resPz = (pseed.Pz()-pgen.Pz())/pgen.Pz()
               resPy = (pseed.Py()-pgen.Py())/pgen.Py()
               histos["h_seed_resE"].Fill(resE)
               histos["h_seed_resPz"].Fill(resPz)
               histos["h_seed_resPy"].Fill(resPy)
               
               ### write out the good seeds
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
               
   histos["h_N_signal"].Fill(Nsignal)       
   histos["h_N_all_seeds"].Fill(Nseeds)       
   histos["h_N_matched_seeds"].Fill(Nmatched)       
   histos["h_N_good_seeds"].Fill(Ngood)       
   
   cnv = TCanvas("","",2000,2000)
   cnv.SaveAs(pdfname+")")
   print("Event: %g --> Nsignal=%g, Nseeds=%g, Nmatched=%g, Ngood=%g" % (n,Nsignal,Nseeds,Nmatched,Ngood))   
   # quit()
   tT.Fill()
   if(n%10==0 and n>0): print("  processed %d events" % n)
   n+=1
print("Total events processed: ",n)

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_chi2ndf_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_xz_sig"].Draw()
histos["h_chi2ndf_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_xz_bkg"].Draw("same")
cnv.cd(2)
histos["h_chi2ndf_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_chi2ndf_yz_sig"].Draw()
histos["h_chi2ndf_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_chi2ndf_yz_bkg"].Draw("same")
cnv.SaveAs("../output/pdf/chi2ndf.pdf")

cnv = TCanvas("","",1000,500)
cnv.Divide(2,1)
cnv.cd(1)
histos["h_residuals_xz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_xz_sig"].Draw()
histos["h_residuals_xz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_xz_bkg"].Draw("same")
cnv.cd(2)
histos["h_residuals_yz_sig"].SetLineColor(ROOT.kRed);   histos["h_residuals_yz_sig"].Draw()
histos["h_residuals_yz_bkg"].SetLineColor(ROOT.kBlack); histos["h_residuals_yz_bkg"].Draw("same")
cnv.SaveAs("../output/pdf/res.pdf")

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
cnv.SaveAs("../output/pdf/svd.pdf")

cnv = TCanvas("","",1500,500)
cnv.Divide(3,1)
cnv.cd(1); histos["h_seed_resE"].Draw("hist")
cnv.cd(2); histos["h_seed_resPz"].Draw("hist")
cnv.cd(3); histos["h_seed_resPy"].Draw("hist")
cnv.SaveAs("../output/pdf/seedsres.pdf")

cnv = TCanvas("","",1000,1000)
histos["h_N_signal"].SetLineColor(ROOT.kBlack); histos["h_N_signal"].Draw()
histos["h_N_all_seeds"].SetLineColor(ROOT.kBlue); histos["h_N_all_seeds"].Draw("same")
histos["h_N_matched_seeds"].SetLineColor(ROOT.kGreen); histos["h_N_matched_seeds"].Draw("same")
histos["h_N_good_seeds"].SetLineColor(ROOT.kRed); histos["h_N_good_seeds"].Draw("same")
cnv.SaveAs("../output/pdf/seedmult.pdf")

tF.cd()
tT.Write()
tF.Write()
tF.Close()





