#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ast
import config as cfg
import Contour as cnt
import itertools
import ROOT
from ROOT import TH1D, TH2D, TCanvas, TPolyLine, TPolyMarker, TRandom, TPad, TLatex

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)
# ROOT.gStyle.SetPadLeftMargin(0.02)
# ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetLineWidth(1)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")


### chip geometry:
x_translation = {1000:67.73}
y_translation = {1000:0.61872}
x_chip_sensitive = 29.94176
y_chip_sensitive = 13.76256
nx_pixels_sensitive = 300 #9 #30 #1024
ny_pixels_sensitive = 150 #4 #15 #512
nhits = 4000 #5

pixelsize = (x_chip_sensitive/nx_pixels_sensitive)*(y_chip_sensitive/ny_pixels_sensitive)

### very important for the conventions
### that is, assuming the pixels are counted from 0 or from 1
### while anyway the assumption is that they are counted from the bottom left corner of the chip
countcellsfromzero = False ## in this simple study

xpixsize = x_chip_sensitive/nx_pixels_sensitive
ypixsize = y_chip_sensitive/ny_pixels_sensitive
pi = ROOT.TMath.Pi()

colors = [1,2,4,6,9,28,15,46,38]

### histograms datastructure
histos = {}


def GetGlobalCornerID(pixel,ncornerx,ncornery):
   ixpix = pixel["x"] if(countcellsfromzero) else pixel["x"]-1
   iypix = pixel["y"] if(countcellsfromzero) else pixel["y"]-1
   cornerid = ixpix + iypix*(nx_pixels_sensitive+1) + ncornerx + ncornery*(nx_pixels_sensitive+1)
   return cornerid


def GetAllGlobalCorners(h):
   text = {}
   for by in range(1,h.GetNbinsY()+1):
      for bx in range(1,h.GetNbinsX()+1):
         ix = bx-1 if(countcellsfromzero) else bx
         iy = by-1 if(countcellsfromzero) else by
         pixel = {"x":ix,"y":iy}
         cornerid = GetGlobalCornerID(pixel,0,0)
         cornerx = h.GetXaxis().GetBinLowEdge(bx)
         cornery = h.GetYaxis().GetBinLowEdge(by)
         text.update({cornerid:[cornerx,cornery]})
         # print(ix,iy,cornerid,[cornerx,cornery])
   return text


def DrawCornerIds(h):
   text = GetAllGlobalCorners(h)
   if(nx_pixels_sensitive*ny_pixels_sensitive>400): return
   for cid,cxy in text.items():
      s = TLatex()
      s.SetTextColor(ROOT.kRed)
      s.SetTextSize(0.015 if(nx_pixels_sensitive*ny_pixels_sensitive<10*5) else 0.01)
      s.DrawLatex(cxy[0],cxy[1],str(cid))


def suffix(layerid,detid):
   suf = "_L"+str(layerid)+"_C"+str(detid)
   return suf


### book histos
def Book(layerid,detid):
   suf = suffix(layerid,detid)
   histos.update( { "h_hit_flags"+suf                   : TH2D("h_hit_flags"+suf,";x [mm];y [mm];Is signal pixel",nx_pixels_sensitive,0,x_chip_sensitive, ny_pixels_sensitive,-y_chip_sensitive/2.,+y_chip_sensitive/2.) } )
   histos.update( { "h_hits"+suf                        : TH2D("h_hits"+suf,";x [mm];y [mm];Number of pixels",nx_pixels_sensitive,0,x_chip_sensitive, ny_pixels_sensitive,-y_chip_sensitive/2.,+y_chip_sensitive/2.) } )
   histos.update( { "h_performance"+suf                 : TH1D("h_performance"+suf,";;Multiplicity",12,0,12) } )
   histos.update( { "h_edep"+suf                        : TH1D("h_edep"+suf,";Energy deposition [keV];Number of pixels",100,0,100) } )
   histos.update( { "h_ntrksperpix"+suf                 : TH1D("h_ntrksperpix"+suf,";Particles per pixel;Number of pixels",7,0,7) } )
   histos.update( { "h_ntrksperpix_gam"+suf             : TH1D("h_ntrksperpix_gam"+suf,";Particles per pixel;Number of pixels",7,0,7) } )
   histos.update( { "h_etrksedepdiff"+suf               : TH1D("h_etrksedepdiff"+suf,";#SigmaE_{trks}-E_{dep} [keV];Number of pixels",100,0,5000) } )
   histos.update( { "h_etrks_vs_etrksedepdiff"+suf      : TH2D("h_etrks_vs_etrksedepdiff"+suf,";#SigmaE_{trks} [keV];#SigmaE_{trks}-E_{dep} [keV];Number of pixels",100,0,5000 ,100,0,5000) } )
   histos.update( { "h_etrks_vs_etrksedepdiff_zoom"+suf : TH2D("h_etrks_vs_etrksedepdiff_zoom"+suf,";#SigmaE_{trks} [keV];#SigmaE_{trks}-E_{dep} [keV];Number of pixels",100,0,500 ,100,0,500) } )
   histos.update( { "h_etrks_vs_etrksedepdiff_wide"+suf : TH2D("h_etrks_vs_etrksedepdiff_wide"+suf,";#SigmaE_{trks} [keV];#SigmaE_{trks}-E_{dep} [keV];Number of pixels",100,0,16000000 ,100,0,16000000) } )

   histos.update( { "h_contour_area"+suf                : TH1D("h_contour_area"+suf,";Cluster area/Pixel area;Number of clusters",100,0,10) } )
   
   histos["h_performance"+suf].GetXaxis().SetBinLabel(1, "All hits")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(2, "All clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(3, "1hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(4, "2hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(5, "3hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(6, "4hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(7, "5hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(8, "6hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(9, "7hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(10,"8hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(11,"9hit-clusters")
   histos["h_performance"+suf].GetXaxis().SetBinLabel(12,"#geq10hit-clusters")
   

### fill histos
def FillRandom(layerid,detid,nhits):
   rnd = TRandom()
   rnd.SetSeed()
   suf = suffix(layerid,detid)
   for n in range(nhits):
      xtrk = rnd.Uniform(0,x_chip_sensitive)
      ytrk = rnd.Uniform(-y_chip_sensitive/2,+y_chip_sensitive/2)
      histos["h_hits"+suf].Fill(xtrk,ytrk)


### get live pixels
def LivePixels(h):
   pixels = []
   for bx in range(1,h.GetNbinsX()+1):
      for by in range(1,h.GetNbinsY()+1):
         hits = h.GetBinContent(bx,by)
         if(hits>0): pixels.append( {"x":bx,"y":by,"hits":hits} )
   return len(pixels),pixels


### get the geometrical center position of the cluster
def GeoPosition(cluster,h):
   x = 0
   y = 0
   n = len(cluster)
   for pixel in cluster:
      x += h.GetXaxis().GetBinCenter( pixel["x"] )
      y += h.GetYaxis().GetBinCenter( pixel["y"] )
   return x/n, y/n


### get the center of gravity position of the cluster
def GrvPosition(cluster,h):
   x = 0
   y = 0
   n = 0
   for pixel in cluster:
      n += pixel["hits"]
      x += h.GetXaxis().GetBinCenter( pixel["x"] ) * pixel["hits"]
      y += h.GetYaxis().GetBinCenter( pixel["y"] ) * pixel["hits"]
   return x/n, y/n

      
### clusterisation recursion
def RecursiveClustering(cluster,pivot,pixels):
   if(pivot in pixels): 
      cluster.append(pivot) ## add this pixel to the cluster
      pixels.remove(pivot)  ## kill pixel from live pixels list
   for pixel in pixels[:]:
      dx = abs(pixel["x"]-pivot["x"])
      dy = abs(pixel["y"]-pivot["y"])
      if((dx+dy)<=2 and dx<=1 and dy<=1):
         nextpivot = pixel
         RecursiveClustering(cluster,nextpivot,pixels)


### get all clusters recursively
def GetAllClusters(pixels,h):
   clusters = []
   positions = []
   while len(pixels)>0: ## loop as long as there are live pixels in the list
      pixel = pixels[0] ## this is the pivot pixel for the cluster recursion
      cluster = []
      RecursiveClustering(cluster,pixel,pixels)
      clusters.append( cluster )
      posx,posy = GrvPosition(cluster,h)
      positions.append( [posx,posy] )
   return len(clusters),clusters,positions
   

### print cluster properties
def PrintClusters(clusters,positions,threshold=0):
   for i in range(len(clusters)):
      npixels = len(clusters[i])
      if(npixels>threshold):
         if(len(positions)>0): print("cluster with ",npixels,"pixels: ",clusters[i],"--> position =",positions[i])
         else:                 print("cluster with ",npixels,"pixels: ",clusters[i])


def GetPixCorners(pixel,h):
   corners = {}
   xdn = h.GetXaxis().GetBinLowEdge(pixel["x"])
   xup = h.GetXaxis().GetBinUpEdge(pixel["x"])
   ydn = h.GetYaxis().GetBinLowEdge(pixel["y"])
   yup = h.GetYaxis().GetBinUpEdge(pixel["y"])
   corners.update({"dd":[xdn,ydn]})
   corners.update({"du":[xdn,yup]})
   corners.update({"uu":[xup,yup]})
   corners.update({"ud":[xup,ydn]})
   return corners
   

def GetOuterCorners(cluster,h):
   allcorners = []
   # if(len(cluster)>2): print("CLuster size:",len(cluster))
   cornerids = {}
   cnames = {"dd":[0,0],"du":[0,1],"uu":[1,1],"ud":[1,0]}
   removedcornerids = []
   for pixel in cluster:
      pixcorners = GetPixCorners(pixel,h)
      for cname,corner in cnames.items():
         cornerid = GetGlobalCornerID(pixel,corner[0],corner[1])
         if(cornerid in cornerids):
            del cornerids[cornerid]           ## this is enoough if only 2 pixels share one corner. otherwise:
            removedcornerids.append(cornerid) ## need to recall this cornerid since 3 or 4 mixels share one corner
         else:
            if(cornerid not in removedcornerids):
               cornerids.update({cornerid:pixcorners[cname]})
   # print(GrvPosition(cluster,h),"--> cornerids keys:",cornerids.keys())
   # print(cornerids)
   xycorners = cornerids.values()
   # if(len(cluster)>2): print("Corners:",len(cluster),len(xycorners))
   return xycorners


def GetOuterContour(cluster,h):
   xpix,ypix  = GeoPosition(cluster,h)
   outcorners = GetOuterCorners(cluster,h)
   corners  = {}
   ## go clockwise
   for corner in outcorners:
      dx = corner[0]-xpix
      dy = corner[1]-ypix
      angle = -999
      if(dx!=0):
         absang = abs(math.atan(dy/dx))
         if(dx<0):
            if(dy>=0): angle = absang
            else:      angle = 2*pi-absang
         if(dx>0):
            if(dy>=0): angle = pi-absang
            else:      angle = pi+absang
      else:
         if(dy<0): angle = 3*pi/2
         else:     angle = +pi/2
      corners.update({angle:corner})
   cwcorners = dict(sorted(corners.items()))
   scwcorners = "["
   for txt,corner in cwcorners.items(): scwcorners += str(corner)+", "
   scwcorners += "]"
   # print("Sorted corners (for cluster size=",len(cluster),"): ",len(outcorners),scwcorners)
   srtcorners = []
   for angl,corner in cwcorners.items(): srtcorners.append(corner)
   contour = cnt.Contour(srtcorners,False)
   return srtcorners,contour


### get clusters contours
def GetClusterContour(cluster,h,col=ROOT.kRed):
   corners,contour = GetOuterContour(cluster,h)
   xlist = []
   ylist = []
   # for name,corner in corners.items():
   for corner in corners:
      xlist.append(corner[0])
      ylist.append(corner[1])
   xlist.append(xlist[0])
   ylist.append(ylist[0])
   x = array.array('d', xlist)
   y = array.array('d', ylist)
   n = len(x)
   contourline = TPolyLine(n,x,y)
   contourline.SetLineColor(col)
   contourline.SetLineWidth(1)
   return contourline,{"xcont":x,"ycont":y},contour


def GetClusterBareArea(cluster):
   barearea = 0
   for pix in cluster: barearea += xpixsize*ypixsize
   return barearea


def GetClusterShapeArea(contourdata):
   ## area = |Sum(x[i]*y[i+1]) - (y[i]*x[i+1])|/2
   ## with x[n+1] being x(1) and y[n+1] being y[1]
   shapearea = 0
   for i in range(len(contourdata["xcont"])):
      xi = contourdata["xcont"][i]
      yi = contourdata["ycont"][i]
      if(i<ncorners-1):
         xi1 = contourdata["xcont"][i+1]
         yi1 = contourdata["ycont"][i+1]
         shapearea += abs(xi*yi1-yi*xi1)/2
      else: 
         xi1 = contourdata["xcont"][0]
         yi1 = contourdata["ycont"][0]
         shapearea += abs(xi*yi1-yi*xi1)/2
   return shapearea


def GetClusterAxes(center,contourdata):
   elongation = 0
   dmaxp = -1e10
   dmaxn = -1e10
   imaxp = -1
   imaxn = -1
   for i in range(len(contourdata["xcont"])):
      x = contourdata["xcont"][i]
      y = contourdata["ycont"][i]
      dx = x-center[0]
      dy = y-center[1]
      d = math.sqrt(dx*dx + dy*dy)
      if(d>dmaxp and dx>0):
         dmaxp = d
         imaxp = i
      if(d>dmaxn and dx<0):
         dmaxn = d
         imaxn = i
   x = array.array('d', [ contourdata["xcont"][imin],contourdata["xcont"][imax] ])
   y = array.array('d', [ contourdata["ycont"][imin],contourdata["ycont"][imax] ])
   n = len(x)
   axis = TPolyLine(n,x,y)
   axis.SetLineColor(col)
   axis.SetLineWidth(1)
   return axis


def GetClusterProperties(cluster,center,contourdata):
   npix = len(cluster)
   barearea = GetClusterBareArea(cluster)
   
   ncorners = len(contourdata["xcont"])
   shapearea = GetClusterShapeArea(contourdata)
   
   elongation = 0
   
   prop = {"npix":npix, "barearea":barearea, "ncorners":ncorners, "shapearea":shapearea, }
   
   return npix,area,
   
   
### get clusters contours
def GetClusterMarker(position,h,col=ROOT.kRed,hsig=None):
   x = position[0]
   y = position[1]
   marker = TPolyMarker()
   issig = False
   if(hsig is not None):
      bx = hsig.GetXaxis().FindBin(x)
      by = hsig.GetYaxis().FindBin(y)
      if(hsig.GetBinContent(bx,by)>0): issig = True
   marker.SetMarkerStyle(20 if(issig) else 24)
   # marker.SetMarkerSize(0.6 if(issig) else 0.08)
   marker.SetMarkerSize(1)
   marker.SetMarkerColor(ROOT.kBlack if(issig) else col)
   if(issig): print("hsig:",hsig.GetName(),"xy=",x,y)
   marker.SetPoint(0,x,y)
   return marker,issig


### get all clusters contours
def GetContours(clusters,positions,h,hsig=None):
   contourspts = []
   contours    = []
   bkgmarkers  = []
   sigmarkers  = []
   rnd = TRandom()
   rnd.SetSeed()
   for i in range(len(clusters)):
      cluster  = clusters[i]
      position = positions[i]
      icol     = int(rnd.Uniform(0,len(colors)))
      col      = colors[icol]
      contourpts,contdata,contour = GetClusterContour(cluster,h,col)
      # marker,issig = GetClusterMarker(position,h,col,hsig)
      marker,issig = GetClusterMarker(position,h,ROOT.kGray+1,hsig)
      contourspts.append(contourpts)
      contours.append(contour)
      if(issig): sigmarkers.append(marker)
      else:      bkgmarkers.append(marker)
   return contourspts,bkgmarkers,sigmarkers,contours


def FillNclustersHist(layerid,detid,npix,ncls):
   suf = suffix(layerid,detid)
   histos["h_performance"+suf].Fill(0,npix)
   histos["h_performance"+suf].Fill(1,ncls)
   for cluster in clusters:
      npix = len(cluster)
      if(npix==1): histos["h_performance"+suf].Fill(2)
      if(npix==2): histos["h_performance"+suf].Fill(3)
      if(npix==3): histos["h_performance"+suf].Fill(4)
      if(npix==4): histos["h_performance"+suf].Fill(5)
      if(npix==5): histos["h_performance"+suf].Fill(6)
      if(npix==6): histos["h_performance"+suf].Fill(7)
      if(npix==7): histos["h_performance"+suf].Fill(8)
      if(npix==8): histos["h_performance"+suf].Fill(9)
      if(npix==9): histos["h_performance"+suf].Fill(10)
      if(npix>9):  histos["h_performance"+suf].Fill(11)


def FillContourHist(layerid,detid,contours):
   suf = suffix(layerid,detid)
   for contour in contours:
      histos["h_contour_area"+suf].Fill(contour.area/pixelsize)


def draw(layerid,detid,contourspts,bmarkers,smarkers,contours):
   suf = suffix(layerid,detid)
   
   cnv = TCanvas("hits","",1000,500)
   ROOT.gPad.SetTicks(1,1)
   histos["h_hits"+suf].Draw("col")
   if(nx_pixels_sensitive*ny_pixels_sensitive<400):
      for contourpts in contourspts: contourpts.DrawClone("same")
      for bmarker    in bmarkers:    bmarker.DrawClone("same")
   for smarker    in smarkers:     smarker.DrawClone("same")
   DrawCornerIds(histos["h_hits"+suf])
         
   # ROOT.gPad.RedrawAxis()
   cnv.SaveAs(storage+"/output/pdf/simpleclustering.pdf(")
   
   cnv = TCanvas("performance","",1000,500)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetLogy()
   histos["h_performance"+suf].Draw("hist text0")
   ROOT.gPad.RedrawAxis()
   cnv.SaveAs(storage+"/output/pdf/simpleclustering.pdf")
   
   cnv = TCanvas("performance","",1000,500)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetLogy()
   histos["h_contour_area"+suf].Draw("hist text0")
   ROOT.gPad.RedrawAxis()
   cnv.SaveAs(storage+"/output/pdf/simpleclustering.pdf)")


##################################################################
##################################################################
##################################################################
   
### actually run
layerid = 0
baredetid = 0
detid = baredetid+1000
suf = suffix(layerid,detid)
Book(layerid,detid)
FillRandom(layerid,detid,nhits)
npix,pixels = LivePixels(histos["h_hits"+suf])
ncls,clusters,positions = GetAllClusters(pixels,histos["h_hits"+suf])
PrintClusters(clusters,positions)
FillNclustersHist(layerid,detid,npix,ncls)
contourspts,bmarkers,smarkers,contours = GetContours(clusters,positions,histos["h_hits"+suf],histos["h_hit_flags"+suf])
FillContourHist(layerid,detid,contours)
draw(layerid,detid,contourspts,bmarkers,smarkers,contours)
