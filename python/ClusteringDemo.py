#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import config as cfg
import ROOT
from ROOT import TH1D, TH2D, TCanvas, TPolyLine, TPolyMarker, TRandom

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.16)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")


### chip geometry:
x_translation = 67.73
y_translation = 0.61872
x_chip_sensitive = 29.94176
y_chip_sensitive = 13.76256
nx_pixels_sensitive = 1024 #30 #1024
ny_pixels_sensitive = 512 #15 #512
nhits = 40

colors = [1,2,4,6,9,28,15,46,38]

### histograms datastructure
histos = {}

### book histos
def Book():
   histos.update( { "h_hits" : TH2D("h_hits",";x [mm];y [mm];Hits",nx_pixels_sensitive,0,x_chip_sensitive, ny_pixels_sensitive,-y_chip_sensitive/2.,+y_chip_sensitive/2.) } )
   histos.update( { "h_edep" : TH1D("h_edep",";Energy deposition [keV];Hits",100,0,200) } )
   histos.update( { "h_performance" : TH1D("h_performance","Clustering performance;;Multiplicity",12,0,12) } )
   histos["h_performance"].GetXaxis().SetBinLabel(1, "All hits")
   histos["h_performance"].GetXaxis().SetBinLabel(2, "All clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(3, "1hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(4, "2hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(5, "3hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(6, "4hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(7, "5hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(8, "6hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(9, "7hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(10,"8hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(11,"9hit-clusters")
   histos["h_performance"].GetXaxis().SetBinLabel(12,"#geq10hit-clusters")

### fill histos
def FillRandom(nhits):
   rnd = TRandom()
   rnd.SetSeed()
   for n in range(nhits):
      xtrk = rnd.Uniform(0,x_chip_sensitive)
      ytrk = rnd.Uniform(-y_chip_sensitive/2,+y_chip_sensitive/2)
      histos["h_hits"].Fill(xtrk,ytrk)


### get data to histo
def FillFromHitsFile(filename):
   with open(filename) as fp:
      Lines = fp.readlines() 
      for line in Lines:
         if("#" in line): continue
         line = line.strip()
         words = line.split()
         ##bxNumber << hit_id << layer_id << det_id << e_dep << ev_weight << cell_x << cell_y
         if(int(words[2])==0 and int(words[3])==1000):
            xtrk = histos["h_hits"].GetXaxis().GetBinCenter( int(words[6]) )
            ytrk = histos["h_hits"].GetYaxis().GetBinCenter( int(words[7]) )
            edep_keV = float(words[4])*1e6
            histos["h_hits"].Fill(xtrk,ytrk)
            histos["h_edep"].Fill(edep_keV)


### get data to histo
def FillFromTracksFile(filename):
   with open(filename) as fp:
      Lines = fp.readlines() 
      for line in Lines:
         if("#" in line): continue
         line = line.strip()
         words = line.split()
         ### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z
         if(abs(int(words[1]))==11 and int(words[3])==1000 and float(words[4])<(x_chip_sensitive+x_translation) and float(words[5])<(y_chip_sensitive+y_translation)):
            xtrk = float(words[4])-x_translation ## x in chip's coordinate system
            ytrk = float(words[5])-y_translation ## y in chip's coordinate system
            histos["h_hits"].Fill(xtrk,ytrk)


### get live pixels
def LivePixels(h):
   pixels = []
   for bx in range(1,h.GetNbinsX()+1):
      for by in range(1,h.GetNbinsY()+1):
         hits = h.GetBinContent(bx,by)
         if(hits>0): pixels.append( {"x":bx,"y":by,"hits":hits} )
   return pixels


### get cluster position
def GeoPosition(cluster,h):
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
      posx,posy = GeoPosition(cluster,h)
      positions.append( [posx,posy] )
   return clusters,positions
   

### print cluster properties
def PrintClusters(clusters,positions,threshold=0):
   for i in range(len(clusters)):
      npixels = len(clusters[i])
      if(npixels>threshold):
         if(len(positions)>0): print("cluster with ",npixels,"pixels: ",clusters[i],"--> position =",positions[i])
         else:                 print("cluster with ",npixels,"pixels: ",clusters[i])
   

### get clusters contours
def GetClusterContour(cluster,h,col=ROOT.kRed):
   corners  = []
   for pixel in cluster:
      xdn = h.GetXaxis().GetBinLowEdge(pixel["x"])
      xup = h.GetXaxis().GetBinUpEdge(pixel["x"])
      ydn = h.GetYaxis().GetBinLowEdge(pixel["y"])
      yup = h.GetYaxis().GetBinUpEdge(pixel["y"])
      corners.append( [xdn,ydn] )
      corners.append( [xdn,yup] )
      corners.append( [xup,yup] )
      corners.append( [xup,ydn] )
      corners.append( [xdn,ydn] )
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
   contour = TPolyLine(n,x,y)
   contour.SetLineColor(col)
   contour.SetLineWidth(1)
   return contour
   
   
### get clusters contours
def GetClusterMarker(positions,h,col=ROOT.kRed):
   marker = TPolyMarker()
   marker.SetPoint(0,positions[0],positions[1])
   marker.SetMarkerColor(col)
   marker.SetMarkerStyle(24)
   return marker


### get all clusters contours
def GetContours(clusters,positions,h):
   contours = []
   markers  = []
   rnd = TRandom()
   rnd.SetSeed()
   for i in range(len(clusters)):
      cluster  = clusters[i]
      position = positions[i]
      icol     = int(rnd.Uniform(0,len(colors)))
      contour = GetClusterContour(cluster,h,colors[icol])
      marker  = GetClusterMarker(position,h,colors[icol])
      contours.append(contour)
      markers.append(marker)
   return contours,markers


##################################################################
##################################################################
##################################################################
### actually run
Book()
# FillRandom(nhits)
FillFromHitsFile("EBeamOnlyWIS_HitsInfo.txt")
# FillFromTracksFile("BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean.txt")
pixels = LivePixels(histos["h_hits"])
histos["h_performance"].Fill(0,len(pixels))
print("N_pixels before:",len(pixels))
clusters,positions = GetAllClusters(pixels,histos["h_hits"])
histos["h_performance"].Fill(1,len(clusters))
for cluster in clusters:
   if(len(cluster)==1):   histos["h_performance"].AddBinContent(2)
   elif(len(cluster)==2): histos["h_performance"].AddBinContent(3)
   elif(len(cluster)==3): histos["h_performance"].AddBinContent(4)
   elif(len(cluster)==4): histos["h_performance"].AddBinContent(5)
   elif(len(cluster)==5): histos["h_performance"].AddBinContent(6)
   elif(len(cluster)==6): histos["h_performance"].AddBinContent(7)
   elif(len(cluster)==7): histos["h_performance"].AddBinContent(8)
   elif(len(cluster)==8): histos["h_performance"].AddBinContent(9)
   elif(len(cluster)==9): histos["h_performance"].AddBinContent(10)
   else:                  histos["h_performance"].AddBinContent(11)
print("N_pixels after:",len(pixels))
PrintClusters(clusters,positions)
contours,markers = GetContours(clusters,positions,histos["h_hits"])

### draw
cnv = TCanvas("cnv","",1000,500)
cnv.SetTicks(1,1)
histos["h_hits"].Draw("colz")
for contour in contours: contour.Draw("same")
for marker  in markers:  marker.Draw("same")
cnv.SaveAs(storage+"/output/pdf/ClusteringDemo.pdf(")

cnv = TCanvas("cnv","",500,500)
cnv.SetTicks(1,1)
cnv.SetLogy()
histos["h_performance"].Draw("hist text0")
cnv.SaveAs(storage+"/output/pdf/ClusteringDemo.pdf")

cnv = TCanvas("cnv","",500,500)
cnv.SetTicks(1,1)
cnv.SetLogy()
histos["h_edep"].Draw("hist")
cnv.SaveAs(storage+"/output/pdf/ClusteringDemo.pdf)")



