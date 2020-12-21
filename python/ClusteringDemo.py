#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ast
import config as cfg
import ROOT
from ROOT import TH1D, TH2D, TCanvas, TPolyLine, TPolyMarker, TRandom, TPad

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
# ROOT.gStyle.SetPadBottomMargin(0.15)
# ROOT.gStyle.SetPadLeftMargin(0.16)

ROOT.gStyle.SetPadLeftMargin(0.02)
ROOT.gStyle.SetPadRightMargin(0.02)
ROOT.gStyle.SetLineWidth(1)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")


### chip geometry:
# x_translation = 67.73
# y_translation = 0.61872
x_chip_sensitive = 29.94176
y_chip_sensitive = 13.76256
nx_pixels_sensitive = 1024 #30 #1024
ny_pixels_sensitive = 512 #15 #512
nhits = 40

nlayers=8
nchips=9
nreallayers=nlayers #int(nlayers/2)
ntotalchips=nchips*2-1 ## one-chip overlap
ipads = { "0_1000":136,"0_1001":135,"0_1002":134,"0_1003":133,"0_1004":132,"0_1005":131,"0_1006":130,"0_1007":129,"0_1008":128,
          #        127          126          125          124          123          122          121          120          
          #        119          118          117          116          115          114          113          112          111
          "1_1000":111,"1_1001":110,"1_1002":109,"1_1003":108,"1_1004":107,"1_1005":106,"1_1006":105,"1_1007":104,"1_1008":103,
          
          "2_1000":102,"2_1001":101,"2_1002":100,"2_1003":99, "2_1004":98, "2_1005":97, "2_1006":96, "2_1007":95, "2_1008":94,
          #        93           92           91           90           89           88           87           86
          #        85           84           83           82           81           80           79           78           77
          "3_1000":77, "3_1001":76, "3_1002":75, "3_1003":74, "3_1004":73, "3_1005":72, "3_1006":71, "3_1007":70, "3_1008":69,

          "4_1000":68, "4_1001":67, "4_1002":66, "4_1003":65, "4_1004":64, "4_1005":63, "4_1006":62, "4_1007":61, "4_1008":60,
          #        59           58           57           56           55           54           53           52
          #        51           50           49           48           47           46           45           44           43
          "5_1000":43, "5_1001":42, "5_1002":41, "5_1003":40, "5_1004":39, "5_1005":38, "5_1006":37, "5_1007":36, "5_1008":35,

          "6_1000":34, "6_1001":33, "6_1002":32, "6_1003":31, "6_1004":30, "6_1005":29, "6_1006":28, "6_1007":27, "6_1008":26,
          #        25           24           23           22           21           20           19           18
          #        17           16           15           14           13           12           11           10            9
          "7_1000":9,  "7_1001":8,  "7_1002":7,  "7_1003":6,  "7_1004":5,  "7_1005":4,  "7_1006":3,  "7_1007":2,  "7_1008":1}

colors = [1,2,4,6,9,28,15,46,38]

### histograms datastructure
histos = {}


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
      histos["h_hits"].Fill(xtrk,ytrk)


def RejectPixelDueToVertex(trks_vtx):
   nbad = 0
   for vtx in trks_vtx:
      x = vtx[0]
      z = vtx[2]
      if(x<0 and (z>3600 and z<4600)): nbad +=  1
   return (nbad==len(trks_vtx))
   

### get data to histo
def FillFromHitsFile(layerid,detid,filename,isigBX=0):
   suf = suffix(layerid,detid)
   with open(filename) as fp:
      Lines = fp.readlines() 
      for line in Lines:
         if("#" in line): continue
         line = line.strip()
         words = line.split()
         ### bxNumber<<hit_id<<layer_id<<det_id<<e_dep<<ev_weight<<cell_x<<cell_y<<[pdgids]<<[energies]<<[trackids]<<[(vx1,vy1,vz1)]
         if(int(words[2])==layerid and int(words[3])==detid):
            bx = int(words[0])
            xtrk = histos["h_hits"+suf].GetXaxis().GetBinCenter( int(words[6]) )
            ytrk = histos["h_hits"+suf].GetYaxis().GetBinCenter( int(words[7]) )
            edep_keV = float(words[4])*1e6
            trks_pdg = words[8].replace("[","").replace("]","").split(",")
            trks_eng = words[9].replace("[","").replace("]","").split(",")
            trks_ids = words[10].replace("[","").replace("]","").split(",")
            trks_vtx = words[11].replace("(","[").replace(")","]")
            trks_pdg = list(map(int,   trks_pdg)) 
            trks_eng = list(map(float, trks_eng))
            trks_ids = list(map(int,   trks_ids))
            trks_vtx = ast.literal_eval(trks_vtx)
            ################################################################
            if(isigBX==0 and RejectPixelDueToVertex(trks_vtx)): continue ###
            ################################################################
            if(isigBX>0 and bx!=isigBX):          continue ### in a signal run, take only one BX data
            if(isigBX>0 and (1 not in trks_ids)): continue ### in a signal run, ignore pixels w/o signal tracks
            if(isigBX>0 and (1 in trks_ids)): histos["h_hit_flags"+suf].Fill(xtrk,ytrk)
            ################################################################
            ngam=0
            for i in range(len(trks_pdg)): ngam += 1 if(trks_pdg[i]==22) else 0
            trkEsum_keV=0
            for E in trks_eng: trkEsum_keV+=E*1e6
            # if(trkEsum_keV<edep_keV):
            #    print("Warning: trks energy sum < edep!")
            #    print("trks energy sum=",trkEsum_keV," edep=",edep_keV)
            diff_keV = trkEsum_keV-edep_keV
            histos["h_edep"+suf].Fill(edep_keV)
            histos["h_ntrksperpix"+suf].Fill(len(trks_pdg))
            histos["h_ntrksperpix_gam"+suf].Fill(ngam)
            histos["h_etrksedepdiff"+suf].Fill( diff_keV )
            histos["h_etrks_vs_etrksedepdiff"+suf].Fill( trkEsum_keV, diff_keV )
            histos["h_etrks_vs_etrksedepdiff_zoom"+suf].Fill( trkEsum_keV, diff_keV )
            histos["h_etrks_vs_etrksedepdiff_wide"+suf].Fill( trkEsum_keV, diff_keV )
            ## exclude pixels which have only photons:
            # if(11 not in trks_pdg and -11 not in trks_pdg): continue
            histos["h_hits"+suf].Fill(xtrk,ytrk)


# ### get data to histo
# def FillFromTracksFile(filename,layerid,detid):
#    suf = suffix(layerid,detid)
#    with open(filename) as fp:
#       Lines = fp.readlines()
#       for line in Lines:
#          if("#" in line): continue
#          line = line.strip()
#          words = line.split()
#          ### bxNumber << pdg << track_id << det_id << xx << yy << eneg << ev_weight << vtx_x << vtx_y << vtx_z
#          if(abs(int(words[1]))==11 and int(words[3])==1000 and float(words[4])<(x_chip_sensitive+x_translation) and float(words[5])<(y_chip_sensitive+y_translation)):
#             xtrk = float(words[4])-x_translation ## x in chip's coordinate system
#             ytrk = float(words[5])-y_translation ## y in chip's coordinate system
#             histos["h_hits"+suf].Fill(xtrk,ytrk)


### get live pixels
def LivePixels(h):
   pixels = []
   for bx in range(1,h.GetNbinsX()+1):
      for by in range(1,h.GetNbinsY()+1):
         hits = h.GetBinContent(bx,by)
         if(hits>0): pixels.append( {"x":bx,"y":by,"hits":hits} )
   return len(pixels),pixels


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
   return len(clusters),clusters,positions
   

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
   marker.SetMarkerSize(0.6 if(issig) else 0.08)
   marker.SetMarkerColor(ROOT.kBlack if(issig) else col)
   if(issig): print("hsig:",hsig.GetName(),"xy=",x,y)
   marker.SetPoint(0,x,y)
   return marker,issig


### get all clusters contours
def GetContours(clusters,positions,h,hsig=None):
   contours = []
   bkgmarkers = []
   sigmarkers = []
   rnd = TRandom()
   rnd.SetSeed()
   for i in range(len(clusters)):
      cluster  = clusters[i]
      position = positions[i]
      icol     = int(rnd.Uniform(0,len(colors)))
      col      = colors[icol]
      contour = GetClusterContour(cluster,h,col)
      # marker,issig = GetClusterMarker(position,h,col,hsig)
      marker,issig = GetClusterMarker(position,h,ROOT.kGray+1,hsig)
      contours.append(contour)
      if(issig): sigmarkers.append(marker)
      else:      bkgmarkers.append(marker)
   return contours,bkgmarkers,sigmarkers


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


def getcnvs(BXname):
   cnvx = 2000
   cnvy = 800
   gapx = 0.00001
   gapy = 0.005
   cnvs = {
      "hits":                        TCanvas("cnv_hits_"+BXname,"",cnvx,cnvy),
      "performance":                 TCanvas("cnv_performance_"+BXname,"",cnvx,cnvy),
      "edep":                        TCanvas("cnv_edep_"+BXname,"",cnvx,cnvy),
      "ntrksperpix":                 TCanvas("cnv_ntrksperpix_"+BXname,"",cnvx,cnvy),
      "etrksedepdiff":               TCanvas("cnv_etrksedepdiff_"+BXname,"",cnvx,cnvy),
      "etrks_vs_etrksedepdiff_zoom": TCanvas("cnv_etrks_vs_etrksedepdiff_zoom_"+BXname,"",cnvx,cnvy),
      "etrks_vs_etrksedepdiff":      TCanvas("cnv_etrks_vs_etrksedepdiff_"+BXname,"",cnvx,cnvy),
      "etrks_vs_etrksedepdiff_wide": TCanvas("cnv_etrks_vs_etrksedepdiff_wide_"+BXname,"",cnvx,cnvy),
      
   }
   for cname,cnv in cnvs.items(): cnv.Divide(ntotalchips,nreallayers,gapx,gapy)
   return cnvs


def getipad(layerid,detid):
   sid = str(layerid)+"_"+str(detid)
   ipad = ipads[sid]
   return ipad


def draw(layerid,detid,contours,bmarkers,smarkers,cnvs):
   suf = suffix(layerid,detid)
   ipad = getipad(layerid,detid)
   
   cnvs["hits"].cd(ipad)
   ROOT.gPad.SetTicks(0,0)
   # histos["h_hits"+suf].Draw("col")
   htmp = histos["h_hits"+suf].Clone("h_hits"+suf+"_tmp")
   htmp.Reset()
   htmp.DrawClone("col")
   # for contour in contours: contour.DrawClone("same")
   for bmarker in bmarkers: bmarker.DrawClone("same")
   for smarker in smarkers: smarker.DrawClone("same")
   print("ipad:",ipad,"nsigvls:",len(smarkers),"nbkgcls:",len(bmarkers))
   ROOT.gPad.RedrawAxis()
   
   cnvs["performance"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetLogy()
   histos["h_performance"+suf].Draw("hist text0")
   ROOT.gPad.RedrawAxis()

   cnvs["edep"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   histos["h_edep"+suf].Draw("hist")
   ROOT.gPad.RedrawAxis()

   cnvs["ntrksperpix"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetLogy()
   histos["h_ntrksperpix"+suf].SetFillColorAlpha(histos["h_ntrksperpix"+suf].GetLineColor(),0.5)
   histos["h_ntrksperpix"+suf].Draw("hist text0")
   histos["h_ntrksperpix_gam"+suf].SetLineColor(ROOT.kRed)
   histos["h_ntrksperpix_gam"+suf].SetFillColorAlpha(histos["h_ntrksperpix_gam"+suf].GetLineColor(),0.5)
   histos["h_ntrksperpix_gam"+suf].Draw("hist same text0")
   ROOT.gPad.RedrawAxis()
   
   cnvs["etrksedepdiff"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   ROOT.gPad.SetLogy()
   histos["h_etrksedepdiff"+suf].Draw("hist")
   ROOT.gPad.RedrawAxis()

   cnvs["etrks_vs_etrksedepdiff_zoom"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   histos["h_etrks_vs_etrksedepdiff_zoom"+suf].Draw("col")
   ROOT.gPad.RedrawAxis()
   
   cnvs["etrks_vs_etrksedepdiff"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   histos["h_etrks_vs_etrksedepdiff"+suf].Draw("col")
   ROOT.gPad.RedrawAxis()

   cnvs["etrks_vs_etrksedepdiff_wide"].cd(ipad)
   ROOT.gPad.SetTicks(1,1)
   histos["h_etrks_vs_etrksedepdiff_wide"+suf].Draw("col")
   ROOT.gPad.RedrawAxis()


def savecnvs(BXname):
   for cname,cnv in cnvs.items():
      cnv.SaveAs(storage+"/output/pdf/"+cname+"_"+BXname+".pdf")


##################################################################
##################################################################
##################################################################

for BX in [1,2,3,4]:
   ### define the canvases
   sBX = str(BX)
   BXname = "BX"+sBX
   cnvs = getcnvs(BXname)
   
   ### actually run
   for layerid in range(nlayers):
      for baredetid in range(nchips):
         detid = baredetid+1000
         print("------------- Starting BX #",BX," layer #",layerid,"and chip #",detid,"-------------")
         suf = suffix(layerid,detid)
         Book(layerid,detid)
         # FillRandom(layerid,detid,nhits)
         FillFromHitsFile(layerid,detid,"list_root_hics_signal_165gev_5000nm_WIS_HitsProperInfo.txt",BX) # --> SIG
         FillFromHitsFile(layerid,detid,"EBeamOnlyWIS_DividedByBX"+sBX+"_HitsProperInfo.txt") # --> BKG
         # FillFromTracksFile("BkgEBeam_SignalHics5000nmProvisional_BX1_trackInfoClean.txt")
         npix,pixels = LivePixels(histos["h_hits"+suf])
         ncls,clusters,positions = GetAllClusters(pixels,histos["h_hits"+suf])
         # PrintClusters(clusters,positions)
         FillNclustersHist(layerid,detid,npix,ncls)
         contours,bmarkers,smarkers = GetContours(clusters,positions,histos["h_hits"+suf],histos["h_hit_flags"+suf])
         draw(layerid,detid,contours,bmarkers,smarkers,cnvs)
   
   ### save the canvases
   savecnvs(BXname)