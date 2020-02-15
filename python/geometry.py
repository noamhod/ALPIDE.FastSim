#!/usr/bin/python
import os
import math
import subprocess
import array
import numpy as np
import ROOT
from ROOT import TFile, TPolyLine3D, TGeoTube, TGeoManager, TGeoMaterial, TGeoMedium, TGeoVolume, TGeoTranslation, TVirtualGeoTrack, TView

# def staves(proc):
#    tfile = TFile("../data/root/"+proc+"_geometry.root","READ")
#    stvs = [ tfile.Get("TPolyLine3D;9"), tfile.Get("TPolyLine3D;8"),
#               tfile.Get("TPolyLine3D;7"), tfile.Get("TPolyLine3D;6"),
#               tfile.Get("TPolyLine3D;5"), tfile.Get("TPolyLine3D;4"),
#               tfile.Get("TPolyLine3D;3"), tfile.Get("TPolyLine3D;2")]
#    return stvs
#
# def dipole(proc):
#    tfile = TFile("../data/root/"+proc+"_geometry.root","READ")
#    dipl = tfile.Get("TPolyLine3D;1")
#    return dipl
   
# def beampipe():
#    manager   = TGeoManager("beampipe", "poza2")
#    material  = TGeoMaterial("Al", 26.98,13,2.7)
#    medium    = TGeoMedium("MED",1,material)
#    topvolume = ROOT.gGeoManager.MakeBox("TOP",medium,100,100,100)
#    ROOT.gGeoManager.SetTopVolume(topvolume)
#    volume    = ROOT.gGeoManager.MakeTube("TUBE",medium, 3.5,4,115);
#    volume.SetLineWidth(2)
#    topvolume.AddNode(volume,1,TGeoTranslation(0,0,315))
#    ROOT.gGeoManager.CloseGeometry()
#    # ROOT.gGeoManager.SetNsegments(80)
#    return topvolume

class GeoLUXE():
    def __init__(self,proc,stracksarr=[],btracksarr=[]):
        self.process = proc
        self.stracks = stracksarr
        self.btracks = btracksarr
        self.geoManager = ROOT.TGeoManager("geoManager","Geometry")
        self.tube = None
        self.dipole = None
        self.stave = None
        self.chcell = None
        self.material  = None
        self.medium    = None
        self.ibeampipe = 1
        self.idipole = 10
        self.ifirststave = 100
        self.ifirstchcell = 1000
        self.ifirstcalolayer = 10000
        self.ifirststrack = 100000
        self.ifirstbtrack = 1000000
        
    def createWorld(self):
        ### beampipe geometry
        beampipe_HalfLength = 200
        beampipe_Rmin = 3.6
        beampipe_Rmax = 4
        beampipe_zCenter = 200
        ### dipole geometry (full block)
        dipole_xHalfWidth = 120/2
        dipole_yHalfHeight = 67.2/2
        dipole_zHalfLength = 100/2
        dipole_zCenter = 152.9
        ### dipole geometry (assembly)
        dipole_horizontal_plate_xHalfWidth = 120/2
        dipole_horizontal_plate_yHalfHeight = (67.2-beampipe_Rmax*10)/2/2
        dipole_horizontal_plate_zHalfLength = 100/2
        dipole_horizontal_plate_zCenter = 152.9
        dipole_horizontal_plate_yAbsCenter = (beampipe_Rmax*10)/2+dipole_horizontal_plate_yHalfHeight
        
        dipole_vertical_plate_xHalfWidth = (120-beampipe_Rmax*10)/2/2
        dipole_vertical_plate_yHalfHeight = (beampipe_Rmax*10)/2
        dipole_vertical_plate_zHalfLength = 100/2
        dipole_vertical_plate_zCenter = 152.9
        dipole_vertical_plate_xAbsCenter = (beampipe_Rmax*10)/2+dipole_vertical_plate_xHalfWidth
        
        ### stave geometry
        stave_xHalfWidth = 27/2 if(self.process=="bppp") else 50/2
        stave_yHalfHeight = 1.5/2
        stave_zHalfLength = 0.01/2
        stave_z1 = 300
        stave_z2 = 310
        stave_z3 = 320
        stave_z4 = 330
        stave_xOffset = 5.7 if(self.process=="bppp") else 14
        stave_xCenter = stave_xOffset+stave_xHalfWidth
        ### Cherencov geometry
        CherCell_xHalfWidth = 1.0
        CherCell_yHalfWidth = 1.0
        CherCell_zHalfWidth = 1.0
        CherCell_zCenter    = 330
        CherCell_xGap       = 0.1
        CherCell_yGap       = 0.1
        CherCell_xMin       = 4.5
        CherCell_xMax       = 54.5
        CherCell_yMin       = -1.
        CherCell_yMax       = +1.
        ### Calorimeter geometry
        CaloLayer_xOffset    = 5.7 if(self.process=="bppp") else 14
        CaloLayer_xHalfWidth = 27.5
        CaloLayer_yHalfWidth = 2.75
        CaloLayer_zHalfWidth = 0.394/2
        CaloLayer_zMin       = 340
        CaloLayer_zGap       = 0.02
        CaloLayer_xCenter    = CaloLayer_xOffset+CaloLayer_xHalfWidth
        CaloLayer_nLayers    = 20
        

        ### volumes        
        # self.material  = TGeoMaterial("Al", 26.98,13,2.7)
        self.material  = TGeoMaterial("Vac", 0,0,0)
        self.medium    = TGeoMedium("MED",1,self.material)
        world          = ROOT.gGeoManager.MakeBox("TOP",self.medium,10000,10000,10000)
        
        # self.dipole = ROOT.gGeoManager.MakeBox("dipole",self.medium,dipole_xHalfWidth,dipole_yHalfHeight,dipole_zHalfLength)
        # self.dipole.SetLineWidth(1)
        # self.dipole.SetLineColor(ROOT.kOrange+2)
        
        self.dipole_hor = ROOT.gGeoManager.MakeBox("dipole_hor",self.medium,dipole_horizontal_plate_xHalfWidth,dipole_horizontal_plate_yHalfHeight,dipole_horizontal_plate_zHalfLength)
        self.dipole_hor.SetLineWidth(1)
        self.dipole_hor.SetLineColor(ROOT.kOrange+2)
        
        self.dipole_ver = ROOT.gGeoManager.MakeBox("dipole_ver",self.medium,dipole_vertical_plate_xHalfWidth,dipole_vertical_plate_yHalfHeight,dipole_vertical_plate_zHalfLength)
        self.dipole_ver.SetLineWidth(1)
        self.dipole_ver.SetLineColor(ROOT.kOrange+2)
        
        self.tube = ROOT.gGeoManager.MakeTube("beampipe",self.medium,beampipe_Rmin,beampipe_Rmax,beampipe_HalfLength)
        self.tube.SetLineWidth(1)
        self.tube.SetLineColor(ROOT.kGray)
        self.stave = ROOT.gGeoManager.MakeBox("stave",self.medium,stave_xHalfWidth,stave_yHalfHeight,stave_zHalfLength)
        self.stave.SetLineWidth(1)
        self.stave.SetLineColor(ROOT.kGreen+2)
        self.chcell = ROOT.gGeoManager.MakeBox("cherenkovcell",self.medium,CherCell_xHalfWidth,CherCell_yHalfWidth,CherCell_zHalfWidth)
        self.chcell.SetLineWidth(1)
        self.chcell.SetLineColor(ROOT.kCyan+2)
        self.calolayer = ROOT.gGeoManager.MakeBox("calolayer",self.medium,CaloLayer_xHalfWidth,CaloLayer_yHalfWidth,CaloLayer_zHalfWidth)
        self.calolayer.SetLineWidth(1)
        self.calolayer.SetLineColor(ROOT.kBlue-1)
        
        ### add nodes to world
        # world.AddNodeOverlap(self.dipole,self.idipole,ROOT.TGeoTranslation(0,0,dipole_zCenter))
        world.AddNodeOverlap(self.dipole_hor,self.idipole,  ROOT.TGeoTranslation(0,+dipole_horizontal_plate_yAbsCenter,dipole_horizontal_plate_zCenter))
        world.AddNodeOverlap(self.dipole_hor,self.idipole+1,ROOT.TGeoTranslation(0,-dipole_horizontal_plate_yAbsCenter,dipole_horizontal_plate_zCenter))
        world.AddNodeOverlap(self.dipole_ver,self.idipole+2,ROOT.TGeoTranslation(+dipole_vertical_plate_xAbsCenter,0,dipole_horizontal_plate_zCenter))
        world.AddNodeOverlap(self.dipole_ver,self.idipole+3,ROOT.TGeoTranslation(-dipole_vertical_plate_xAbsCenter,0,dipole_horizontal_plate_zCenter))

        world.AddNodeOverlap(self.tube,self.ibeampipe,ROOT.TGeoTranslation(0,0,beampipe_zCenter))
        
        world.AddNodeOverlap(self.stave,self.ifirststave+0,ROOT.TGeoTranslation(-stave_xCenter,0,stave_z1))
        world.AddNodeOverlap(self.stave,self.ifirststave+1,ROOT.TGeoTranslation(-stave_xCenter,0,stave_z2))
        world.AddNodeOverlap(self.stave,self.ifirststave+2,ROOT.TGeoTranslation(-stave_xCenter,0,stave_z3))
        world.AddNodeOverlap(self.stave,self.ifirststave+3,ROOT.TGeoTranslation(-stave_xCenter,0,stave_z4))
        if(self.process=="bppp"):
           world.AddNodeOverlap(self.stave,self.ifirststave+4,ROOT.TGeoTranslation(+stave_xCenter,0,stave_z1))
           world.AddNodeOverlap(self.stave,self.ifirststave+5,ROOT.TGeoTranslation(+stave_xCenter,0,stave_z2))
           world.AddNodeOverlap(self.stave,self.ifirststave+6,ROOT.TGeoTranslation(+stave_xCenter,0,stave_z3))
           world.AddNodeOverlap(self.stave,self.ifirststave+7,ROOT.TGeoTranslation(+stave_xCenter,0,stave_z4))
        ### add Cherenkov
        if(self.process=="trident"):
           n = 0
           x = CherCell_xMin+CherCell_xHalfWidth
           while x<CherCell_xMax:
              y = CherCell_yMin
              while y<CherCell_yMax:
                 world.AddNodeOverlap(self.chcell,self.ifirstchcell+n,ROOT.TGeoTranslation(x,y,CherCell_zCenter))
                 n += 1
                 y += CherCell_yGap+CherCell_yHalfWidth
              x += CherCell_xGap+CherCell_xHalfWidth
        ### add calo
        for j in range(CaloLayer_nLayers):
           world.AddNodeOverlap(self.calolayer,self.ifirstcalolayer+j,ROOT.TGeoTranslation(-CaloLayer_xCenter,0,CaloLayer_zMin+j*(CaloLayer_zGap+2*CaloLayer_zHalfWidth)))
           if(self.process=="bppp"): world.AddNodeOverlap(self.calolayer,self.ifirstcalolayer+(j+100),ROOT.TGeoTranslation(+CaloLayer_xCenter,0,CaloLayer_zMin+j*(CaloLayer_zGap+2*CaloLayer_zHalfWidth)))
        
        return world

    def createTracks(self):
       for i in range(len(self.stracks)):
          geotrack_index = self.geoManager.AddTrack(self.ifirststrack+i,11);
          geotrack = self.geoManager.GetTrack(geotrack_index)
          geotrack.SetLineWidth(1)
          geotrack.SetLineColor(ROOT.kBlack)
          points = np.ndarray((self.stracks[i].GetN()*3), 'f', self.stracks[i].GetP()) ## GetP() returns a *float* *buffer* for GetN()*3!
          j = 0
          while j<(len(points)):
             geotrack.AddPoint(points[j+0],points[j+1],points[j+2],0)
             j += 3
       for i in range(len(self.btracks)):
          geotrack_index = self.geoManager.AddTrack(self.ifirstbtrack+i,11);
          geotrack = self.geoManager.GetTrack(geotrack_index)
          geotrack.SetLineWidth(1)
          geotrack.SetLineColor(ROOT.kRed)
          points = np.ndarray((self.btracks[i].GetN()*3), 'f', self.btracks[i].GetP()) ## GetP() returns a *float* *buffer* for GetN()*3!
          j = 0
          while j<(len(points)):
             if(points[j+2]>=200 and points[j+2]<=360):
                if(points[j+2]>points[(j-3)+2]): geotrack.AddPoint(points[j+0],points[j+1],points[j+2],0)
             j += 3
             

    def configureGeoManager(self,world):
        self.createTracks()
        self.geoManager.SetTopVolume(world)
        self.geoManager.CloseGeometry()
        self.geoManager.GetTopVolume().Raytrace()
        for trk in self.stracks: trk.Draw("same")
        for trk in self.btracks: trk.Draw("same")
        # self.geoManager.SetVisLevel(3)
        # self.geoManager.SetVisOption(0)
        # self.geoManager.SetTopVisible()
        # self.geoManager.CheckOverlaps()
        self.geoManager.Export("../output/root/GeoLUXE_"+self.process+".root")
        

    def draw(self,world):
        # for trk in self.stracks: trk.Draw("same") ### this is TPolyLine3D
        self.geoManager.DrawTracks("same")
        world.Draw("same")
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        
# if __name__ =='__main__':
#     geom = GeoLUXE()
#     world = geom.createWorld()
#     geom.configureGeoManager(world)
#     geom.draw(world)