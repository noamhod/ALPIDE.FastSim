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
    def __init__(self,proc,tracksarr=[]):
        self.process = proc
        self.tracks = tracksarr
        self.geoManager = ROOT.TGeoManager("geoManager","Geometry")
        self.tube = None
        self.dipole = None
        self.stave = None
        self.chcell = None
        self.material  = None
        self.medium    = None
        self.idipole = 1
        self.ibeampipe = 2
        self.ifirststave = 100
        self.ifirstchcell = 1000
        self.ifirstcalolayer = 10000
        
    def createWorld(self):
        ### beampipe geometry
        beampipe_HalfLength = 200
        beampipe_Rmin = 3.5
        beampipe_Rmax = 4
        beampipe_zCenter = 200
        ### dipole geometry
        dipole_xHalfWidth = 120/2
        dipole_yHalfHeight = 67.2/2
        dipole_zHalfLength = 100/2
        dipole_zCenter = 152.9
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
        CherCell_xMin       = 4.00
        CherCell_xMax       = 54.0
        CherCell_yMin       = -1.
        CherCell_yMax       = +1.
        ### Calorimeter geometry
        CaloLayer_xOffset    = 5.7 if(self.process=="bppp") else 14
        CaloLayer_xHalfWidth = 27.5
        CaloLayer_yHalfWidth = 2.75
        CaloLayer_zHalfWidth = 0.394/2
        CaloLayer_zMin       = 350
        CaloLayer_zGap       = 0.02
        CaloLayer_xCenter    = CaloLayer_xOffset+CaloLayer_xHalfWidth
        CaloLayer_nLayers    = 20
        

        ### volumes        
        self.material  = TGeoMaterial("Al", 26.98,13,2.7)
        self.medium    = TGeoMedium("MED",1,self.material)
        world          = ROOT.gGeoManager.MakeBox("TOP",self.medium,10000,10000,10000)
        self.tube = ROOT.gGeoManager.MakeTube("beampipe",self.medium,beampipe_Rmin,beampipe_Rmax,beampipe_HalfLength)
        self.tube.SetLineColor(ROOT.kGray+2)
        self.dipole = ROOT.gGeoManager.MakeBox("dipole",self.medium,dipole_xHalfWidth,dipole_yHalfHeight,dipole_zHalfLength)
        self.dipole.SetLineColor(ROOT.kOrange+2)
        self.dipole.SetFillColor(ROOT.kOrange+2)
        self.stave = ROOT.gGeoManager.MakeBox("stave",self.medium,stave_xHalfWidth,stave_yHalfHeight,stave_zHalfLength)
        self.stave.SetLineColor(ROOT.kGreen+2)
        self.stave.SetFillColor(ROOT.kGreen+2)
        self.chcell = ROOT.gGeoManager.MakeBox("cherenkovcell",self.medium,CherCell_xHalfWidth,CherCell_yHalfWidth,CherCell_zHalfWidth)
        self.chcell.SetLineWidth(1)
        self.chcell.SetLineColor(ROOT.kCyan+2)
        self.chcell.SetFillColor(ROOT.kCyan+2)
        self.calolayer = ROOT.gGeoManager.MakeBox("calolayer",self.medium,CaloLayer_xHalfWidth,CaloLayer_yHalfWidth,CaloLayer_zHalfWidth)
        self.calolayer.SetLineWidth(1)
        self.calolayer.SetLineColor(ROOT.kBlue+2)
        self.calolayer.SetFillColor(ROOT.kBlue+2)
        
        ### add nodes to world
        # world.AddNodeOverlap(self.dipole,self.idipole,ROOT.TGeoTranslation(0,0,dipole_zCenter))
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
           x = CherCell_xMin
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
           
        world.AddNodeOverlap(self.tube,self.ibeampipe,ROOT.TGeoTranslation(0,0,beampipe_zCenter))
        
        return world

    def createTracks(self):
       for i in range(len(self.tracks)):
          track_index = ROOT.gGeoManager.AddTrack(i,11);
          track = ROOT.gGeoManager.GetTrack(track_index)
          track.SetLineWidth(1)
          track.SetLineColor(ROOT.kBlack)
          points = np.ndarray((self.tracks[i].GetN()*3), 'f', self.tracks[i].GetP()) ## GetP() returns a *float* *buffer* for GetN()*3!
          j = 0
          while j<(len(points)-3):
             x = points[j+0]
             y = points[j+1]
             z = points[j+2]
             track.AddPoint(x,y,z,0)
             j += 3

    def configureGeoManager(self,world):
        self.geoManager.SetTopVolume(world)
        self.geoManager.CloseGeometry()
        world.SetLineColor(8)
        self.geoManager.SetVisLevel(3)
        self.geoManager.SetVisOption(0)
        self.geoManager.SetTopVisible()

    def draw(self,world):
        # for trk in self.tracks: trk.Draw("same") ### this is TPolyLine3D
        self.createTracks()
        self.geoManager.DrawTracks("same")
        world.Draw("same")
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        
# if __name__ =='__main__':
#     geom = GeoLUXE()
#     world = geom.createWorld()
#     geom.configureGeoManager(world)
#     geom.draw(world)