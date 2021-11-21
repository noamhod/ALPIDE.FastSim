#!/usr/bin/python
import os
import math, time
import h5py
import subprocess
from array import array
import numpy as np
import ROOT
from ROOT import TFile, TTree, TLorentzVector
import glob
import subprocess
from subprocess import call
import argparse

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptFit(0);
ROOT.gStyle.SetOptStat(0);
ROOT.gStyle.SetPadBottomMargin(0.15)
ROOT.gStyle.SetPadLeftMargin(0.13)
# storage =  ROOT.gSystem.ExpandPathName("$STORAGEDIR")
storage =  os.path.expandvars("$STORAGEDIR")
#storage = "/Volumes/Study/Weizmann_PostDoc/TomPtarmiganFiles/"

meMeV = 0.5109989461 ## MeV
meGeV = meMeV/1000.
MeV2GeV = 1./1000.

def beta2(betax,betay,betaz):
   return betax*betax+betay*betay+betaz*betaz

def gamma(betax,betay,betaz):
   return 1./math.sqrt(1.-beta2(betax,betay,betaz))

def elmomentum(bx,by,bz,E0):
   b2 = beta2(bx,by,bz)
   gm = gamma(bx,by,bz)
   E = gm*meGeV
   px = bx*E
   py = by*E
   pz = bz*E
   if((E0-E)/E0>0.00001): print("E=%g, E0=%g" % (E,E0))
   return px,py,pz,E

def gmmomentum(bx,by,bz,E):
   px = bx*E
   py = by*E
   pz = bz*E
   return px,py,pz,E

def readparticles(name,xivalue=7.0,mpids=[]):
   print("reading: ",name)
   filtermpids = (len(mpids)>0)
   particles = {}
   i = 0
   fIn = h5py.File(name, 'r')
   
   ### electrons
   id_value_electron       = fIn['final-state/electron']['id'][()]
   parentid_value_electron = fIn['final-state/electron']['parent_id'][()]
   momentum_value_electron = fIn['final-state/electron']['momentum'][()]
   position_value_electron = fIn['final-state/electron']['position'][()]
   weight_value_electron   = fIn['final-state/electron']['weight'][()]
   
   for j in range(0, len(id_value_electron)):
        vx    = position_value_electron[j][0]
        vy    = position_value_electron[j][1]
        vz    = position_value_electron[j][2]
        time  = position_value_electron[j][3]
        E     = momentum_value_electron[j][0]
        px    = momentum_value_electron[j][1]
        py    = momentum_value_electron[j][2]
        pz    = momentum_value_electron[j][3]
        pdgId = 11
        wgt   = weight_value_electron[j]
        mpid  = str(id_value_electron[j])+"_"+str(pdgId)
        xi    = xivalue
        ### filter out mpids which are not in the stdhep array
        if(filtermpids and mpid not in mpids): continue
    
        particles.update( {mpid:{"E":E, "vx":vx, "vy":vy, "vz":vz, "px":px, "py":py, "pz":pz, "pdgId":pdgId, "wgt":wgt, "time":time, "xi":xi}} )

        
   ### positrons
   id_value_positron       = fIn['final-state/positron']['id'][()]
   parentid_value_positron = fIn['final-state/positron']['parent_id'][()]
   momentum_value_positron = fIn['final-state/positron']['momentum'][()]
   position_value_positron = fIn['final-state/positron']['position'][()]
   weight_value_positron   = fIn['final-state/positron']['weight'][()]
    
   for j in range(0, len(id_value_positron)):
        vx    = position_value_positron[j][0]
        vy    = position_value_positron[j][1]
        vz    = position_value_positron[j][2]
        time  = position_value_positron[j][3]
        E     = momentum_value_positron[j][0]
        px    = momentum_value_positron[j][1]
        py    = momentum_value_positron[j][2]
        pz    = momentum_value_positron[j][3]
        pdgId = -11
        wgt   = weight_value_positron[j]
        mpid  = str(id_value_positron[j])+"_"+str(pdgId)
        xi    = xivalue ### hardcoded now
        ### filter out mpids which are not in the stdhep array
        if(filtermpids and mpid not in mpids): continue
    
        particles.update( {mpid:{"E":E, "vx":vx, "vy":vy, "vz":vz, "px":px, "py":py, "pz":pz, "pdgId":pdgId, "wgt":wgt, "time":time, "xi":xi}} )
        
   '''
   ### photons
   id_value_photon       = fIn['final-state/photon']['id'][()]
   parentid_value_photon = fIn['final-state/photon']['parent_id'][()]
   momentum_value_photon = fIn['final-state/photon']['momentum'][()]
   position_value_photon = fIn['final-state/photon']['position'][()]
   weight_value_photon   = fIn['final-state/photon']['weight'][()]
   xi_value_photon       = fIn['final-state/photon']['xi'][()]
    
   for j in range(0, len(id_value_photon)):
        E     = momentum_value_photon[j][0]
        vx    = position_value_photon[j][0]
        vy    = position_value_photon[j][1]
        vz    = position_value_photon[j][2]
        px    = momentum_value_photon[j][1]
        py    = momentum_value_photon[j][2]
        pz    = momentum_value_photon[j][3]
        pdgId = 22
        wgt   = weight_value_photon[j]
        mpid  = str(id_value_photon[j])+"_"+str(pdgId)
        time  = position_value_photon[j][3]
        xi    = xi_value_photon[j]
        ### filter out mpids which are not in the stdhep array
        if(filtermpids and mpid not in mpids): continue
    
        particles.update( {mpid:{"E":E, "vx":vx, "vy":vy, "vz":vz, "px":px, "py":py, "pz":pz, "pdgId":pdgId, "wgt":wgt, "time":time, "xi":xi}} )
   '''
   fIn.close()
   print("got %g particles" % len(particles))
   return particles


def main():
    
    parser = argparse.ArgumentParser(description='Code to transform h5 format to root')
    parser.add_argument('-p', action="store", required=True, dest="processName", type=str, default="glaser")
    parser.add_argument('-t', action="store", required=True, dest="phase", type=str, default="phase0")
    parser.add_argument('-x', action="store", required=True, dest="xi", type=str, default="3.0")
    parser.add_argument('-g', action="store_true", required=False, dest="needGLaser")
    argus = parser.parse_args()

     
    phase   = argus.phase
    xiStr   = argus.xi 
    xiInput = float(argus.xi)
    photon  = argus.needGLaser
    process = argus.processName
    
    if process=="glaser":
        indir = "brem-laser"
    else:
        indir = "e-laser"
        
    path    = storage+"/data/h5/"+process+"/"+phase+"/"+xiStr+"/"
    
    #### replace this while running on DESY
    targetdir  = storage+"/data/root/raw/"+process+"/"+phase+"/"+xiStr+"/"
    p         = subprocess.Popen("mkdir -p "+targetdir, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err  = p.communicate()


    tf = TFile( targetdir+'/raw_'+process+'.root', 'recreate' )

    tt_out    = TTree( 'tt','tt' )
    vx_out    = ROOT.std.vector( float )()
    vy_out    = ROOT.std.vector( float )()
    vz_out    = ROOT.std.vector( float )()
    px_out    = ROOT.std.vector( float )()
    py_out    = ROOT.std.vector( float )()
    pz_out    = ROOT.std.vector( float )()
    E_out     = ROOT.std.vector( float )()
    pdgId_out = ROOT.std.vector( int )()
    mpid_out  = ROOT.std.vector( str )()
    wgt_out   = ROOT.std.vector( float )()
    time_out  = ROOT.std.vector( float )()
    xi_out    = ROOT.std.vector( float )()
    tt_out.Branch('vx', vx_out)
    tt_out.Branch('vy', vy_out)
    tt_out.Branch('vz', vz_out)
    tt_out.Branch('px', px_out)
    tt_out.Branch('py', py_out)
    tt_out.Branch('pz', pz_out)
    tt_out.Branch('E',  E_out)
    tt_out.Branch('wgt',  wgt_out)
    tt_out.Branch('pdgId',pdgId_out)
    tt_out.Branch('mpid',mpid_out)
    tt_out.Branch('time',time_out)
    tt_out.Branch('xi',xi_out)

    print("path=",path)
    print("targetdir=",targetdir)
    fIns = glob.glob(path+"/*.h5")
    print(fIns)   
 
    counter = 0
    for name in fIns:
        ### clear output tree branches
        mpid_out.clear()
        pdgId_out.clear()
        wgt_out.clear()
        vx_out.clear()
        vy_out.clear()
        vz_out.clear()
        px_out.clear()
        py_out.clear()
        pz_out.clear()
        E_out.clear()
        xi_out.clear()
        time_out.clear()

        ### read files
        particles    = readparticles(name, xiInput)
        
        ### loop over all output particles
        for MP_ID,particle in particles.items():
            if(counter%10000==0): print("Processed: ", counter)
            px0    = particle["px"]
            py0    = particle["py"]
            pz0    = particle["pz"]
            Energy = particle["E"]
            wgt0   = particle["wgt"]
            pdgId0 = particle["pdgId"]
            vx0    = particle["vx"]*1.e-1 ## mm to cm
            vy0    = particle["vy"]*1.e-1 ## mm to cm
            vz0    = particle["vz"]*1.e-1 ## mm to cm
            xi0    = particle["xi"]
            t0     = particle["time"]
            mpid_out.push_back(str(MP_ID))
            wgt_out.push_back(wgt0)  
            pdgId_out.push_back(int(pdgId0))  
            vx_out.push_back(vx0)
            vy_out.push_back(vy0)
            vz_out.push_back(vz0)
            px_out.push_back(px0)
            py_out.push_back(py0)
            pz_out.push_back(pz0)
            E_out.push_back(Energy)
            time_out.push_back(t0)
            xi_out.push_back(xi0)
            counter += 1
            
        tt_out.Fill()
        
    
    tt_out.Write()
    tf.Write()
    tf.Write()
    tf.Close()


if __name__=="__main__":
    intime = time.time()
    main()
    print("----- the time taken ", time.time() - intime, " s")
