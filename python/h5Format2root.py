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


storage =  os.path.expandvars("$STORAGEDIR")

meMeV = 0.5109989461 ## MeV
meGeV = meMeV/1000.
MeV2GeV = 1./1000.

def main():
    
    parser = argparse.ArgumentParser(description='Code to transform h5 format to root')
    parser.add_argument('-process', action="store", dest="processName", type=str, default="glaser")
    parser.add_argument('-phase', action="store", dest="phase", type=str, default="phase0")
    parser.add_argument('-xi', action="store", dest="xi", type=str, default="3.0")
    parser.add_argument('-g', action="store_true", dest="needGLaser")
    parser.add_argument('-polarization', action="store", dest="polarizationName", default="gpc")
    
    argus = parser.parse_args()
     
    phase        = argus.phase
    xiStr        = argus.xi 
    xiInput      = float(argus.xi)
    photon       = argus.needGLaser
    process      = argus.processName
    polarization = argus.polarizationName
    
    if process=="glaser":
        indir = "brem-laser"
    else:
        indir = "e-laser"
        
    
    
    ### for running on DESY
    #path      =  "/nfs/dust/luxe/group/MCProduction/Signal/ptarmigan-v0.8.1/"+indir+"/"+phase+"/"+polarization+"/"+xiStr
    #storage   = "TomPtarmiganFiles"
    
    
    #### for running locally
    path      = storage+"/data/h5/"+process+"/"+phase+"/"+polarization+"/"+xiStr+"/"
    
    targetdir = storage+"/data/root/raw/"+process+"/"+phase+"/"+polarization+"/"+xiStr+"/"
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
 
    counter  = 0
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
        
        #### input file
        fIn = h5py.File(name, 'r')
        print("reading: ",name)
        
        ### electrons
        id_value_electron       = fIn['final-state/electron']['id'][()]
        parentid_value_electron = fIn['final-state/electron']['parent_id'][()]
        momentum_value_electron = fIn['final-state/electron']['momentum'][()]
        position_value_electron = fIn['final-state/electron']['position'][()]
        weight_value_electron   = fIn['final-state/electron']['weight'][()]
        
        for j in range(0, len(id_value_electron)):
            if(counter%1000==0): print("Processed: ", counter)
            vx0    = position_value_electron[j][0]*1.e-1 ## mm to cm
            vy0    = position_value_electron[j][1]*1.e-1 ## mm to cm
            vz0    = position_value_electron[j][2]*1.e-1 ## mm to cm
            t0     = position_value_electron[j][3]
            Energy = momentum_value_electron[j][0]
            px0    = momentum_value_electron[j][1]
            py0    = momentum_value_electron[j][2]
            pz0    = momentum_value_electron[j][3]
            pdgId0 = 11
            wgt0   = weight_value_electron[j]
            MP_ID  = str(id_value_electron[j])+"_"+str(pdgId0)
            xi0    = xiInput
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
                

                
        ### positrons
        id_value_positron       = fIn['final-state/positron']['id'][()]
        parentid_value_positron = fIn['final-state/positron']['parent_id'][()]
        momentum_value_positron = fIn['final-state/positron']['momentum'][()]
        position_value_positron = fIn['final-state/positron']['position'][()]
        weight_value_positron   = fIn['final-state/positron']['weight'][()]
            
        for j in range(0, len(id_value_positron)):
            if(counter%1000==0): print("Processed: ", counter)
            vx0    = position_value_positron[j][0]*1.e-1 ## mm to cm
            vy0    = position_value_positron[j][1]*1.e-1 ## mm to cm
            vz0    = position_value_positron[j][2]*1.e-1 ## mm to cm
            t0     = position_value_positron[j][3]
            Energy = momentum_value_positron[j][0]
            px0    = momentum_value_positron[j][1]
            py0    = momentum_value_positron[j][2]
            pz0    = momentum_value_positron[j][3]
            pdgId0 = -11
            wgt0   = weight_value_positron[j]
            MP_ID  = str(id_value_positron[j])+"_"+str(pdgId0)
            xi0    = xiInput
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
