//#include <TClonesArray.h> 

#include "iostream"

#include <TDatabasePDG.h>
#include <TFile.h>
#include "GenMUONLMR.h" 
#include "TRandom.h"
#include "TParticle.h"

ClassImp(GenMUONLMR)
  
  
GenMUONLMR::GenMUONLMR(Double_t energy, Int_t kLowEnergy) :  
                                                   fNMuMin(2), 
						   fGenSingleProc(-1),
						   fCosTheta(0x0), 
						   fRhoLineShape(0x0),  
						   fHMultMu(0x0), 
						   fHNProc(0x0) { 
    //
    // default constructor 
    //
    // initialize pt and y distributions according to a fit to 
    // Pythia simulation at sqrt(s) = 7 TeV
    for (Int_t ipart=0; ipart < fgkNpart; ipart++) fScaleMult[ipart] = 1; 
    fScaleMult[kPionLMR] = 0; // set pion multiplicity to zero 
    fScaleMult[kKaonLMR] = 0; // set kaon multiplicity to zero
    Int_t pdg[9] = {211, 321, 221, 113, 223, 333, 331, 443, 443}; 
    const char* fptname[9] = {"fPtPion","fPtKaon","fPtEta","fPtRho","fPtOmega","fPtPhi","fPtEtaPrime","fPtJPsi","fPtJPsi"};
    const char* fyname[9] = {"fYPion","fYKaon","fYEta","fYRho","fYOmega","fYPhi","fYEtaPrime","fYJPsi","fYJpsi_Sch"}; 
    const char* fnname[9] = {"fMultPion","fMultKaon","fMultEta","fMultRho","fMultOmega","fMultPhi","fMultEtaPrime","fMultJPsi","fMultJPsi"};
    const char* fdname[2] = {"fDecPion","fDecKaon"};
    Double_t ctau[2] = {7.8045, 3.712};  
    Double_t ptparam[9][9];
    Double_t yparam[9][9];
    Double_t nparam[9][9];
    
    // parameters for 7 TeV generation
    if (energy==7.0) {
      printf ("GenMUONLMR: using pp parameterization at 7 TeV\n");  
      Double_t ptparam7000[9][9] = {{1,0.427,2.52,0,0,0,0,0,0}, // pions from Pythia
				    {1,0.58,2.57,0,0,0,0,0,0},  // kaons from Pythia
				    {1,0.641,2.62,0,0,0,0,0,0}, // eta from Pythia
				    {1,1.44,3.16,0,0,0,0,0,0},  // rho+omega from ALICE muon  
				    {1,1.44,3.16,0,0,0,0,0,0},  // rho+omega from ALICE muon  
				    {1,1.16,2.74,0,0,0,0,0,0},  // phi from ALICE muon  
				    {1,0.72,2.5,0,0,0,0,0,0},  // etaPrime from Pythia    
      				    {1,0.72,2.5,0,0,0,0,0,0},   // J/psi - meaningless     
      				    {1,0.72,2.5,0,0,0,0,0,0}};  // J/psi - meaningless     

      Double_t yparam7000[9][9] = {{1,0.8251,3.657,0,0,0,0,0,0}, // pions from pythia
				   {1,1.83,2.698,0,0,0,0,0,0},   // kaons from pythia
				   {1,1.169,3.282,0,0,0,0,0,0},  // eta from pythia
				   {1,1.234,3.264,0,0,0,0,0,0},  // rho from pythia
				   {1,1.311,3.223,0,0,0,0,0,0},  // omega from pythia
				   {1,2.388,2.129,0,0,0,0,0,0},  // phi from pythia
				   {1,1.13,3.3,0,0,0,0,0,0},    // eta prime from pythia
      				   {1,1.13,3.3,0,0,0,0,0,0},     // J/psi - meaningless
      				   {1,1.13,3.3,0,0,0,0,0,0}};    // J/psi - meaningless

      // multiplicity parameters from pythia
      Double_t nparam7000[9][9] = {{353.582, 6.76263, 1.66979, 998.445, 9.73281, 12.6704, 175.187, 29.08, 40.2531},
				   {1.e4,  0.2841, 0,0,0,0,0,0,0},
				   {1.e4,  0.2647, 0,0,0,0,0,0,0},
				   {7055,  0.1786, 0,0,0,0,0,0,0},
				   {7500,  0.1896, 0,0,0,0,0,0,0},
				   {5.e4,  1.167,  0,0,0,0,0,0,0}, 
				   {2.9e4, 0.714,  0,0,0,0,0,0,0},
				   {2.9e4, 0.714,  0,0,0,0,0,0,0}, // J/psi - meaningless
				   {2.9e4, 0.714,  0,0,0,0,0,0,0}}; // J/psi - meaningless

      for (Int_t i=0; i<fgkNpart; i++) { 
	for (Int_t j=0; j<9; j++) {
	  ptparam[i][j] = ptparam7000[i][j];
	  yparam[i][j] = yparam7000[i][j];
	  nparam[i][j] = nparam7000[i][j];	
	}
      }
    }  
    
    // parameters for 2.76 generation
    // pt params has been determined as <pt>ALICE_2.76 = <pt>ALICE_7 * <pt>PYTHIA_2.76 / <pt>PYTHIA_7
    if (energy == 2.76){
      printf ("GenMUONLMR: using pp parameterization at 2.76 TeV\n");  
      Double_t ptparam2760[9][9] = {{1,0.1665,8.878,0,0,0,0,0,0},  // pions from Pythia
				    {1,0.1657,8.591,0,0,0,0,0,0},  // kaons from Pythia
				    {1,0.641,2.62,0,0,0,0,0,0},    // eta from ALICE 7 TeV
				    {1,1.3551,3.16,0,0,0,0,0,0},   // rho with <pt> scaled 
				    {1,1.3551,3.16,0,0,0,0,0,0},   // omega with <pt> scaled 
				    {1,1.0811,2.74,0,0,0,0,0,0},   // phi with <pt> scaled 
				    {1,0.72,2.5,0,0,0,0,0,0},     // etaPrime from ALICE 7 TeV
      				    {1,0.72,2.5,0,0,0,0,0,0},      // J/psi - meaningless
      				    {1,0.72,2.5,0,0,0,0,0,0}};     // J/psi - meaningless

      Double_t yparam2760[9][9] = {{1,0.8251,3.657,0,0,0,0,0,0},      // pions from pythia
				   {1,1.83,2.698,0,0,0,0,0,0},        // kaons from pythia
				   {1,0.011,3.474,0,0,0,0,0,0},       // eta from pythia
				   {1,-0.01,3.409,0,0,0,0,0,0},       // rho from pythia
				   {1,-0.037,3.294,0,0,0,0,0,0},      // omega from pythia
				   {1,-0.016,2.717,0,0,0,0,0,0},      // phi from pythia
				   {1,-0.010,3.312,0,0,0,0,0,0},     // eta prime from pythia  
				   {1,-0.010,3.312,0,0,0,0,0,0},      // J/psi - meaningless 
				   {1,-0.010,3.312,0,0,0,0,0,0}};     // J/psi - meaningless 

      
      Double_t nparam2760[9][9] = {{1.e6,535.,0,0,0,0,0,0,0},                                 // pions
				   {1.e5,70,0,0,0,0,0,0,0},                                   // kaons
				   {1.e4,0.351,0,0,0,0,0,0,0},                                // eta
				   {1.e4,0.2471,0,0,0,0,0,0,0},                               // rho
				   {1.e4,0.2583,0,0,0,0,0,0,0},                               // omega
				   {1.e5,1.393,0,0,0,0,0,0,0},                                // phi
				   {1.e4,0.9005,0,0,0,0,0,0,0},                              // etaPrime
				   {1.e4,0.9005,0,0,0,0,0,0,0},                              // J/psi meaningless
				   {1.e4,0.9005,0,0,0,0,0,0,0}};                              // J/psi meaningless

      
      for (Int_t i=0; i<fgkNpart; i++) { 
	for (Int_t j=0; j<9; j++) {
	  ptparam[i][j] = ptparam2760[i][j];
	  yparam[i][j] = yparam2760[i][j];
	  nparam[i][j] = nparam2760[i][j];	
	}
      }
    }
    
    TDatabasePDG* pdgp = TDatabasePDG::Instance();

    for (Int_t i=0; i<fgkNpart; i++) { 
      printf("\n"); //problem with new aliroot version
      fPDG[i] = pdg[i]; 
      if (i>1) { // if (i!=0) se non includi i kaoni
	fMult[i] = new TF1(fnname[i],"[0]*exp(-[1]*x)",0,30);
	fMult[i]->SetParameters(nparam[i][0],nparam[i][1]);  
      }
      else { 
	//	fMult[i] = new TF1(fnname[i],"gaus(0)+gaus(3)+gaus(6)",0,150);
	//	for (Int_t j=0; j<9; j++) fMult[i]->SetParameter(j,nparam[i][j]);
	fMult[i] = new TF1("fnname[i]","[0]*TMath::Poisson(x,[1])",0,30);
	fMult[i]->SetParameters(nparam[i][0],nparam[i][1]);  
      }
      
      if(kLowEnergy == 0){
	fPt[i] = new TF1(fptname[i],GenMUONLMR::PtDistr,0,20,3);
	fPt[i]->SetParameters(ptparam[i][0], ptparam[i][1], ptparam[i][2]);  
	fY[i] = new TF1(fyname[i],GenMUONLMR::YDistr,-10,10,3);
	fY[i]->SetParameters(yparam[i][0], yparam[i][1], yparam[i][2]); 
      }
      else {
        if(i != 7 && i != 8){
	fPt[i] = new TF1(fptname[i],GenMUONLMR::PtDistrLowEnergy,0,20,3);
	fPt[i]->SetParameters(ptparam[i][0], ptparam[i][1], pdgp->GetParticle(fPDG[i])->Mass());  
	printf("mass=%f\n",pdgp->GetParticle(fPDG[i])->Mass());
	fPt[i]->FixParameter(2,pdgp->GetParticle(fPDG[i])->Mass());
	} else {
	fPt[i] = new TF1(fptname[i],GenMUONLMR::PtDistrPythia6,0,10,4);
	fPt[i]->SetParameters(ptparam[i][0], ptparam[i][1], ptparam[i][2], ptparam[i][3]); 
      }
	if(i != 8){
      	fY[i] = new TF1(fyname[i],GenMUONLMR::YDistrLowEnergy,-10.,10.,3);
	fY[i]->SetParameters(yparam[i][0], yparam[i][1], yparam[i][2]);
	} else {
	fY[i] = new TF1(fyname[i],GenMUONLMR::YDistrSch,-10.,10.,5);
	fY[i]->SetParameters(yparam[i][0], yparam[i][1], yparam[i][2], yparam[i][3], yparam[i][4]);

	}
      }
    }
   
    for(Int_t i = 0; i<2; i++){
      fDecay[i] = new TF1(fdname[i],"exp(-x/[0])",0,150);
      fDecay[i]->SetParameter(0,ctau[i]);
    }
    
    for (Int_t ipart = 0; ipart < fgkNpart; ipart++) { 
      fParticle[ipart] = new TParticle(); 
      fParticle[ipart]->SetPdgCode(fPDG[ipart]); 
    }
    
    TDatabasePDG *pdgdb = TDatabasePDG::Instance(); 
    Double_t mumass = pdgdb->GetParticle(13)->Mass();
    fMu[0] = new TParticle(); 
    fMu[0]->SetPdgCode(-13); 
    fMu[0]->SetCalcMass(mumass); 
    fMu[1] = new TParticle(); 
    fMu[1]->SetPdgCode(13); 
    fMu[1]->SetCalcMass(mumass); 
    
    // function for polarized theta distributions
    // fCosTheta = new TF1 ("fCosTheta","1+[0]*x*x",-1,1);
//     fCosTheta->SetParameter(0,1);




    
    // Dalitz decays 
    Int_t nbins = 1000;
    Double_t xmin = 0, xmax = 2; 
    fDalitz[0] = new TH1F("hDalitzEta","",nbins,xmin,xmax);
    fDalitz[1] = new TH1F("hDalitzOmega","",nbins,xmin,xmax);
    fDalitz[2] = new TH1F("hDalitzEtaPrime","",nbins,xmin,xmax);
    
    Double_t meta   = pdgdb->GetParticle("eta")->Mass(); 
    Double_t momega = pdgdb->GetParticle("omega")->Mass(); 
    Double_t metaPrime = pdgdb->GetParticle("eta'")->Mass(); 
    Double_t mpi0   = pdgdb->GetParticle("pi0")->Mass(); 
    Double_t md3 = 0, mres = 0; 
    
    for (Int_t index = 0; index < 3; index++) { 
      if (index == 0) { 
	mres = meta; 
	md3 = 0; 
      }
      else if (index == 1) { 
	mres = momega; 
	md3 = mpi0; 
      }
      else if (index == 2) { 
	mres = metaPrime; 
	md3 = 0; 
      }
      Double_t delta   = md3 * md3 / (mres * mres);
      Double_t epsilon = mumass * mumass / (mres * mres);
      nbins = fDalitz[index]->GetNbinsX();
      xmin = fDalitz[index]->GetXaxis()->GetXmin(); 
      Double_t deltax =  fDalitz[index]->GetBinWidth(1);
      Double_t xd = xmin - deltax/2.; 
      for (Int_t ibin = 0; ibin< nbins; ibin++) { 
	Double_t dalval = 0; 
	xd += deltax; 
	if (xd > 4. *epsilon) { 
	  Double_t bracket = TMath::Power(1. + xd/(1. - delta),2)      
	    - 4. * xd / ((1. - delta) * (1. - delta));
	  if (bracket > 0) { 
	    dalval = TMath::Power(bracket,1.5) /xd *
	      TMath::Sqrt(1 - 4 * epsilon / xd) * (1 + 2 * epsilon / xd) * 
	      FormFactor(xd * mres * mres, index);
	    fDalitz[index]->Fill(xd,dalval); 
	  }
	}
      }
    }
    
    fRhoLineShape = new TF1("fRhoLineShape",RhoLineShapeNew,0,2,2); 
    fHMultMu = new TH1D("fHMultMu","Muon multiplicity",20,-0.5,19.5); 
    fHNProc = new TH1D("fHNProc","Number of gen. evts. per process in 4 pi",9,-0.5,8.5); 
  }

//-----------------------------------------------------------

GenMUONLMR::GenMUONLMR(GenMUONLMR &gen) :  
						    fNMuMin(gen.fNMuMin), 
						    fGenSingleProc(gen.fGenSingleProc),
						    fCosTheta(gen.fCosTheta), 
						    fRhoLineShape(gen.fRhoLineShape),  
						    fHMultMu(gen.fHMultMu), 
						    fHNProc(gen.fHNProc) {  
  for (Int_t i=0; i < fgkNpart; i++) { 
    fPDG[i] = gen.fPDG[i]; 
    fScaleMult[i] = gen.fScaleMult[i]; 
    fPt[i] = (TF1*) gen.fPt[i]->Clone(); 
    fY[i] = (TF1*) gen.fY[i]->Clone();  
    fMult[i] = (TF1*) gen.fMult[i]->Clone(); 
    fParticle[i] = (TParticle*) gen.fParticle[i]->Clone(); 
  }
  
  for(Int_t i = 0; i<2; i++) fDecay[i] = (TF1*) gen.fDecay[i]->Clone(); 
  for(Int_t i = 0; i<3; i++) fDalitz[i] = (TH1F*) gen.fDalitz[i]->Clone(); 
  for(Int_t i = 0; i<2; i++) fMu[i] = (TParticle*) gen.fMu[i]->Clone(); 
}

//-----------------------------------------------------------

GenMUONLMR& GenMUONLMR::operator=(const GenMUONLMR &gen) {
  fNMuMin = gen.fNMuMin; 
  fGenSingleProc = gen.fGenSingleProc; 
  fCosTheta = (TF1*) gen.fCosTheta->Clone();  
  fRhoLineShape = (TF1*) gen.fRhoLineShape->Clone();
  fHMultMu = (TH1D*) gen.fHMultMu->Clone();
  fHNProc = (TH1D*) gen.fHNProc->Clone();  
  
  for (Int_t i=0; i < fgkNpart; i++) { 
    fPDG[i] = gen.fPDG[i]; 
    fScaleMult[i] = gen.fScaleMult[i]; 
    fPt[i] = (TF1*) gen.fPt[i]->Clone(); 
    fY[i] = (TF1*) gen.fY[i]->Clone();  
    fMult[i] = (TF1*) gen.fMult[i]->Clone(); 
    fParticle[i] = (TParticle*) gen.fParticle[i]->Clone(); 
  }
  
  for(Int_t i = 0; i<2; i++) fDecay[i] = (TF1*) gen.fDecay[i]->Clone(); 
  for(Int_t i = 0; i<3; i++) fDalitz[i] = (TH1F*) gen.fDalitz[i]->Clone(); 
  for(Int_t i = 0; i<2; i++) fMu[i] = (TParticle*) gen.fMu[i]->Clone(); 
  return *this; 
}

 
//-----------------------------------------------------------

GenMUONLMR::~GenMUONLMR()
{
  // Default destructor
  for (Int_t i=0; i<8; i++) { 
    delete fPt[i]; 
    delete fY[i]; 
    delete fMult[i]; 
    delete fParticle[i]; 
  }    
  
  for (Int_t i=0; i<2; i++) { 
    delete fDecay[i]; 
    delete fMu[i]; 
  }

  for (Int_t i=0; i<3; i++) delete fDalitz[i]; 

  delete fCosTheta; fCosTheta = 0;  
  delete fRhoLineShape; fRhoLineShape = 0;  
  delete fHMultMu; fHMultMu = 0;   
  delete fHNProc;  fHNProc = 0;   
}

//-----------------------------------------------------------

void GenMUONLMR::FinishRun(){ 
  // save some histograms to an output file 
  Int_t nbins = fHNProc->GetNbinsX(); 
  for (Int_t ibin=1; ibin <= nbins; ibin++) printf ("ibin = %d nEvProc = %g\n",
						      ibin,fHNProc->GetBinContent(ibin));
    TFile *fout = new TFile("GenMUONLMR_histos.root","recreate"); 
  fHMultMu->Write(); 
  fHNProc->Write(); 
  fout->Close(); 
}

//-----------------------------------------------------------

Double_t GenMUONLMR::YDistr(Double_t *px, Double_t *par){ 
  // function for rapidity distribution: plateau at par[0] +
  // gaussian tails centered at par[1] and with par[2]=sigma  
  Double_t func = 0.;
  Double_t x = TMath::Abs(px[0]);
  if (x<par[1]) func = par[0]; 
  else { 
    Double_t z = (x-par[1])/(par[2]); 
    func = par[0] * TMath::Exp(-0.5 * z * z); 
  }
  
  
  return func; 
}

Double_t GenMUONLMR::YDistrLowEnergy(Double_t *px, Double_t *par){ 
  // gaussian centered at par[1] and with par[2]=sigma  
  Double_t func = 0.;
  
  Double_t x = px[0];
  Double_t z = (x-par[1])/(par[2]); 
  func = par[0] * TMath::Exp(-0.5 * z * z); 
  

  return func; 
}

Double_t GenMUONLMR::YDistrSch(Double_t *px, Double_t *par){ 
  // param. obtained transforming dsigma/dy = (1-TMath::Abs(x))**d (Schuler, Tesi Fleuret)
  // following formula gives ycms. 
  Double_t func = 0.;
  
  Double_t ycms = par[4];
  Double_t ymax = TMath::ASinH(par[2]/(2*par[1]))+ycms;
  Double_t ymin = TMath::ASinH(-par[2]/(2*par[1]))+ycms;
  
  Double_t x = px[0];  
  if(x>ymin && x<ymax){
    func= par[0]*TMath::Power((1-2*par[1]/par[2]*TMath::Abs(sinh(x-ycms))),par[3])*2*par[1]/par[2]*cosh(x-ycms);
  } else func = 0;
  
  return func; 
}

//-----------------------------------------------------------

Double_t GenMUONLMR::PtDistr(Double_t *px, Double_t *par){
  // pt distribution: power law 
  Double_t x = px[0];

  Double_t func = 0.;
  func = par[0] * x / TMath::Power((1+(x/par[1])*(x/par[1])),par[2]); 

  return func; 
}

Double_t GenMUONLMR::PtDistrLowEnergy(Double_t *px, Double_t *par){
  // pt distribution: mt scaling law 
  Double_t x = px[0];

  Double_t func = 0.;
  func = par[0] * x * TMath::Exp(-pow((x*x + par[2] * par[2]),0.5)/par[1]);

  return func; 
}

Double_t GenMUONLMR::PtDistrPythia6(Double_t *px, Double_t *par){
  Double_t x = px[0];

  Double_t func = 0.;
  func = par[0]*x/TMath::Power((1+TMath::Power(x/par[1],par[3])),par[2]);
  return func; 
}

//-----------------------------------------------------------

Double_t GenMUONLMR::CosThetaParam(Double_t *px, Double_t *par) {
  //costheta distribution for dalitz decay
  
  TDatabasePDG *pdgdb = TDatabasePDG::Instance(); 
  Double_t mumass = pdgdb->GetParticle(13)->Mass();
  Double_t cosTheta2 = px[0]*px[0];
  Double_t mass = par[0];
  
  return (1+cosTheta2)+ TMath::Power(2*mumass/mass,2)*(1-cosTheta2);
  
}
//-----------------------------------------------------------

void GenMUONLMR::Generate() {
  //
  // generate the low mass resonances and their decays according to  
  // the multiplicity parameterized by pythia and BR from PDG  
  // rapidity distributions parametrized from pythia 
  // pt distributions from data (or pythia for etaprime) 
  //
  Double_t pxPushed[100], pyPushed[100], pzPushed[100], ePushed[100]; 
  Int_t nmuons = -1, npartPushed = 0, pdgPushed[100]; 
  Double_t polar[3]= {0,0,0};  // Polarisation of the parent particle (for GEANT tracking)
//   Double_t origin0[3];         // Origin of the generated parent particle (for GEANT tracking)
//   // Calculating vertex position per event
//   for (Int_t j=0;j<3;j++) origin0[j]=fOrigin[j];
//   if(fVertexSmear==kPerEvent) {
//     Vertex();
//     for (Int_t j=0;j<3;j++) origin0[j]=fVertex[j];
//   }
  
  TParticle *mother; 
  TDatabasePDG* pdg = TDatabasePDG::Instance();

  Double_t pt, y, phi, mass, px, py, pz, ene, mt; 

  const Int_t nproc = 11; 
  Int_t idRes[nproc] = {kEtaLMR, kEtaLMR, kRhoLMR, kOmegaLMR, kOmegaLMR, kPhiLMR, kEtaPrimeLMR, kPionLMR, kKaonLMR,kJPsi,kJPsiSch}; 
  Double_t BR[nproc] = {5.8e-6, 3.1e-4, 4.55e-5, 7.28e-5, 1.3e-4, 2.86e-4, 1.04e-4, 1, 0.6344, 0.05, 0.05};
  //  Double_t BR[nproc] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  Int_t idDec[nproc] = {0, 1, 0, 0, 1, 0, 1, 2, 2, 0, 0}; // 0:2body, 1:Dalitz, 2:pi/K 
  Int_t mult[nproc] = {0,0,0,0,0,0,0,0,0,0,0}; 
  //printf("fNMuMin=%d\n",fNMuMin);
  //  while (nmuons < fNMuMin) { 

    nmuons = 0; 
    npartPushed = 0; 
    for (Int_t iproc=0; iproc<nproc; iproc++) { 
      if (fGenSingleProc == -1) { 
	mult[iproc] = Int_t(fMult[idRes[iproc]]->GetRandom()*fScaleMult[idRes[iproc]]); 
      }
      else { 
	if (iproc==fGenSingleProc) { 
	  mult[iproc] = 1; 
	  BR[iproc] = 1;
	} 
	else { 
	  mult[iproc] = 0; 
	  BR[iproc] = 0;
	}
      }
    }
    
    if (fGenSingleProc == -1) { 
      mult[1] = mult[0]; 
      mult[4] = mult[3]; 
    }
    //printf("nproc=%d\n",nproc);    
    for (Int_t iproc = 0; iproc < nproc; iproc++) { 
      //printf("mult[%d]=%d BR[%d]=%f\n",iproc,mult[iproc],iproc,BR[iproc]); 
      for (Int_t imult=0; imult<mult[iproc]; imult++) { 
	if (gRandom->Rndm() < BR[iproc]) { 
	  fHNProc->Fill(iproc); 
	  Int_t ipart = idRes[iproc]; 
	  pt  = fPt[ipart]->GetRandom(); 
	  y   = fY[ipart]->GetRandom(); 
	  phi = gRandom->Rndm() * 2 * TMath::Pi(); 
	  mass = pdg->GetParticle(fPDG[ipart])->Mass(); 
	  px  = pt * TMath::Cos(phi); 
	  py  = pt * TMath::Sin(phi); 
	  mt  = TMath::Sqrt(pt * pt + mass * mass);
	  pz  = mt * TMath::SinH(y); 
	  ene = mt * TMath::CosH(y); 
	
	  mother = fParticle[ipart]; 
	  mother->SetMomentum(px,py,pz,ene); 
	  mother->SetCalcMass(mass);
	  //	  fParticle[ipart]->SetMomentum(px,py,pz,ene);  // gu
	  //	  if (!KinematicSelection(mother,0)) continue; 
	  
	  Bool_t hasDecayed = kTRUE;
	  if (idDec[iproc] == 0) Decay2Body(mother);
	  else if (idDec[iproc] == 1) DalitzDecay(mother); 
	  else DecayPiK(mother,hasDecayed); 
	  if (!hasDecayed) continue; 
// 	  Bool_t isMu0Acc = KinematicSelection(fMu[0],1); 
// 	  Bool_t isMu1Acc = KinematicSelection(fMu[1],1); 
// 	  Bool_t isMuFromPiKAcc = kTRUE; 

// 	  if (idDec[iproc] == 2) isMuFromPiKAcc = (mother->GetPdgCode()>0) ? isMu0Acc : isMu1Acc;
// 	  // mother 
// 	  if ((idDec[iproc]  < 2 && (isMu0Acc || isMu1Acc)) || 
// 	      (idDec[iproc] == 2 && isMuFromPiKAcc)) { 
// 	    pdgPushed[npartPushed] = mother->GetPdgCode(); 
// 	    pxPushed[npartPushed] = mother->Px(); 
// 	    pyPushed[npartPushed] = mother->Py(); 
// 	    pzPushed[npartPushed] = mother->Pz();
// 	    ePushed[npartPushed] = mother->Energy(); 
// 	    npartPushed++; 
// 	    if (isMu0Acc && (idDec[iproc] < 2 || mother->GetPdgCode() > 0)) { 
// 	      pdgPushed[npartPushed] = fMu[0]->GetPdgCode(); 
// 	      pxPushed[npartPushed] = fMu[0]->Px(); 
// 	      pyPushed[npartPushed] = fMu[0]->Py(); 
// 	      pzPushed[npartPushed] = fMu[0]->Pz();
// 	      ePushed[npartPushed] = fMu[0]->Energy(); 
// 	      npartPushed++; 
// 	      nmuons++; 
// 	    }
	    
// 	    if (isMu1Acc && (idDec[iproc] < 2 || mother->GetPdgCode() < 0)) { 
// 	      pdgPushed[npartPushed] = fMu[1]->GetPdgCode(); 
// 	      pxPushed[npartPushed] = fMu[1]->Px(); 
// 	      pyPushed[npartPushed] = fMu[1]->Py(); 
// 	      pzPushed[npartPushed] = fMu[1]->Pz();
// 	      ePushed[npartPushed] = fMu[1]->Energy(); 
// 	      npartPushed++; 
// 	      nmuons++; 
// 	    }
// 	  }
	  
	} // end if BR
      } // end loop on multiplicity 
    }  // end loop on process 
    fHMultMu->Fill(nmuons); 
    //  } // keep on generating until at least a muon is created in the event
  
//   Int_t ntmother = 0, ntchild =0; 
//   for (Int_t ipart = 0; ipart < npartPushed; ipart++) { 
//     if (TMath::Abs(pdgPushed[ipart]) != 13) { // particle is not a muon, hence it's a mother
//       PushTrack(0,-1,pdgPushed[ipart],
// 		pxPushed[ipart],pyPushed[ipart],pzPushed[ipart],ePushed[ipart],
// 		origin0[0],origin0[1],origin0[2],0.,
// 		polar[0],polar[1],polar[2],
// 		kPPrimary,ntmother,1,11);
//       KeepTrack(ntmother); 
//     }
//     else { 
//       PushTrack(1,ntmother,pdgPushed[ipart],
// 		pxPushed[ipart],pyPushed[ipart],pzPushed[ipart],ePushed[ipart],
// 		origin0[0],origin0[1],origin0[2],0.,
// 		polar[0],polar[1],polar[2],
// 		kPDecay,ntchild,1,1);
//       KeepTrack(ntchild); 
//     }
//   }
//   SetHighWaterMark(ntchild); 
//   GenEventHeader* header = new GenEventHeader("LMR");
//   header->SetPrimaryVertex(fVertex);
//   header->SetNProduced(fNprimaries);
//   AddHeader(header); 
}

//------------------------------------------------------------------

void GenMUONLMR::Decay2Body(TParticle *mother){ 
  // performs decay in two muons of the low mass resonances
  Double_t md1 = fMu[0]->GetMass(); 
  Int_t pdg = mother->GetPdgCode(); 
  Double_t mres =0; 
  // if mother is a rho, extract the mass from its line shape
  // otherwise consider the resonance mass 
  if (pdg == 113) mres = fRhoLineShape->GetRandom(); 
  else mres = mother->GetCalcMass(); 
  //  while (mres < md1 + md2) mres =  fDsigmaDm[res]->GetRandom();
  // energies and momenta in rest frame 
  Double_t e1 = mres / 2.;
  Double_t p1 = TMath::Sqrt((e1 + md1)*(e1 - md1)); 
  // orientation in decaying particle rest frame
  Double_t costheta = gRandom->Rndm() * 2 - 1;
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2. * TMath::Pi() * gRandom->Rndm(); 
  Double_t px1      = p1 * sintheta * TMath::Cos(phi); 
  Double_t py1      = p1 * sintheta * TMath::Sin(phi); 
  Double_t pz1      = p1 * costheta; 

  // boost muons into lab frame 

  TLorentzVector vmother, v1, v2;
  //  TLorentzVector boosted1, boosted2;   
  vmother.SetPxPyPzE(mother->Px(),mother->Py(),mother->Pz(),mother->Energy());
  v1.SetPxPyPzE(px1,py1,pz1,e1); 
  v2.SetPxPyPzE(-px1,-py1,-pz1,e1); 

  TVector3 betaParent = (1./vmother.E())*vmother.Vect(); // beta = p/E
  v1.Boost(betaParent);
  v2.Boost(betaParent);

  fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  fMu[1]->SetMomentum(v2.Px(),v2.Py(),v2.Pz(),v2.E());
} 

//------------------------------------------------------------------

void GenMUONLMR::DecayPiK(TParticle *mother, Bool_t &hasDecayed){ 
  // performs decays of pions and kaons
  Double_t md1 = fMu[0]->GetMass(); 
  // extract the mass from the resonance's line shape
  Double_t mres = mother->GetMass(); 
  // choose the pi/k sign, assuming 50% probabilities for both signs
  Int_t sign = (gRandom->Rndm() > 0.5) ? 1 : -1;
  mother->SetPdgCode(sign * TMath::Abs(mother->GetPdgCode())); 

  // energies and momenta in rest frame 
  Double_t e1 = (mres*mres + md1*md1)/(2*mres);
  Double_t p1 = TMath::Sqrt((e1 + md1)*(e1 - md1)); 
  // orientation in decaying particle rest frame
  Double_t costheta = gRandom->Rndm() * 2 - 1;
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2. * TMath::Pi() * gRandom->Rndm(); 
  Double_t px1      = p1 * sintheta * TMath::Cos(phi); 
  Double_t py1      = p1 * sintheta * TMath::Sin(phi); 
  Double_t pz1      = p1 * costheta;  

  // boost muons into lab frame 
  TLorentzVector vmother, v1;
  vmother.SetPxPyPzE(mother->Px(),mother->Py(),mother->Pz(),mother->Energy());
  v1.SetPxPyPzE(px1,py1,pz1,e1); 

  TVector3 betaParent = (1./vmother.E())*vmother.Vect(); // beta = p/E
  v1.Boost(betaParent);  
  if (mother->GetPdgCode()>0) fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  else fMu[1]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());

  Int_t idmother = -1; 
  if (TMath::Abs(mother->GetPdgCode())== 211) idmother = 0; 
  if (TMath::Abs(mother->GetPdgCode())== 321) idmother = 1; 
  Double_t gammaRes = mother->Energy()/mres;
  Double_t zResCM = fDecay[idmother]->GetRandom();
  Double_t zResLab = gammaRes*zResCM;  
  if(zResLab > 0.938) hasDecayed = 0; // 0.938: distance from IP to absorber + lambda_i
  else hasDecayed = 1;

} 

//-------------------------------------------------------------------

void GenMUONLMR::DalitzDecay(TParticle *mother){
  //
  // perform dalitz decays of eta, omega and etaprime 
  //
  //in the rest frame of the virtual photon:
  Double_t mres = mother->GetCalcMass(); 
  Double_t mumass  = fMu[0]->GetMass(); 
  Double_t md3  = 0;  // unless differently specified, third particle is a photon 
  if (mother->GetPdgCode() == 223) md3 = 0.134977; // if mother is an omega, third particle is a pi0
  Int_t index = -1; 
  if (mother->GetPdgCode() == 221) index = 0;  // eta
  else if (mother->GetPdgCode() == 223) index = 1; // omega  
  else if (mother->GetPdgCode() == 331) index = 2; // etaPrime  
  Int_t flag = 0; 
  Double_t xd=0, mvirt2=0; 
  Double_t countIt = 0;
  while (flag==0) {  
    xd       = fDalitz[index]->GetRandom(); 
    mvirt2   = xd * mres * mres;   // mass of virtual photon 
    // check kinematics 
    if (mres - md3 > TMath::Sqrt(mvirt2) && TMath::Sqrt(mvirt2)/2. > mumass) flag=1;
    if (++countIt>1E11) {
      mvirt2 =  mres * mres * 0.998; 
      break;
    }
  }  
 
  //
  //        Generate muons in virtual photon rest frame. 
 
  //

  Double_t e1 = TMath::Sqrt(mvirt2)/2.; // energy of mu1 in the virtual photon frame
  Double_t psquare = (e1 + mumass)*(e1 - mumass); 
  if (psquare<0) {
    printf("Error in GenMUONLMR::DalitzDecay: sqrt of psquare = %f put to 0\n",psquare); 
    psquare = 0;
  }
  Double_t p1 = TMath::Sqrt(psquare);
  
  //theta angle between the pos. muon and the resonance 
  // function for polarized theta distributions
  fCosTheta = new TF1("fCosTheta",CosThetaParam,-1,1,1);
  fCosTheta-> SetParameter(0,TMath::Sqrt(mvirt2));
  Double_t costheta = fCosTheta->GetRandom();
  if (costheta>1)  costheta = 1; 
  if (costheta<-1) costheta = -1; 
  Double_t sintheta = TMath::Sqrt((1. + costheta)*(1. - costheta));
  Double_t phi      = 2 * TMath::Pi() * gRandom->Rndm();
  Double_t sinphi   = TMath::Sin(phi);
  Double_t cosphi   = TMath::Cos(phi);

  // fill 4-vectors of leptons in the virtual photon frame

  Double_t px1 = p1*sintheta*cosphi; 
  Double_t py1 = p1*sintheta*sinphi; 
  Double_t pz1 = p1*costheta; 
  Double_t px2 = -p1*sintheta*cosphi; 
  Double_t py2 = -p1*sintheta*sinphi; 
  Double_t pz2 = -p1*costheta; 
  Double_t e2  = e1; 

  TLorentzVector v1(px1,py1,pz1,e1); 
  TLorentzVector v2(px2,py2,pz2,e2);
  fMu[0]->SetMomentum(px1,py1,pz1,e1); 
  fMu[1]->SetMomentum(px2,py2,pz2,e2); 

  // calculate components of non-dilepton in CMS of parent resonance 

  Double_t e3 = (mres * mres + md3 * md3 - mvirt2) / (2.*mres);
  Double_t psquare3 = (e3 + md3)*(e3 - md3); 
  if (psquare3<0) {
    printf("Error in GenMUONLMR::DalitzDecay: sqrt of psquare3 = %f put to 0\n",psquare3); 
    psquare3 = 0;
  }
  Double_t p3 = TMath::Sqrt(psquare3);
  Double_t costheta2 = 2.* gRandom->Rndm() - 1.;   // angle between virtual photon and resonance
  if (costheta2>1)  costheta2 = 1; 
  if (costheta2<-1) costheta2 = -1; 
  Double_t sintheta2 = TMath::Sqrt((1. + costheta2)*(1. - costheta2));
  Double_t phi2      = 2 * TMath::Pi() * gRandom->Rndm();
  Double_t sinphi2   = TMath::Sin(phi2);
  Double_t cosphi2   = TMath::Cos(phi2);
  Double_t px3 = p3*sintheta2*cosphi2; 
  Double_t py3 = p3*sintheta2*sinphi2; 
  Double_t pz3 = p3*costheta2; 
  TLorentzVector v3(px3,py3,pz3,e3); 
  
  Double_t evirt = mres - e3; 
  Double_t pxvirt = -px3;
  Double_t pyvirt = -py3;
  Double_t pzvirt = -pz3;
  TLorentzVector vvirt(pxvirt,pyvirt,pzvirt,evirt); 

  // rotate lepton vectors into coordinate system of virtual photon
  TVector3 v3virt = (vvirt.Vect()).Unit(); 
  v1.RotateUz(v3virt);
  v2.RotateUz(v3virt);

  TVector3 betaVirt = (1./evirt) * vvirt.Vect(); // virtual photon beta in res frame


  // boost the muons in the frame where the resonance is at rest 

  v1.Boost(betaVirt); 
  v2.Boost(betaVirt); 

  // boost muons and third particle in lab frame

  TLorentzVector vmother(mother->Px(), mother->Py(), mother->Pz(), mother->Energy());  
  TVector3 resBetaLab = (1./vmother.E())*vmother.Vect(); // eta beta in lab frame
  v1.Boost(resBetaLab); 
  v2.Boost(resBetaLab); 
  v3.Boost(resBetaLab); 
  vvirt.Boost(resBetaLab); 

  fMu[0]->SetMomentum(v1.Px(),v1.Py(),v1.Pz(),v1.E());
  fMu[1]->SetMomentum(v2.Px(),v2.Py(),v2.Pz(),v2.E());
//   part3->SetMomentum(v3.Px(),v3.Py(),v3.Pz(),v3.E());


}

//------------------------------------------------------------------

Double_t GenMUONLMR::FormFactor(Double_t q2, Int_t decay){ 
  //  Calculates the form factor for Dalitz decays A->B+l+l
  //  Returns: |F(q^2)|^2
  //
  //  References: L.G. Landsberg, Physics Reports 128 No.6 (1985) 301-376. 
 
  Double_t ff2, mass2;
  Double_t n2, n4, m2; 
  // Lepton-G
  
  Double_t lambda2inv = 0;
  switch (decay) { 
  case 0:   // eta -> mu mu gamma  
  // eta   -> l+ l- gamma: pole approximation
    lambda2inv = 1.95; 
    mass2 = fParticle[kEtaLMR]->GetMass() * fParticle[kEtaLMR]->GetMass(); 
    if (q2 < mass2) ff2 = TMath::Power(1./(1.-lambda2inv*q2),2);
    else ff2 = 0; 
    break;
  case 1:   // omega -> mu mu pi0 
    // omega -> l+ l- pi0: pole approximation
    mass2 = fParticle[kOmegaLMR]->GetMass() * fParticle[kOmegaLMR]->GetMass(); 
    lambda2inv = 2.26; 
    if (q2 < mass2) ff2 = TMath::Power(1./(1.-lambda2inv*q2),2);
    else ff2 = 0; 
    break;
  case 2:   // etaPrime -> mu mu gamma 
    mass2 = fParticle[kEtaPrimeLMR]->GetMass() * fParticle[kEtaPrimeLMR]->GetMass(); 
    // eta'  -> l+ l- gamma: Breit-Wigner fitted to data
    n2 = 0.764 * 0.764; 
    n4 = n2 * n2; 
    m2 = 0.1020 * 0.1020;
    if (q2 < mass2) ff2 = n4 / (TMath::Power(n2-q2,2) + m2 * n2); 
    else ff2 = 0; 
    break;
  default:
    printf ("FormFactor: Decay not found\n"); 
    return 0; 
    break; 
  }
  return ff2; 
}

//____________________________________________________________

Double_t GenMUONLMR::RhoLineShapeNew(Double_t *x, Double_t *para){
  //new parameterization implemented by Hiroyuki Sako (GSI)
  Double_t mass = *x;
  double r, GammaTot;
  Double_t mRho    = TDatabasePDG::Instance()->GetParticle("rho0")->Mass();
  Double_t mPi     = TDatabasePDG::Instance()->GetParticle("pi0")->Mass();
  Double_t mMu     = TDatabasePDG::Instance()->GetParticle("mu-")->Mass();
  Double_t Gamma0  = TDatabasePDG::Instance()->GetParticle("rho0")->Width();

  const double Norm = 0.0744416*1.01;  

  // 0.0744416 at m = 0.72297
  // is the max number with Norm=1 (for rho)
  
  double mThreshold = 2.*mPi;

  const double T = 0.170; // Assumption of pi+ temperature [GeV/c^2]
  //const double T = 0.11; // Taken from fit to pi+ temperature [GeV/c^2]
  // with Reference: LEBC-EHS collab., Z. Phys. C 50 (1991) 405

  if (mass < mThreshold) {
    r = 0.;
    return r;
  }

  double k = sqrt(0.25*mass*mass-(mThreshold/2)*(mThreshold/2));
  double k0 = sqrt(0.25*mRho*mRho-(mThreshold/2)*(mThreshold/2));

  GammaTot = (k/k0)*(k/k0)*(k/k0)*(mRho/mass)*(mRho/mass)*Gamma0;

  double FormFactor2 = 1/((mass*mass-mRho*mRho)*(mass*mass-mRho*mRho)+
			  mass*mass*GammaTot*GammaTot);

  r = pow(mass,1.5)*pow((1-mThreshold*mThreshold/(mass*mass)),1.5)*
    ((mass*mass+2*mMu*mMu)/(mass*mass))*(pow((mass*mass-4*mMu*mMu),0.5)/mass)*FormFactor2
    *exp(-mass/T)/Norm;

  return r;
}
