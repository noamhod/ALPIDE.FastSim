#ifndef GenMUONLMR_h
#define GenMUONLMR_h

#include <TH1F.h> 
#include <TH1D.h> 
#include <TF1.h> 
#include <TParticle.h> 
#include <TLorentzVector.h> 
 
 
class GenMUONLMR  { 
 public:
  enum parttype_t {kPionLMR, kKaonLMR, kEtaLMR, kRhoLMR, kOmegaLMR, kPhiLMR, kEtaPrimeLMR, kJPsi,kJPsiSch};
  GenMUONLMR(Double_t energy=7.0, Int_t kLowEnergy = kFALSE); 
  GenMUONLMR(GenMUONLMR &gen); 
  GenMUONLMR &operator=(const GenMUONLMR &gen);  
  virtual ~GenMUONLMR(); 
  static Double_t PtDistr(Double_t *x, Double_t *par); 
  static Double_t YDistr(Double_t *x, Double_t *par); 
  static Double_t PtDistrLowEnergy(Double_t *x, Double_t *par); 
  static Double_t YDistrLowEnergy(Double_t *x, Double_t *par); 
  static Double_t PtDistrPythia6(Double_t *x, Double_t *par);  
  static Double_t YDistrSch(Double_t *x, Double_t *par);  
  static Double_t CosThetaParam(Double_t *x, Double_t *par); 
  void SetPtParams(Int_t iproc, Double_t p1, Double_t p2, Double_t p3, Double_t p4) {fPt[iproc]->SetParameters(p1,p2,p3,p4);}
  void SetYParams(Int_t iproc, Double_t p1, Double_t p2, Double_t p3, Double_t p4, Double_t p5=0.) {fY[iproc]->SetParameters(p1,p2,p3,p4,p5);}
  void Decay2Body(TParticle *mother);
  void DalitzDecay(TParticle *mother);
  void DecayPiK(TParticle *mother, Bool_t &hadDecayed);
  Double_t FormFactor(Double_t q2, Int_t decay); 
  void Generate(); 
  TParticle* GetMuon(Int_t i) {return fMu[i];} 
  void SetNMuMin(Int_t nmin) {fNMuMin = nmin;}
  void GenerateSingleProcess(Int_t whichproc) { fGenSingleProc = whichproc;}
  void SetScaleMultiplicity(Int_t ipart, Double_t scale) { fScaleMult[ipart] = scale; } 
  static Double_t RhoLineShapeNew(Double_t *, Double_t *); 
  void FinishRun(); 
  TF1* GetRapidity(Int_t iproc) { return fY[iproc]; }
  TF1* GetPt(Int_t iproc) { return fPt[iproc]; }
  TF1* GetMultiplicity(Int_t iproc) { return fMult[iproc]; }
 private: 
  static const Int_t fgkNpart = 9; // number of particles to be generated 
  Int_t fNMuMin;                   // min. number of muons to accept the event for writing
  Int_t fGenSingleProc;            // flag to generate a single process (1) or the whole cocktail (0)
  Int_t fPDG[fgkNpart];                   // pdg code of particle to be generated 
  Double_t fScaleMult[fgkNpart];          // multiplicity scaling factor (w.r.t. pythia@7TeV)
  TF1 *fPt[fgkNpart];                     // pt distribution
  TF1 *fY[fgkNpart];                      // rapidity distribution
  TF1 *fMult[fgkNpart];                   // multiplicity distribution 
  TF1 *fDecay[2];                  // fDecay[0] = pion, fDecay[1] = kaon
  TH1F *fDalitz[3];                // Dalitz decay form factor for eta, omega, etaprime
  TF1 *fCosTheta;                  // function for polarized theta distributions
  TF1 *fRhoLineShape;              // rho line shape 
  TParticle* fParticle[fgkNpart];         // TPaticle object for the particles to be generated
  TParticle* fMu[2];               // fMu[0] = mu+    fMu[1] = mu-
  TH1D *fHMultMu;                  // muon multiplicity 
  TH1D *fHNProc;                   // number of events generated per process
  ClassDef(GenMUONLMR, 1)       // low mass dimuons parametric generator
}; 

#endif
