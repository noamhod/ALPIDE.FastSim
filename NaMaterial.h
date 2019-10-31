#ifndef NAMAERIAL_H
#define NAMAERIAL_H

#include "TMaterial.h"

//====================================================================
class NaMaterial :public TMaterial {
public:
  enum {kNELossPar=2};
  NaMaterial();
  NaMaterial(const char* name, const char *title, 
	     Float_t a, Float_t z, Float_t dens, Float_t radl=0, 
	     Float_t absl=0, Float_t *elbuff=0);
  virtual ~NaMaterial();
  //
  virtual void Dump() const;
  virtual void Print(Option_t* option="") const;
  //
  Bool_t IsELossProvided() const {return fELossPar[0]>0;}
  Float_t GetELoss(Float_t p, float mass) const;
  Float_t GetELoss2ETP(Float_t p, float mass) const;
  
  const Float_t *GetELossPars() const {return fELossPar;}
  Float_t GetELossPar(int i) const {return (i<kNELossPar && i>=0) ? fELossPar[i]:0.;}
  //
 protected:
  Float_t fELossPar[kNELossPar];  // Params: mean excitation energy I and plasma energy (GeV): http://pdg.lbl.gov/2018/AtomicNuclearProperties/HTML/
  //
  ClassDef(NaMaterial,1) // Material Object
};

//==========================================================================
class NaMixture :public NaMaterial {
public:
  NaMixture();
  NaMixture(const char* name, const char *title,
	    Float_t a, Float_t z, Float_t dens,
	    Float_t radl=0, Float_t absl=0, Float_t *elbuff=0);
  virtual void SetComponents(Int_t nmixt, Float_t* a,Float_t* z,Float_t* w);
  virtual ~NaMixture();
  //
  virtual void Print(Option_t* option="") const;
  //
 protected:
  Int_t fNMix;    // number of components to mix (<>0 a la Geant3 mixture)
  Float_t* fAMix; // [fNMix] A of components
  Float_t* fZMix; // [fNMix] Z ..
  Float_t* fWMix; // [fNMix] Weights ...
  //
  ClassDef(NaMixture,1) // Mixture Class
};


#endif
