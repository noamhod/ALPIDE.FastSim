#include "NaMaterial.h"
#include "TrackPar.h"

//==========================================================================

ClassImp(NaMaterial)

NaMaterial::NaMaterial()
{ 
// Defaul Constructor
}
//----------------------------------------------------------------
//
NaMaterial::NaMaterial(const char* name, const char* title, 
		       Float_t a, Float_t z, Float_t dens, Float_t radl, 
		       Float_t absl, Float_t *elbuff) 
  : TMaterial(name,title,a,z,dens,radl,absl) 
{
  // Constructor
  if (elbuff) // ELoss Parameters are provided
    for (int i=0;i<kNELossPar;i++) fELossPar[i] = elbuff[i];
  else
    for (int i=0;i<kNELossPar;i++) fELossPar[i] = -1.;
  //
}
//----------------------------------------------------------------
//
NaMaterial::~NaMaterial() {}

//
//----------------------------------------------------------------
//
void NaMaterial::Print(Option_t* option) const
{
//Print material information
  Float_t a=0,z=0,rho=0,rl=0,il=0;
  //  
  a = GetA(); z = GetZ(); rho = GetDensity(); 
  rl = GetRadLength(); il = GetInterLength();
  //
  if (option[0] != 'h' && option[0] != 'H' ) //header is not suppressed
    printf("%-20s %7s %7s %7s %10s %10s %10s %10s %10s\n",
	   "Material","   A   ","   Z   ","Density","  Rad.L  ",
	   "  Inter.L "," ELossC0 "," ELossC1 "," ELossC2 ");
  printf("%-20s %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e\n",
	 GetName(),a,z,rho,rl,il,fELossPar[0],
	 fELossPar[1],fELossPar[2]);
  //
}
//----------------------------------------------------------------
//
Float_t NaMaterial::GetELoss(Float_t p, float m) const
{
  // return energy lost (in MeV cm^2/g) at momentum p for mass p
  float bg = p/m;
  if (!IsELossProvided()) {
    return TrackPar::BetheBlochSolid(bg);
  }
  // real Bethe-Bloch
  float b2g2 = bg*bg, b2 = b2g2/(1.+b2g2), g = TMath::Sqrt(b2g2/b2);
  const float me = 0.511e-3;
  const float k = 0.307075e-3;
  float me2m = me/m;
  float T = (2*me*b2g2) / (1 + me2m*(2*g + me2m) );
  float densEff = 0.f; // density effect
  float plasmaE = GetELossPars()[1], meanExcitE = GetELossPars()[0];
  if (plasmaE>0) { 
    densEff = TMath::Log(plasmaE/meanExcitE) + 0.5*TMath::Log(b2g2) - 0.5;
  }
  float SP = ((k*GetZ())/(GetA()*b2))*(0.5*TMath::Log( (2*me*b2g2*T) / (meanExcitE*meanExcitE) ) - b2 - densEff);
  return SP;
  
}
//
Float_t NaMaterial::GetELoss2ETP(Float_t p, float m) const
{
  // return ratio between real BB and TrackPar energy loss 
  float bg = p/m;
  if (!IsELossProvided()) {
    return 1.;
  }
  // real Bethe-Bloch
  float b2g2 = bg*bg, b2 = b2g2/(1.+b2g2), g = TMath::Sqrt(b2g2/b2);
  const float me = 0.511e-3;
  const float k = 0.307075e-3;
  float me2m = me/m;
  float T = (2*me*b2g2) / (1 + me2m*(2*g + me2m) );
  float densEff = 0.f; // density effect
  float plasmaE = GetELossPars()[1], meanExcitE = GetELossPars()[0];
  if (plasmaE>0) { 
    densEff = TMath::Log(plasmaE/meanExcitE) + 0.5*TMath::Log(b2g2) - 0.5;
  }
  float SP = ((k*GetZ())/(GetA()*b2))*(0.5*TMath::Log( (2*me*b2g2*T) / (meanExcitE*meanExcitE) ) - b2 - densEff);
  return SP/TrackPar::BetheBlochSolid(bg);
  
}

//----------------------------------------------------------------
//
void NaMaterial::Dump() const 
{
  Print();
}
//----------------------------------------------------------------
//

//==========================================================================

ClassImp(NaMixture)

NaMixture::NaMixture()
{ 
// Defaul Constructor
  fNMix = 0;
  fAMix = 0;
  fZMix = 0;
  fWMix = 0;
  //
}
//----------------------------------------------------------------
//
NaMixture::NaMixture(const char* name, const char* title, 
		     Float_t a, Float_t z, Float_t dens,
		     Float_t radl,Float_t absl, Float_t *elbuff) 
  : NaMaterial(name,title,a,z,dens,radl,absl,elbuff) 
{
  // Constructor
  fNMix = 0;
  fAMix = fZMix = fWMix = 0;
}
//----------------------------------------------------------------
//
NaMixture::~NaMixture() 
{
  if (fAMix) delete[] fAMix;
  if (fZMix) delete[] fZMix;
  if (fWMix) delete[] fWMix;
}
//

void NaMixture::SetComponents(Int_t nmixt, Float_t* a,Float_t* z,Float_t* w)
{
  // Sets the components of the mixture (Geant3 conventions preserved)
  if (fAMix || fZMix || fWMix) { 
    Error("SetComponents","Components arrays were already initialized %s",GetName());
    return;
  }
  fNMix = TMath::Abs(nmixt);
  if ( fNMix<1 ) {
    Error("SetComponents","Number of components is 0 for %s",GetName());
    return;
  }
  fAMix = new Float_t[fNMix];
  fZMix = new Float_t[fNMix];
  fWMix = new Float_t[fNMix];
  Float_t amol = 0.;
  for (int i=0;i<fNMix;i++) {
    fAMix[i] = a[i];
    fZMix[i] = z[i];
    fWMix[i] = w[i];
    amol += fWMix[i]*fAMix[i];
  }
  if (amol<=0.) {Error("SetComponents","total weigth for %s is <=0",GetName()); return;}
  if (nmixt<0 ) for (int i=0;i<fNMix;i++) fWMix[i] *= fAMix[i]/amol; //use 'proportion by weight'
  //
}

//----------------------------------------------------------------
//
void NaMixture::Print(Option_t* option) const
{
//Print mixture information
  Float_t a=0,z=0,rho=0,rl=0,il=0;
  //  
  a = GetA(); z = GetZ(); rho = GetDensity(); 
  rl = GetRadLength(); il = GetInterLength();
  //
  if (option[0] != 'h' && option[0] != 'H' ) //header is not suppressed
    printf("%-20s %7s %7s %7s %10s %10s %10s %10s %10s\n",
	   "Material","   A   ","   Z   ","Density","  Rad.L  ",
	   "  Inter.L "," ELossC0 "," ELossC1 "," ELossC2 ");
  printf("%-20s %7.3f %7.3f %7.3f %10.3e %10.3e %10.3e %10.3e %10.3e %6s %6s %6s\n",
	 GetName(),a,z,rho,rl,il,fELossPar[0],
	 fELossPar[1],fELossPar[2],"  A  ","  Z  ","  W  ");
  for (int i=0;i<TMath::Abs(fNMix);i++) 
    printf("%107s %6.2f %6.2f %6.4f\n"," ",fAMix[i],fZMix[i],fWMix[i]);
  //
}
//----------------------------------------------------------------
//
