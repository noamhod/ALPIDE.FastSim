#ifndef KMCPROBEFWD_H
#define KMCPROBEFWD_H

#include "TrackPar.h"

class KMCLayerFwd;




class KMCProbeFwd: public TObject {

 public:
  enum {kBitKilled=BIT(14)};
  enum {kNDOF=5,kMaxITSLr=32};
  enum {kY2=0,kZ2=2,kSnp2=5,kTgl2=9,kPtI2=14};
  enum {kY,kZ,kSnp,kTgl,kPtI};
  //
  KMCProbeFwd();
  KMCProbeFwd(double *xyz, double *pxyz, Int_t sign, double errLoose=-1);  
  KMCProbeFwd(const KMCProbeFwd& src);
  KMCProbeFwd& operator=(const KMCProbeFwd& src);
  void ImposeKinematics(const double* xyzLab,const double* cosinesLab, double en, double mass, int charge);
  //
  Int_t GetTrID()     const   {return int(GetUniqueID())-1;}
  void  SetTrID(int id)       {SetUniqueID(id+1);}
  //
  TrackPar* GetTrack()    const {return (TrackPar*)&fTrack;}
  Double_t* GetCovariance()            const {return (double*)fTrack.GetCovariance();}
  //
  void      Reset(); 
  void      ResetCovariance(float err=-1.);
  Bool_t    Init(const double *xyz, const double *pxyz, Int_t sign, double errLoose=-1);
  virtual   Bool_t IsSortable()                     const {return kTRUE;}
  virtual   Int_t  Compare(const TObject* obj)      const;
  void      Kill(Bool_t v=kTRUE)                          {SetBit(kBitKilled,v);}
  Bool_t    IsKilled()                              const {return TestBit(kBitKilled);}
  //
  void      SetWeight(double w=1.)                      {fWeight = w;}
  void      SetMass(double m=0.14)                      {fMass = m;}
  Double_t  GetMass()                             const {return fMass;}
  Double_t  GetWeight()                           const {return fWeight;}
  Double_t  GetChi2()                             const {return fChi2;}
  Double_t  GetChi2ITS()                          const {return fChi2ITS;}

  void      SetChi2(double chi2)                        {fChi2 = chi2;}
  void      SetChi2ITS(double chi2)                     {fChi2ITS = chi2;}
  void      AddChi2(double chi2)                        {fChi2 += chi2;}
  void      AddChi2ITS(double chi2)                     {fChi2ITS += chi2;}
  void      SetInnerLrChecked(Int_t n)                  {if (n<fgNITSLayers) fInnLrCheck = n;}
  Int_t     GetNITSHits()                         const {return fNHitsITS;}
  Int_t     GetNMSHits()                          const {return fNHitsMS;}
  Int_t     GetNTRHits()                          const {return fNHitsTR;}
  Int_t     GetNFakeITSHits()                     const {return fNHitsITSFake;}
  Int_t     GetNHits()                            const {return fNHits;}
  Int_t     GetInnerLayerChecked()                const {return fInnLrCheck;}
  UInt_t&   GetHitsPatt()                               {return fHits;}
  UInt_t&   GetFakesPatt()                              {return fFakes;}

  Double_t  GetNormChi2(Bool_t penalize=kFALSE)      const;
  Double_t  GetNormChi2ITS(Bool_t penalize=kFALSE)   const;
  void      AddHit(const KMCLayerFwd*lr , double chi2, Int_t clID=-1);
  void      ResetHit(Int_t lr);
  Bool_t    IsHit(Int_t lr)                       const {return (lr<fgNITSLayers) ? IsWBit(fHits,lr)  : kFALSE;}
  Bool_t    IsHitFake(Int_t lr)                   const {return (lr<fgNITSLayers) ? IsWBit(fFakes,lr) : kFALSE;}

  Bool_t   ApplyMSEL(double x2X0, double xTimesRho);
  Bool_t   CorrectForMeanMaterial(double xOverX0, double xTimesRho, Bool_t modeMC=kFALSE, Bool_t anglecorr=kFALSE);
  Bool_t   PropagateToZBxByBz(double z, double maxDZ=1.0, Double_t xOverX0=0., Double_t xTimesRho=0., Bool_t modeMC=kFALSE);
  Bool_t   PropagateToZBxByBz(double z, const double *bxyz);
  Bool_t   PropagateToDCA(KMCProbeFwd* partner);
  Bool_t   Update(Double_t cov[3]);

  Double_t GetR()                        const {double x=GetX(),y=GetY(),r=x*x+y*y; return r>0?TMath::Sqrt(r):0;}
  Double_t GetX()                        const {return  fTrack.GetY();}
  Double_t GetY()                        const {return  fTrack.GetZ();} // Y of local frame
  Double_t GetZ()                        const {return NegDir() ? -fTrack.GetX():fTrack.GetX();}
  void     GetXYZ(double *xyz)           const;
  void     GetPXYZ(double *pxyz)         const;  
  double   GetP()                        const {return fTrack.GetP();}
  Double_t GetSigmaX2()                  const {return fTrack.GetSigmaY2();}
  Double_t GetSigmaY2()                  const {return fTrack.GetSigmaZ2();}
  Double_t GetSigmaXY()                  const {return fTrack.GetSigmaZY();}  
  Double_t GetSigmaP2()                  const;
  Double_t GetSigmaPX2()                 const;
  Double_t GetSigmaPY2()                 const;
  Double_t GetSigmaPZ2()                 const;
  Double_t GetXLoc()                     const {return fTrack.GetX();}
  Double_t GetYLoc()                     const {return fTrack.GetY();}
  Double_t GetZLoc()                     const {return fTrack.GetZ();}
  Double_t GetPredictedChi2(Double_t* p, Double_t* cov) const {return fTrack.GetPredictedChi2(p,cov);}
  Double_t GetAlpha()                    const {return fTrack.GetAlpha();}
  Double_t GetCharge()                   const {return fTrack.Charge();}

  //
  static void   SetWBit(UInt_t &patt,UInt_t bit)               {patt |= 0x1<<bit;}
  static void   ResetWBit(UInt_t &patt,UInt_t bit)             {patt &= ~(0x1<<bit);}
  static Bool_t IsWBit(const UInt_t &patt,const UInt_t bit)    {return patt&(0x1<<bit);}
  static void   SetNITSLayers(Int_t n)                         {fgNITSLayers = n;}
  static int    GetNITSLayers()                                {return fgNITSLayers;}
  //
  static Double_t GetMissingHitPenalty()                        {return fgMissingHitPenalty;}
  static void     SetMissingHitPenalty(double p=2.)             {fgMissingHitPenalty = p;}
  //
  static void Lab2Trk(const double *vLab, double *vTrk); 
  static void Trk2Lab(const double *vTrk, double *vLab); 
  //
  // protected:
  Bool_t Update(Double_t p[2],Double_t cov[3]);
  Bool_t NegDir()                     const {return TMath::Abs(fTrack.GetAlpha())>TMath::Pi()/2;}
  Bool_t IsZero(double val, double tol=1e-9) const {return TMath::Abs(val)<tol;}
  virtual void  Print(Option_t* option = "") const;
  //
 protected:
  Double_t fWeight; // weight (for decay?)
  Double_t fMass;   // particle mass
  Double_t fChi2;   // total chi2
  Double_t fChi2ITS;// total chi2 ITS
  UInt_t   fHits;   // pattern on hits (max 32!)
  UInt_t   fFakes;  // pattern of fakes among hits
  Int_t    fNHits;    // total hits
  Int_t    fNHitsITS; // total ITS hits
  Int_t    fNHitsMS;  // total MS hits
  Int_t    fNHitsTR;  // total TR hits
  Int_t    fNHitsITSFake; // number of fake ITS hits
  UInt_t   fInnLrCheck;   // lowest active layer where update was checked
  Int_t    fClID[kMaxITSLr]; // id's of attached clusters
  TrackPar fTrack;  // track params
  //
  static Int_t    fgNITSLayers;
  static Double_t fgMissingHitPenalty;  //
  ClassDef(KMCProbeFwd,1)
};
/*
// RS: swap adapted for the field along lab X axis
//_______________________________________________________________________
inline void  KMCProbeFwd::Lab2Trk(const double *vLab, double *vTrk)
{
  // convert alice coordinates to modified
  vTrk[0] = vLab[2];
  vTrk[1] = vLab[1];
  vTrk[2] =-vLab[0];
}

//_______________________________________________________________________
inline void  KMCProbeFwd::Trk2Lab(const double *vTrk, double *vLab)
{
  // convert modified coordinates to Lab ones
  vLab[0] =-vTrk[2];
  vLab[1] = vTrk[1];
  vLab[2] = vTrk[0];
}
*/

// RS: swap adapted for the field along lab Y axis
//_______________________________________________________________________
inline void  KMCProbeFwd::Lab2Trk(const double *vLab, double *vTrk)
{
  // convert alice coordinates to modified
  vTrk[0] = vLab[2];
  vTrk[1] = vLab[0];
  vTrk[2] = vLab[1];
}

//_______________________________________________________________________
inline void  KMCProbeFwd::Trk2Lab(const double *vTrk, double *vLab)
{
  // convert modified coordinates to Lab ones
  vLab[0] = vTrk[1];
  vLab[1] = vTrk[2];
  vLab[2] = vTrk[0];
}


//_______________________________________________________________________
inline void KMCProbeFwd::GetXYZ(double *xyz) const
{
  // track position in Lab coordinates
  double xyzTrk[3];
  fTrack.GetXYZ(xyzTrk);
  Trk2Lab(xyzTrk, xyz);
}

//_______________________________________________________________________
inline void KMCProbeFwd::GetPXYZ(double *pxyz) const
{
  // track position in Lab coordinates
  double pxyzTrk[3];
  fTrack.GetPxPyPz(pxyzTrk);
  Trk2Lab(pxyzTrk, pxyz);
}

//_______________________________________________________________________
inline Bool_t KMCProbeFwd::PropagateToZBxByBz(double z, const double *bxyz)
{
  //
  //  return fTrack.PropagateToBxByBz(NegDir() ? -z:z,bxyz);
  //  printf("PropagateTo(%f, %f)\n",NegDir() ? -z:z,bxyz[2]);
  Bool_t res = fTrack.PropagateToBxByBz(z,bxyz);
  // if (!res) fTrack.Print("");
  return res;
  //
}

//_______________________________________________________________________
inline Bool_t KMCProbeFwd::Update(Double_t p[2],Double_t cov[3])
{
  // update with measurement in the plane normal to Z axis (in Lab X,Y axis) 
  return fTrack.Update(p,cov);
}

//_______________________________________________________________________
inline Bool_t KMCProbeFwd::Update(Double_t cov[3])
{
  // dummy update with measurement errors in the plane normal to Z axis (in Lab X,Y axis) 
  //  printf("Before update %f %f %f\n",cov[0],cov[1],cov[2]); fTrack.Print();
  double p[2] = {fTrack.GetY(), fTrack.GetZ()};
  if (!fTrack.Update(p,cov)) return kFALSE;
  //  printf("After update: \n"); fTrack.Print();
  return kTRUE;  
}

//_______________________________________
inline Double_t KMCProbeFwd::GetNormChi2(Bool_t penalize) const
{
  // normalized chi2, penilized for missing hits
  if (fNHits<3) return 0;
  double chi2 = fChi2;
  if (penalize) {
    int nMiss = fgNITSLayers - fNHitsITS - fInnLrCheck;
    chi2 +=  fgMissingHitPenalty*nMiss;
  }
  return chi2/( (fNHits<<1)-kNDOF);
}

//_______________________________________
inline Double_t KMCProbeFwd::GetNormChi2ITS(Bool_t penalize) const
{
  // normalized chi2, penalized for missing hits
  double chi2 = fChi2ITS;
  if (penalize) {
    int nMiss = fgNITSLayers - fNHitsITS - fInnLrCheck;
    chi2 += fgMissingHitPenalty*nMiss;
  }
  if (fNHitsITS<3) return chi2;
  //
  return chi2/( (fNHitsITS<<1)-kNDOF);
}

//_______________________________________
inline void KMCProbeFwd::ResetHit(Int_t lr) {
  // note: lr is active layer ID
  if (lr>=fgNITSLayers) return; 
  if (IsWBit(fHits,lr))  {fNHitsITS--;     ResetWBit(fHits,lr);}
  if (IsWBit(fFakes,lr)) {fNHitsITSFake--; ResetWBit(fFakes,lr);}
}

//_______________________________________
inline Bool_t KMCProbeFwd::CorrectForMeanMaterial(double xOverX0, double xTimesRho, Bool_t modeMC, Bool_t anglecorr)
{
  Bool_t res;
  if (modeMC) {
    if (!(res=ApplyMSEL(xOverX0,xTimesRho))) {
      printf("%s\n",Form("Failed to apply MS, %f",xOverX0));
      fTrack.Print();
      return kFALSE;
    }
    else xTimesRho = 0; // don't apply e.loss in the ETP
  }
  res = fTrack.CorrectForMeanMaterial(xOverX0, xTimesRho, fMass, anglecorr);
  if (!res) {
    printf("%s\n",Form("Failed to CorrectForMeanMaterial: %f %f",xOverX0,xTimesRho));
    fTrack.Print();
    return res;
  }
  return kTRUE;
}

#endif




