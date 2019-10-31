#ifndef KMCLAYERFWD_H
#define KMCLAYERFWD_H

#include <TNamed.h>
#include <TClonesArray.h>
#include "KMCClusterFwd.h"
#include "KMCProbeFwd.h"
#include "NaMaterial.h"

class KMCLayerFwd : public TNamed {
public:
  enum {kTypeNA=-1,kVTX,kITS,kMS,kTRIG,kABS,kDUMMY,kBitVertex=BIT(15)};
  enum {kMaxAccReg = 5};
  KMCLayerFwd(const char *name);
  Float_t GetZ()         const {return fZ;}
  Float_t GetRMin()      const {return fRMin[0];}
  Float_t GetRMax()      const {return fRMax[fNAccReg-1];}
  Float_t GetRMin(int i)      const {return fRMin[i];}
  Float_t GetRMax(int i)      const {return fRMax[i];}
  int GetNAccRegions() const {return fNAccReg;}
  void SetNAccRegions(int n=1) {
    if (n<1) n=1;
    if (n>kMaxAccReg) n = kMaxAccReg;
    fNAccReg = n;
  }
  
  Float_t GetX2X0()      const {return fx2X0;}
  Float_t GetXTimesRho() const {return fXRho;}
  int GetAccRegion(float r) const {
    for (int rid=0;rid<fNAccReg;rid++) if (r>=fRMin[rid] && r<fRMax[rid]) return rid;
    return -1;
  }

  Float_t GetXResId(int i)  const {return fXRes[i];}
  Float_t GetXRes(float r=-1)  const {
    int id = r<0 ? 0 : GetAccRegion(r);
    return (id<0) ? fXRes[0] : fXRes[id];
  }

  Float_t GetYResId(int i)  const {return fYRes[i];}
  Float_t GetYRes(float r=-1)  const {
    int id = r<0 ? 0 : GetAccRegion(r);
    return (id<0) ? fYRes[0] : fYRes[id];
  }

  Float_t GetLayerEff()  const {
    return fEff;
  }
  
  Float_t GetThickness() const {return fThickness;}
  Int_t   GetActiveID()  const {return fActiveID;}
  Int_t   GetID()        const {return GetUniqueID();}
  void    SetID(int id)        {SetUniqueID(id);}
  //
  void    SetZ(Float_t v)         {fZ = v;}
  void    SetRMin(Float_t v, int i=0)      {fRMin[i] = v;}
  void    SetRMax(Float_t v, int i=0)      {fRMax[i] = v;}
  void    SetXRes(Float_t v, int i=0)      {fXRes[i] = v;}
  void    SetYRes(Float_t v, int i=0)      {fYRes[i] = v;}
  void    SetLayerEff(Float_t v)  {fEff = v;}

  void    SetX2X0(Float_t v)      {fx2X0 = v;}
  void    SetXTimesRho(Float_t v) {fXRho = v;}

  void    SetThickness(Float_t v) {fThickness = v;}
  void    SetActiveID(Int_t v)    {fActiveID = v;}
  void    SetDead(Bool_t v)       {fIsDead = v;}
  void    SetType(Int_t tp)       {fType = tp;}
  //
  Bool_t  IsDead()       const {return fIsDead;}
  Bool_t  IsITS()        const {return fType==kITS;}
  Bool_t  IsMS()         const {return fType==kMS;}
  Bool_t  IsTrig()       const {return fType==kTRIG;}
  Bool_t  IsAbs()        const {return fType==kABS;}
  Bool_t  IsDummy()      const {return fType==kDUMMY;}
  Int_t   GetType()      const {return fType;}
  Bool_t  IsVertex()     const {return TestBit(kBitVertex);}
  //
  Int_t                 AddBgCluster(double x,double y,double z, Int_t id);
  KMCClusterFwd*        GetBgCluster(Int_t i)    const {return (KMCClusterFwd*)fClBg[i];}
  TClonesArray*         GetBgClusters()          const {return (TClonesArray*)&fClBg;}
  KMCClusterFwd*        GetMCCluster()           const {return (KMCClusterFwd*)&fClMC;}
  KMCClusterFwd*        GetCorCluster()          const {return (KMCClusterFwd*)&fClCorr;}
  //
  void                  SetAnProbe(KMCProbeFwd& prb)   {fTrCorr = prb;}
  KMCProbeFwd*          GetAnProbe()             const {return (KMCProbeFwd*)&fTrCorr;}
  Int_t                 GetNMCTracks()           const {return fTrMC.GetEntries();}
  TClonesArray*         GetMCTracks()            const {return (TClonesArray*)&fTrMC;}
  KMCProbeFwd*          GetMCTrack(Int_t it)     const {return (KMCProbeFwd*)fTrMC[it];}
  Int_t                 GetNBgClusters()         const {return fClBg.GetEntries();}
  KMCProbeFwd*          AddMCTrack(KMCProbeFwd* src=0);
  KMCProbeFwd*          GetWinnerMCTrack();
  //
  Double_t              GetSig2EstX()            const {return fSig2EstX;}
  Double_t              GetSig2EstY()            const {return fSig2EstY;}
  void                  SetSig2EstX(double v)          {fSig2EstX=v;}
  void                  SetSig2EstY(double v)          {fSig2EstY=v;}
  void                  SortBGClusters()               {fClBg.Sort();}
  //
  void   Reset();
  void   ResetMC() { fTrMC.Clear();}
  void   ResetBgClusters() { fClBg.Clear(); }
  void   ResetMCTracks()   { fTrMC.Clear(); }
  virtual void  Print(Option_t* option = "") const;
  //
  void SetMaterial(NaMaterial* m) { fMaterial = m; }
  const NaMaterial* GetMaterial() const {return fMaterial;}
  float GetELoss2ETP(Float_t p, float m) const {return fMaterial ? fMaterial->GetELoss2ETP(p,m) : 1.0;}
  //
  static Double_t GetDefEff()   {return fgDefEff;}
  static void     SetDefEff(double eff=1) {fgDefEff = eff>1. ? 1.: (eff<0? 0:eff);}
  //
  //
 protected:
  //
  Float_t fZ; 
  Float_t fRMin[kMaxAccReg];
  Float_t fRMax[kMaxAccReg];
  Float_t fThickness;
  Float_t fx2X0;
  Float_t fXRho;    // x*density
  Float_t fXRes[kMaxAccReg]; 
  Float_t fYRes[kMaxAccReg];   
  Float_t fEff;
  Bool_t  fIsDead;
  Int_t   fNAccReg;
  Int_t   fType;   // its, tpc etc
  Int_t   fActiveID;   // active layer id
  Float_t fSig2EstX;
  Float_t fSig2EstY;
  //
  KMCClusterFwd   fClCorr;     // ideal cluster
  KMCClusterFwd   fClMC;       // MC cluster (from MS scattered track)
  TClonesArray    fClBg;       // bg clusters for MC
  //
  KMCProbeFwd     fTrCorr;     // ideal track
  TClonesArray    fTrMC;       // MC tracks
  //
  NaMaterial*     fMaterial;
  //
  static Double_t fgDefEff;
  ClassDef(KMCLayerFwd,2);
};

//_________________________________________________________________
inline Int_t  KMCLayerFwd::AddBgCluster(double x,double y,double z, int id)
{
  int n=GetNBgClusters(); 
  new (fClBg[n]) KMCClusterFwd(x,y,z, id); 
  return ++n;
}

//______________________________________________________
class BeamPipe : public KMCLayerFwd {
public:
  BeamPipe(char* name) :  KMCLayerFwd(name) {}
  //
  void SetRadius(float r)             {SetZ(r);}
  //
  Float_t GetInnerRadius()      const {return fZ-fThickness/2;}
  Float_t GetOuterRadius()      const {return fZ+fThickness/2;}

protected:

  ClassDef(BeamPipe,1);
};



#endif
