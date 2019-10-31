#ifndef KMCCLUSTERFWD_H
#define KMCCLUSTERFWD_H

#include <TObject.h>

class KMCClusterFwd : public TObject {
 public:
  //
  enum {kBitKilled=BIT(14)};
  //
 KMCClusterFwd(Float_t x=0, Float_t y=0, Float_t z=0, int id=-1) : fX(x),fY(y),fZ(z) {SetTrID(id);}
  KMCClusterFwd(KMCClusterFwd &src);
  KMCClusterFwd& operator=(const KMCClusterFwd& src);
  virtual ~KMCClusterFwd() {}
  void  Reset()               {Clear(); Kill();}
  Int_t GetTrID()    const    {return int(GetUniqueID())-1;}
  void  SetTrID(int id)       {SetUniqueID(id+1);}
  //
  Double_t GetX()    const   {return fX;}
  Double_t GetY()    const   {return fY;}
  Double_t GetZ()    const   {return fZ;}
  Double_t GetXLab() const   {return -fZ;}
  Double_t GetYLab() const   {return fY;}
  Double_t GetZLab() const   {return fX;}
  void     SetX(double v)    {fX = v;}
  void     SetY(double v)    {fY = v;}
  void     SetZ(double v)    {fZ = v;}
  //
  void    Kill(Bool_t v=kTRUE)          {SetBit(kBitKilled,v);}
  Bool_t  IsKilled()              const {return TestBit(kBitKilled);}
  void    Set(Float_t x, Float_t y, Float_t z, int id=-1) {fX=x; fY=y; fZ=z; SetTrID(id); ResetBit(kBitKilled);}
  virtual Bool_t  IsSortable() const {return kTRUE;}
  virtual Bool_t  IsEqual(const TObject* obj) const;
  virtual Int_t   Compare(const TObject* obj) const;
  virtual void Print(Option_t *opt = 0) const;
  //
protected:
  Float_t fX; // tracking coordinates
  Float_t fY; 
  Float_t fZ; 
  //
  ClassDef(KMCClusterFwd,1);
};

#endif
