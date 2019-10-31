#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <TGeoGlobalMagField.h>

//==========================================================================
class MagField: public TVirtualMagField
{
 public:
  enum {kMaxReg=5}; // max number of field regions may be changed
  MagField(UInt_t id);
  virtual ~MagField() {}
  virtual void Field(const Double_t *xyz, Double_t *bxyz);
  //
  void SetNReg(int n) {fNReg = n;}
  int GetNReg() const {return fNReg;}
  const double* GetZMin() const {return fZMin;}
  const double* GetZMax() const {return fZMax;}
  const double* GetBVals(int ir) const {return &fBVal[ir][0];} 
  void SetZMin(int nreg, double zmin) { fZMin[nreg] = zmin; }
  void SetZMax(int nreg, double zmin) { fZMax[nreg] = zmin; }
  void SetBVals(int nreg, int index, double val) { fBVal[nreg][index] = val; }

 protected:
  int fNReg = 0;
  double fZMin[kMaxReg]; // min z of each field region
  double fZMax[kMaxReg]; // max z of each field region
  double fBVal[kMaxReg][3]; // field values
  //
  ClassDef(MagField, 1) // custom magfield
};


#endif
