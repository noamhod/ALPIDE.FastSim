#ifndef KMCUTILS_H
#define KMCUTILS_H

#include <TGeoGlobalMagField.h>
#include <TF3.h>

//==========================================================================
class MagField: public TVirtualMagField
{
 public:
  enum {kMaxReg=5}; // max number of field regions may be changed
  MagField(UInt_t id);
  virtual ~MagField() {}
  virtual void Field(const Double_t *xyz, Double_t *bxyz);
  
  // double fField(const Double_t *xyz);
 
  //
  void SetNReg(int n) {fNReg = n;}
  int GetNReg() const {return fNReg;}
  const double* GetZMin() const {return fZMin;}
  const double* GetZMax() const {return fZMax;}
  const double* GetYMin() const {return fYMin;}
  const double* GetYMax() const {return fYMax;}
  const double* GetXMin() const {return fXMin;}
  const double* GetXMax() const {return fXMax;}
  const double* GetBVals(int ir) const {return &fBVal[ir][0];} 
  void SetZMin(int nreg, double zmin) { fZMin[nreg] = zmin; }
  void SetZMax(int nreg, double zmax) { fZMax[nreg] = zmax; }
  void SetYMin(int nreg, double ymin) { fYMin[nreg] = ymin; }
  void SetYMax(int nreg, double ymax) { fYMax[nreg] = ymax; }
  void SetXMin(int nreg, double xmin) { fXMin[nreg] = xmin; }
  void SetXMax(int nreg, double xmax) { fXMax[nreg] = xmax; }
  void SetBVals(int nreg, int index, double val) { fBVal[nreg][index] = val; }

 protected:
  int fNReg = 0;
  double fZMin[kMaxReg]; // min z of each field region
  double fZMax[kMaxReg]; // max z of each field region
  double fYMin[kMaxReg]; // min y of each field region
  double fYMax[kMaxReg]; // max y of each field region
  double fXMin[kMaxReg]; // min x of each field region
  double fXMax[kMaxReg]; // max x of each field region
  double fBVal[kMaxReg][3]; // field values
  TF3* fBValNonUniform[kMaxReg][3]; // non uniform field values
  //
  ClassDef(MagField, 1) // custom magfield
};


#endif
