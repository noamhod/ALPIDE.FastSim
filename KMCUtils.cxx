#include "KMCUtils.h"
#include <TMath.h>
//---------------------------------
ClassImp(MagField)


MagField::MagField(UInt_t id) {
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;) {
    fZMin[i] = 1e9;
    fZMax[1] = 1e-9;
    for (int j=3;j--;) fBVal[i][j] = 0;
  }
}

//__________________________________________________
void MagField::Field(const Double_t *xyz, Double_t *bxyz) 
{
  bxyz[0]=bxyz[1]=bxyz[2]=0.;
  for (int ireg=0;ireg<fNReg;ireg++) {
    if((xyz[2] >fZMin[ireg]) && (xyz[2] < fZMax[ireg])) {
      for (int i=3;i--;) bxyz[i] = fBVal[ireg][i];
      return;
    }
  }
  return;
}
