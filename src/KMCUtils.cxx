#include "KMCUtils.h"
#include <TMath.h>
//---------------------------------
ClassImp(MagField)


MagField::MagField(UInt_t id) {
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;) {
    fZMin[i] = 1e9;
    // fZMax[1] = 1e-9;
    fZMax[i] = 1e-9;
	 
    fXMin[i] = 1e9;
    fXMax[i] = 1e-9;
	 
    fYMin[i] = 1e9;
    fYMax[i] = 1e-9;
	 
    for (int j=3;j--;) fBVal[i][j] = 0;
  }
}

//---------------------------------
void MagField::Field(const Double_t *xyz, Double_t *bxyz)
{
  bxyz[0]=bxyz[1]=bxyz[2]=0.;
  for (int ireg=0;ireg<fNReg;ireg++)
  {
    bool okz = ((xyz[2] >fZMin[ireg]) && (xyz[2] < fZMax[ireg]));
    bool oky = ((xyz[1] >fYMin[ireg]) && (xyz[1] < fYMax[ireg]));
    bool okx = ((xyz[0] >fXMin[ireg]) && (xyz[0] < fXMax[ireg]));
	 if(okx and oky and okz)
	 {
      for (int i=3;i--;) bxyz[i] = fBVal[ireg][i];
      return;
    }
  }
  return;
}
