#include "KMCUtils.h"
#include <TMath.h>
//---------------------------------
ClassImp(MagField)

/// constructor for uniform magnetic field
MagField::MagField(UInt_t id) {
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;){
    for (int j=3;j--;) fBVal[i][j] = 0;
    for (int j=3;j--;) fBValNonUniform[i][j] = NULL;
  }
}

/// constructor for the non-uniform magnetic field
MagField::MagField(UInt_t id, double dipoleConst, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, std::string *functionForm){
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;){
    for (int j=3;j--;) fBVal[i][j] = 0;
    for (int j=3;j--;) fBValNonUniform[i][j] = NULL;
    std::string str = std::to_string(i);
    fBValNonUniform[i][1] = new TF3(("B1_RegId_"+str+"_yComp").c_str(),(*(functionForm+1)).c_str(),xmin,xmax,ymin,ymax,zmin,zmax); //// x component uses 0, y component uses 1 and z component uses 2
    fBValNonUniform[i][1]->SetParameter(0,dipoleConst);
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

    if(okx && oky && okz){
      for (int i=3;i--;)
		  {
        if(fBValNonUniform[ireg][i]==NULL) bxyz[i] = 0;
        else
        {
          bxyz[i] = fBValNonUniform[ireg][i]->Eval(xyz[0],xyz[1],xyz[2]);
        }
		  }
      return;
    }
  }
  return;
}
