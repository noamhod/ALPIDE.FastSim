#include "KMCUtils.h"
#include <TMath.h>
//---------------------------------
ClassImp(MagField)


MagField::MagField(UInt_t id) {
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;){
    fZMin[i] = 1e9;
    fZMax[i] = 1e-9;
	 
    fXMin[i] = 1e9;
    fXMax[i] = 1e-9;
	 
    fYMin[i] = 1e9;
    fYMax[i] = 1e-9;
	 
    for (int j=3;j--;) fBVal[i][j] = 0;
    for (int j=3;j--;) fBValNonUniform[i][j] = NULL;
    double Ly  = 59.92/10; // cm
    double Lz  = 1440/10; // cm
    double Lx  = (980-2*(202+70))/10; // cm
    double zDipoleCenter = 2050/10; // cm
    double B0y = -0.95*10; // kG
    
    double xmin = -Lx/2.;
    double xmax = +Lx/2.;
    double ymin = -Ly/2.;
    double ymax = +Ly/2.;
    double zmin = zDipoleCenter-Lz/2.;
    double zmax = zDipoleCenter+Lz/2.;
    if(i==0){
        fBValNonUniform[i][1] = new TF3("B1y","[0]*1/((1+exp(([1]-(z - [9]))/[2]))*(1+exp(((z - [9])-[3])/[4])))*1/((1+exp(([5]-x)/[6]))*(1+exp((x-[7])/[8])))",xmin,xmax,ymin,ymax,zmin,zmax);
        fBValNonUniform[i][1]->SetParameter(0,B0y);
        fBValNonUniform[i][1]->SetParameter(1,-616.0/10);
        fBValNonUniform[i][1]->SetParameter(2,28.66/10);
        fBValNonUniform[i][1]->SetParameter(3,622.0/10);
        fBValNonUniform[i][1]->SetParameter(4,28.91/10);
        fBValNonUniform[i][1]->SetParameter(5,-165.0/10);
        fBValNonUniform[i][1]->SetParameter(6,7.7/10);
        fBValNonUniform[i][1]->SetParameter(7,165.0/10);
        fBValNonUniform[i][1]->SetParameter(8,7.7/10);
        fBValNonUniform[i][1]->SetParameter(9,zDipoleCenter);
    }
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
    if(okx and oky and okz){
      // for (int i=3;i--;) bxyz[i] = fBVal[ireg][i];
      for (int i=3;i--;){
        if(fBValNonUniform[ireg][i]==NULL) bxyz[i] = 0;
        else bxyz[i] = fBValNonUniform[ireg][i]->Eval(xyz[0],xyz[1],xyz[2]);
//         if(xyz[2]>(2050/10.0 - 1440/20.0) && xyz[2]<(2050.0/10.0+1440.0/20.0))
//             std::cout << "field " << i << " at xyz=(" << xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ") is " << bxyz[i] << std::endl;
      }
      return;
    }
  }
  return;
}
