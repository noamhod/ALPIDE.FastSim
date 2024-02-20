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
  }
}

/// constructor for the non-uniform magnetic field
MagField::MagField(UInt_t id, UInt_t nreg, std::vector<double> dipoleConst,
						 std::vector<double> xmin, std::vector<double> xmax,
						 std::vector<double> ymin, std::vector<double> ymax,
						 std::vector<double> zmin, std::vector<double> zmax,
						 std::vector<std::vector<std::string>> functionForm)
{
  SetUniqueID(id);
  fNReg = 0;
  for(unsigned int r=0 ; r<nreg ; r++)
  {
	 std::string reg = std::to_string(r);
    for(unsigned int x=0 ; x<3 ; x++)
	 {
      std::string axs = std::to_string(x);
      if(functionForm[r][x]!="NONE")
		{
      	fBValNonUniform[r][x] = new TF3(("B1_RegId_"+reg+"_dim_"+axs).c_str(),(functionForm[r][x]).c_str(),xmin[r],xmax[r],ymin[r],ymax[r],zmin[r],zmax[r]);
			fBValNonUniform[r][x]->SetParameter(0,dipoleConst[r]);
		}
		else
		{
      	fBValNonUniform[r][x] = new TF3(("B1_RegId_"+reg+"_dim_"+axs).c_str(),"0",xmin[r],xmax[r],ymin[r],ymax[r],zmin[r],zmax[r]);
			fBValNonUniform[r][x]->SetParameter(0,0);
		}
		SetBVals(r,x,dipoleConst[r]);
		SetXMin(r,xmin[r]);
		SetXMax(r,xmax[r]);
		SetYMin(r,ymin[r]);
		SetYMax(r,ymax[r]);
		SetZMin(r,zmin[r]);
		SetZMax(r,zmax[r]);
    }
  }
}



/// constructor for the non-uniform magnetic field
MagField::MagField(UInt_t id, double dipoleConst, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, std::string *functionForm){
  SetUniqueID(id);
  fNReg = 0;
  for (int i=kMaxReg;i--;){
    for (int j=3;j--;){
      std::string reg = std::to_string(i);
      std::string dim = std::to_string(j);
      if(*(functionForm+j)!="NONE"){
        fBValNonUniform[i][j] = new TF3(("B1_RegId_"+reg+"_dim_"+dim).c_str(),(*(functionForm+j)).c_str(),xmin,xmax,ymin,ymax,zmin,zmax);
        fBValNonUniform[i][j]->SetParameter(0,dipoleConst);
      }
      else fBValNonUniform[i][j]=NULL;
    }
  }
}

void MagField::SetFunctionForm(int nreg, int index, std::string funcForm, double dipoleConst, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
   std::string reg = std::to_string(nreg);
   std::string dim = std::to_string(index);
	functionFormStr[nreg][index] = funcForm;
	if(fBValNonUniform[nreg][index]) delete fBValNonUniform[nreg][index];
   fBValNonUniform[nreg][index] = new TF3(("B1_RegId_"+reg+"_dim_"+dim).c_str(),(funcForm).c_str(),xmin,xmax,ymin,ymax,zmin,zmax);
   fBValNonUniform[nreg][index]->SetParameter(0,dipoleConst);
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
