#if !defined(__CINT__) || defined(__MAKECINT__)
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "KMCLayerFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TView.h"
#include "TPolyLine3D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "GenMUONLMR.h"
#include "TStopwatch.h"
#include "TTree.h"
#include "TTreeStream.h"
#include <iterator> 
#include <map>
#include <sstream>
#endif

using namespace std;

typedef map<TString, TH1D* > TMapTSTH1D;
typedef map<TString, TH2D* > TMapTSTH2D;
typedef map<TString, int >   TMapTSi;
typedef map<int,TString>     TMapiTS;

TMapTSTH1D histos1;
TMapTSTH2D histos2;

TString storage =  gSystem->ExpandPathName("$STORAGEDIR");

KMCDetectorFwd* det = 0;
double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;

int index_offset_bkg = 100000;
int index_offset_sig = 10000000;

vector<TString> sides{"Eside","Pside"};
vector<TString> layersnames;
vector<double>  layersz;
vector<double> zlayer;
TMapiTS layers;
TMapTSi szlayers;
TMapTSi silayers;
TMapiTS islayers;

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 50; //27.12;   // cm
double Rbeampipe = 2.413; // cm
double RoffsetBfield = 5.7; // cm
double xPsideL = -RoffsetBfield-Lstave;
double xPsideR = -RoffsetBfield;       
double xEsideL = +RoffsetBfield;       
double xEsideR = +RoffsetBfield+Lstave;
double yUp = +Hstave/2.;
double yDn = -Hstave/2.;

//// dipole geometry
double xWdipole = 33.0;
double yHdipole = 10.8;
double z1dipole = 100;
double z2dipole = 202.9;
double zDipoleExit = z2dipole;

// dipole field
double B  = 1.0;
double LB = 1.029; // meters

double zEL1I = -999;
double zEL1O = -999;
double zPL1I = -999;
double zPL1O = -999;

double zEL2I = -999;
double zEL2O = -999;
double zPL2I = -999;
double zPL2O = -999;

double zEL3I = -999;
double zEL3O = -999;
double zPL3I = -999;
double zPL3O = -999;

double zEL4I = -999;
double zEL4O = -999;
double zPL4I = -999;
double zPL4O = -999;

int iEL1I = -999;
int iEL1O = -999;
int iPL1I = -999;
int iPL1O = -999;

int iEL2I = -999;
int iEL2O = -999;
int iPL2I = -999;
int iPL2O = -999;

int iEL3I = -999;
int iEL3O = -999;
int iPL3I = -999;
int iPL3O = -999;

int iEL4I = -999;
int iEL4O = -999;
int iPL4I = -999;
int iPL4O = -999;

double xMinEI = -999;
double xMinEO = -999;
double xMinPI = -999;
double xMinPO = -999;

double xMaxEI = -999;
double xMaxEO = -999;
double xMaxPI = -999;
double xMaxPO = -999;

double zLastLayer  = -999;
double zFirstLayer = -999;
TString LastLayer  = "";
TString FirstLayer = "";

void setParametersFromDet(TString side)
{
	cout << "====================================" << endl;
	cout << "============DEFINITIONS=============" << endl;
	cout << "====================================" << endl;
	KMCLayerFwd* layer_outer = (side=="Eside") ? det->GetLayer("EL1O") : det->GetLayer("PL1O");
	KMCLayerFwd* layer_inner = (side=="Eside") ? det->GetLayer("EL1I") : det->GetLayer("PL1I");
	
	Hstave = layer_outer->GetYMax()-layer_outer->GetYMin();
	Lstave = layer_outer->GetXMax()-layer_outer->GetXMin();
	cout << "Hstave=" << Hstave << ", Lstave=" << Lstave << endl;
	
	xMinEI = (side=="Eside") ? layer_inner->GetXMin() : -999;
	xMinEO = (side=="Eside") ? layer_outer->GetXMin() : -999;
	xMinPI = (side=="Pside") ? layer_inner->GetXMin() : -999;
	xMinPO = (side=="Pside") ? layer_outer->GetXMin() : -999;
	cout << "xMinEI=" << xMinEI << ", xMinEO=" << xMinEO << ", xMinPI=" << xMinPI << ", xMinPO=" << xMinPO << endl;
	
	xMaxEI = (side=="Eside") ? layer_inner->GetXMax() : -999;
	xMaxEO = (side=="Eside") ? layer_outer->GetXMax() : -999;
	xMaxPI = (side=="Pside") ? layer_inner->GetXMax() : -999;
	xMaxPO = (side=="Pside") ? layer_outer->GetXMax() : -999;
	cout << "xMaxEI=" << xMaxEI << ", xMaxEO=" << xMaxEO << ", xMaxPI=" << xMaxPI << ", xMaxPO=" << xMaxPO << endl;
	
	xPsideL = (side=="Eside") ? layer_outer->GetXMin() : -999;
	xPsideR = (side=="Eside") ? layer_inner->GetXMax() : -999;
	xEsideL = (side=="Pside") ? layer_inner->GetXMin() : -999;
	xEsideR = (side=="Pside") ? layer_outer->GetXMax() : -999;
	cout << "xPsideL=" << xPsideL << ", xPsideR=" << xPsideR << ", xEsideL=" << xEsideL << ", xEsideR=" << xEsideR << endl;
	
	yUp = +Hstave/2.;
	yDn = -Hstave/2.;
	cout << "yUp=" << yUp << ", yDn=" << yDn << endl;
	
	//// get the Bfield from the setup
	TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
	MagField* fldm = (MagField*) fld;
	const double BfieldKG = fldm->GetBVals(0,1);  // region 0, the B field is only in the y direction (0,B,0), hence the index 1
	const double* BfieldXminObj = fldm->GetXMin(); // region 0
	const double* BfieldXmaxObj = fldm->GetXMax(); // region 0
	const double* BfieldYminObj = fldm->GetYMin(); // region 0
	const double* BfieldYmaxObj = fldm->GetYMax(); // region 0
	const double* BfieldZminObj = fldm->GetZMin(); // region 0
	const double* BfieldZmaxObj = fldm->GetZMax(); // region 0
	const std::string BFunction = fldm->GetFunctionForm(0,1); /// 1 is for y component of the field
	double BfieldValTesla = BfieldKG/10; /// the B field is only in the y direction (0,B,0), hence the index 1
	double BfieldXmin = BfieldXminObj[0]; // region 0
	double BfieldXmax = BfieldXmaxObj[0]; // region 0
	double BfieldYmin = BfieldYminObj[0]; // region 0
	double BfieldYmax = BfieldYmaxObj[0]; // region 0
	double BfieldZmin = BfieldZminObj[0]; // region 0
	double BfieldZmax = BfieldZmaxObj[0]; // region 0
	cout << "BfieldValTesla=" << BfieldValTesla << ", BfieldXmin=" << BfieldXmin << ", BfieldXmax=" << BfieldXmax << ", BfieldYmin=" << BfieldYmin << ", BfieldYmax=" << BfieldYmax << ", BfieldZmin=" << BfieldZmin << ", BfieldZmax=" << BfieldZmax << endl;
	cout << "Bfunction: " << BFunction << std::endl;
	xWdipole = BfieldXmax-BfieldXmin;
	yHdipole = BfieldYmax-BfieldYmin;
	z1dipole = BfieldZmin;
	z2dipole = BfieldZmax;
	zDipoleExit = z2dipole;
	B  = BfieldValTesla;
	LB = z2dipole-z1dipole;
	cout << "xWdipole=" << xWdipole << ", yHdipole=" << yHdipole << ", z1dipole=" << z1dipole << ", z2dipole=" << z2dipole << endl;
	
	zEL1I = (side=="Eside") ? det->GetLayer("EL1I")->GetZ():-999; iEL1I = (side=="Eside") ? det->GetLayer("EL1I")->GetID():-999;
	zEL1O = (side=="Eside") ? det->GetLayer("EL1O")->GetZ():-999; iEL1O = (side=="Eside") ? det->GetLayer("EL1O")->GetID():-999;
	zPL1I = (side=="Pside") ? det->GetLayer("PL1I")->GetZ():-999; iPL1I = (side=="Pside") ? det->GetLayer("PL1I")->GetID():-999;
	zPL1O = (side=="Pside") ? det->GetLayer("PL1O")->GetZ():-999; iPL1O = (side=="Pside") ? det->GetLayer("PL1O")->GetID():-999;
	cout << "zEL1I=" << zEL1I << ", zEL1O=" << zEL1O << ", zPL1I=" << zPL1I << ", zPL1O=" << zPL1O << endl;

	zEL2I = (side=="Eside") ? det->GetLayer("EL2I")->GetZ():-999;  iEL2I = (side=="Eside") ? det->GetLayer("EL2I")->GetID():-999;
	zEL2O = (side=="Eside") ? det->GetLayer("EL2O")->GetZ():-999;  iEL2O = (side=="Eside") ? det->GetLayer("EL2O")->GetID():-999;
	zPL2I = (side=="Pside") ? det->GetLayer("PL2I")->GetZ():-999;  iPL2I = (side=="Pside") ? det->GetLayer("PL2I")->GetID():-999;
	zPL2O = (side=="Pside") ? det->GetLayer("PL2O")->GetZ():-999;  iPL2O = (side=="Pside") ? det->GetLayer("PL2O")->GetID():-999;
	cout << "zEL2I=" << zEL2I << ", zEL2O=" << zEL2O << ", zPL2I=" << zPL2I << ", zPL2O=" << zPL2O << endl;

	zEL3I = (side=="Eside") ? det->GetLayer("EL3I")->GetZ():-999; iEL3I = (side=="Eside") ? det->GetLayer("EL3I")->GetID():-999;
	zEL3O = (side=="Eside") ? det->GetLayer("EL3O")->GetZ():-999; iEL3O = (side=="Eside") ? det->GetLayer("EL3O")->GetID():-999;
	zPL3I = (side=="Pside") ? det->GetLayer("PL3I")->GetZ():-999; iPL3I = (side=="Pside") ? det->GetLayer("PL3I")->GetID():-999;
	zPL3O = (side=="Pside") ? det->GetLayer("PL3O")->GetZ():-999; iPL3O = (side=="Pside") ? det->GetLayer("PL3O")->GetID():-999;
	cout << "zEL3I=" << zEL3I << ", zEL3O=" << zEL3O << ", zPL3I=" << zPL3I << ", zPL3O=" << zPL3O << endl;

	zEL4I = (side=="Eside") ? det->GetLayer("EL4I")->GetZ():-999; iEL4I = (side=="Eside") ? det->GetLayer("EL4I")->GetID():-999;
	zEL4O = (side=="Eside") ? det->GetLayer("EL4O")->GetZ():-999; iEL4O = (side=="Eside") ? det->GetLayer("EL4O")->GetID():-999;
	zPL4I = (side=="Pside") ? det->GetLayer("PL4I")->GetZ():-999; iPL4I = (side=="Pside") ? det->GetLayer("PL4I")->GetID():-999;
	zPL4O = (side=="Pside") ? det->GetLayer("PL4O")->GetZ():-999; iPL4O = (side=="Pside") ? det->GetLayer("PL4O")->GetID():-999;
	cout << "zEL4I=" << zEL4I << ", zEL4O=" << zEL4O << ", zPL4I=" << zPL4I << ", zPL4O=" << zPL4O << endl;
	
	// IP (vertex) --> start of dipol --> end of dipole and then the layers
	if(side=="Eside")
	{
		layersnames = {"EL1I","EL1O","EL2I","EL2O","EL3I","EL3O","EL4I","EL4O"};
		zlayer   = {0,z1dipole,z2dipole,zEL1I,zEL1O,zEL2I,zEL2O,zEL3I,zEL3O,zEL4I,zEL4O};
		layersz  = {zEL1I,zEL1O,zEL2I,zEL2O,zEL3I,zEL3O,zEL4I,zEL4O};
		layers   = {{iEL1I,"EL1I"},{iEL1O,"EL1O"},{iEL2I,"EL2I"},{iEL2O,"EL2O"},{iEL3I,"EL3I"},{iEL3O,"EL3O"},{iEL4I,"EL4I"},{iEL4O,"EL4O"}};
		szlayers = {{"EL1I",zEL1I},{"EL1O",zEL1O},{"EL2I",zEL2I},{"EL2O",zEL2O},{"EL3I",zEL3I},{"EL3O",zEL3O},{"EL4I",zEL4I},{"EL4O",zEL4O}};
		silayers = {{"EL1I",iEL1I},{"EL1O",iEL1O},{"EL2I",iEL2I},{"EL2O",iEL2O},{"EL3I",iEL3I},{"EL3O",iEL3O},{"EL4I",iEL4I}, {"EL4O",iEL4O}};
		islayers = {{iEL1I,"EL1I"},{iEL1O,"EL1O"},{iEL2I,"EL2I"},{iEL2O,"EL2O"},{iEL3I,"EL3I"},{iEL3O,"EL3O"},{iEL4I,"EL4I"},{iEL4O,"EL4O"}};
	}
	if(side=="Pside")
	{
		layersnames = {"PL1I","PL1O","PL2I","PL2O","PL3I","PL3O","PL4I","PL4O"};
		zlayer   = {0,z1dipole,z2dipole,zPL1I,zPL1O,zPL2I,zPL2O,zPL3I,zPL3O,zPL4I,zPL4O};
		layersz  = {zPL1I,zPL1O,zPL2I,zPL2O,zPL3I,zPL3O,zPL4I,zPL4O};
		layers   = {{iPL1I,"PL1I"},{iPL1O,"PL1O"},{iPL2I,"PL2I"},{iPL2O,"PL2O"},{iPL3I,"PL3I"},{iPL3O,"PL3O"},{iPL4I,"PL4I"},{iPL4O,"PL4O"}};
		szlayers = {{"PL1I",zPL1I},{"PL1O",zPL1O},{"PL2I",zPL2I},{"PL2O",zPL2O},{"PL3I",zPL3I},{"PL3O",zPL3O},{"PL4I",zPL4I},{"PL4O",zPL4O}};
		silayers = {{"PL1I",iPL1I},{"PL1O",iPL1O},{"PL2I",iPL2I},{"PL2O",iPL2O},{"PL3I",iPL3I},{"PL3O",iPL3O},{"PL4I",iPL4I}, {"PL4O",iPL4O}};
		islayers = {{iPL1I,"PL1I"},{iPL1O,"PL1O"},{iPL2I,"PL2I"},{iPL2O,"PL2O"},{iPL3I,"PL3I"},{iPL3O,"PL3O"},{iPL4I,"PL4I"},{iPL4O,"PL4O"}};
	}
	
	zLastLayer  = (side=="Eside") ? zEL4I : zPL4I;
	zFirstLayer = (side=="Eside") ? zPL1O : zPL1O;
	LastLayer  = (side=="Eside") ? "EL4I" : "PL4I";
	FirstLayer = (side=="Eside") ? "PL1O" : "PL1O";
	cout << "zLastLayer=" << zLastLayer << " ("<<LastLayer<<"), zFirstLayer=" << zFirstLayer << " ("<<FirstLayer<<")"<< endl;
	
	cout << "====================================" << endl;
	cout << "====================================" << endl;
}


int acceptcls(double x, double y, double z, double step=0.1)
{
	if(x<0 && (abs(z-zEL1I)>step && abs(z-zEL1O)>step) && (abs(z-zEL2I)>step && abs(z-zEL2O)>step) && (abs(z-zEL3I)>step && abs(z-zEL3O)>step) && (abs(z-zEL4I)>step && abs(z-zEL4O)>step)) return 0;
	if(x>0 && (abs(z-zPL1I)>step && abs(z-zPL1O)>step) && (abs(z-zPL2I)>step && abs(z-zPL2O)>step) && (abs(z-zPL3I)>step && abs(z-zPL3O)>step) && (abs(z-zPL4I)>step && abs(z-zPL4O)>step)) return 0;
	if(x<0 && (x<xMinEO || x>xMaxEI)) return 0;
	if(x>0 && (x<xMinPI || x>xMaxPO)) return 0;
   if(y>yUp    || y<yDn)    return 0;
   return 1;
}

Color_t trkcol(double E)
{
   if     (E>=14)           return kBlack;
   else if(E>=12 and E<14.) return kRed;
   else if(E>=10 and E<12.) return 95;
   else if(E>=8. and E<10.) return 91;
   else if(E>=7. and E<8. ) return 80;
   else if(E>=6. and E<7. ) return 71;
   else if(E>=5. and E<6. ) return 65;
   else if(E>=4. and E<5. ) return 60;
   else if(E>=3. and E<4. ) return 53;
   else if(E>=2. and E<3. ) return 51;
   else if(E>=1. and E<2. ) return 223;
   else if(E>=.5 and E<1. ) return 224;
   else                     return kGray;
   return kRed;
}

TLegend* trkcolleg()
{
   TLegend* leg = new TLegend(0.12,0.30,0.50,0.60);
   leg->SetFillStyle(4000); // will be transparent
   leg->SetFillColor(0);
   leg->SetTextFont(42);
   leg->SetBorderSize(0);
   TLine* l1  = new TLine(); l1->SetLineColor(trkcol(14));   leg->AddEntry(l1,  "#it{E}[GeV]#geq14",    "l");
   TLine* l2  = new TLine(); l2->SetLineColor(trkcol(12));   leg->AddEntry(l2,  "12#leq#it{E}[GeV]<14", "l");
   TLine* l3  = new TLine(); l3->SetLineColor(trkcol(10));   leg->AddEntry(l3,  "10#leq#it{E}[GeV]<12", "l");
   TLine* l4  = new TLine(); l4->SetLineColor(trkcol(8));    leg->AddEntry(l4,  "8#leq#it{E}[GeV]<10",  "l");
   TLine* l5  = new TLine(); l5->SetLineColor(trkcol(7));    leg->AddEntry(l5,  "7#leq#it{E}[GeV]<8",   "l");
   TLine* l6  = new TLine(); l6->SetLineColor(trkcol(6));    leg->AddEntry(l6,  "6#leq#it{E}[GeV]<7",   "l");
   TLine* l7  = new TLine(); l7->SetLineColor(trkcol(5));    leg->AddEntry(l7,  "5#leq#it{E}[GeV]<6",   "l");
   TLine* l8  = new TLine(); l8->SetLineColor(trkcol(4));    leg->AddEntry(l8,  "4#leq#it{E}[GeV]<5",   "l");
   TLine* l9  = new TLine(); l9->SetLineColor(trkcol(3));    leg->AddEntry(l9,  "3#leq#it{E}[GeV]<4",   "l");
   TLine* l10 = new TLine(); l10->SetLineColor(trkcol(2));   leg->AddEntry(l10, "2#leq#it{E}[GeV]<3",   "l");
   TLine* l11 = new TLine(); l11->SetLineColor(trkcol(1));   leg->AddEntry(l11, "1#leq#it{E}[GeV]<2",   "l");
   TLine* l12 = new TLine(); l12->SetLineColor(trkcol(0.5)); leg->AddEntry(l12, "0.5#leq#it{E}[GeV]<1", "l");
   TLine* l13 = new TLine(); l13->SetLineColor(trkcol(0.1)); leg->AddEntry(l13, "#it{E}[GeV]<0.5",      "l");
   return leg;
}

TPolyLine3D* TrackLine3d(const KMCProbeFwd* source, Double_t zMax, Double_t step=1, Color_t col=kBlack)
{	 
	 double xyz[3];
	 source->GetXYZ(xyz);
    double zCurr = round(xyz[2]); //source->GetZ();
    int nZ = (zMax - zCurr)/step + 1;
    if(nZ<2)
	 {
       printf("bad limits\n");
       return 0;
    }
    KMCProbeFwd tmp(*source);
    double xp[nZ],yp[nZ],zp[nZ];
    xp[0] = xyz[0]; // x-vertex
    yp[0] = xyz[1]; // y-vertex
    zp[0] = xyz[2]; // z-vertex
	 int nz = 0;
    for(int iz=1 ; iz<nZ ; iz++)
	 {
       if(!det->PropagateToZBxByBz(&tmp, TMath::Min(tmp.GetZ()+step, zMax), step)) break;
       tmp.GetXYZ(xyz);
       xp[iz] = xyz[0];
       yp[iz] = xyz[1];
       zp[iz] = xyz[2];
       nz++;
    }
    TPolyLine3D *polyline = new TPolyLine3D(nz+1);
	 polyline->SetLineColor(col);
    for(int i=0 ; i<nz+1 ; i++) polyline->SetPoint(i,xp[i],yp[i],zp[i]);
    return polyline;
}

bool islayer(double z, int layerindex=-1, double stepsize=1)
{
	if(layerindex>=0)
	{
		double dz = abs(zlayer[layerindex]-z);
		if(dz<stepsize) return true;
	}
	else
	{
		for(int j=0 ; j<(int)zlayer.size() ; ++j)
		{
			double dz = abs(zlayer[j]-z);
			if(dz<stepsize/2.)
			{
				return true;
			}
		}
	}
	return false;
}

TPolyMarker3D* TrackMarker3d(const KMCProbeFwd* source, double zmin, double zmax, double zstep, Color_t col, bool doPrint=false, bool fullrange=false)
{
    KMCProbeFwd tmp(*source);
	 int nZ = (int)(zmax-zmin)/zstep;
    double xp[nZ],yp[nZ],zp[nZ];
    double xyz[3];
    tmp.GetXYZ(xyz);
    xp[0] = xyz[0];
    yp[0] = xyz[1];
    zp[0] = xyz[2];
	 int nz = 0;
    for(int iz=1;iz<nZ;iz++)
	 {
		 if(!det->PropagateToZBxByBz(&tmp, tmp.GetZ()+zstep, zstep)) break;
       tmp.GetXYZ(xyz);
       xp[iz] = xyz[0];
       yp[iz] = xyz[1];
       zp[iz] = xyz[2];
       nz++;
    }
    TPolyMarker3D *polymarker = new TPolyMarker3D(zlayer.size());
	 polymarker->SetMarkerColor(col);
	 int n = 0;
    for(int i=0;i<nz+1;i++) {
       if(!fullrange && !islayer(zp[i],-1,zstep)) continue;
       polymarker->SetPoint(n,xp[i],yp[i],zp[i]);
       n++;
    }
	 
	 if(doPrint)
	 {
		/// in acceptance of layer 4?
		double xr,yr,zr;
		for(int i=polymarker->GetN()-1 ; i>=0 ; --i)
		{
			 polymarker->GetPoint(i,xr,yr,zr);
			 cout << "in TrackMarker3d: r["<<i<<"]={"<<xr<<","<<yr<<","<<zr<<"}" << endl;
		}
	 }
	 
    return polymarker;
}

TPolyLine3D* GetLayer(TString side, TString io, double z, Color_t col)
{
   Int_t n=5;
	
	Double_t xIE[] = { xMinEI,xMinEI,xMaxEI,xMaxEI,xMinEI };
	Double_t xIP[] = { xMinPI,xMinPI,xMaxPI,xMaxPI,xMinPI };
	
	Double_t xOE[] = { xMinEO,xMinEO,xMaxEO,xMaxEO,xMinEO };
	Double_t xOP[] = { xMinPO,xMinPO,xMaxPO,xMaxPO,xMinPO };
	
   Double_t y[] = {yDn,yUp,yUp,yDn,yDn};
   
	Double_t zC[] = {z,z,z,z,z};
	
	TPolyLine3D* polyline = 0;
   if(side=="P" && io=="I") polyline = new TPolyLine3D(n,xIP,y,zC);
   if(side=="E" && io=="I") polyline = new TPolyLine3D(n,xIE,y,zC);
   if(side=="P" && io=="O") polyline = new TPolyLine3D(n,xOP,y,zC);
   if(side=="E" && io=="O") polyline = new TPolyLine3D(n,xOE,y,zC);
   polyline->SetLineColor(col);
   return polyline;
}

TPolyLine3D* GetLayerFront(TString side, TString io, double z, Color_t col)
{
   Int_t n=4;
	
	Double_t xIE[] = { xMinEI,xMinEI,xMaxEI,xMaxEI };
	Double_t xIP[] = { xMinPI,xMinPI,xMaxPI,xMaxPI };
	
	Double_t xOE[] = { xMinEO,xMinEO,xMaxEO,xMaxEO };
	Double_t xOP[] = { xMinPO,xMinPO,xMaxPO,xMaxPO };
	
   Double_t y[] = {yUp,yDn,yDn,yUp};
   
	Double_t zC[] = {z,z,z,z};
   
	TPolyLine3D* polyline = 0;
   if(side=="P" && io=="I") polyline = new TPolyLine3D(n,xIP,y,zC);
   if(side=="E" && io=="I") polyline = new TPolyLine3D(n,xIE,y,zC);
   if(side=="P" && io=="O") polyline = new TPolyLine3D(n,xOP,y,zC);
   if(side=="E" && io=="O") polyline = new TPolyLine3D(n,xOE,y,zC);
   polyline->SetLineColor(col);
   return polyline;
}

TPolyLine3D* GetDipole(Color_t col)
{
   TPolyLine3D* polyline = new TPolyLine3D();
   polyline->SetPoint(0,-xWdipole/2,-yHdipole/2,z1dipole);
   polyline->SetPoint(1,-xWdipole/2,+yHdipole/2,z1dipole);
   polyline->SetPoint(2,+xWdipole/2,+yHdipole/2,z1dipole);
   polyline->SetPoint(3,+xWdipole/2,-yHdipole/2,z1dipole);
   polyline->SetPoint(4,-xWdipole/2,-yHdipole/2,z1dipole);

   polyline->SetPoint(5,-xWdipole/2,-yHdipole/2,z2dipole); // go up

   polyline->SetPoint(6,-xWdipole/2,+yHdipole/2,z2dipole); // move
   polyline->SetPoint(7,-xWdipole/2,+yHdipole/2,z1dipole); // go down
   polyline->SetPoint(8,-xWdipole/2,+yHdipole/2,z2dipole); // up again

   polyline->SetPoint(9,+xWdipole/2,+yHdipole/2,z2dipole); // move
   polyline->SetPoint(10,+xWdipole/2,+yHdipole/2,z1dipole); // go down
   polyline->SetPoint(11,+xWdipole/2,+yHdipole/2,z2dipole); // up again

   polyline->SetPoint(12,+xWdipole/2,-yHdipole/2,z2dipole); // move
   polyline->SetPoint(13,+xWdipole/2,-yHdipole/2,z1dipole); // go down
   polyline->SetPoint(14,+xWdipole/2,-yHdipole/2,z2dipole); // up again

   polyline->SetPoint(15,-xWdipole/2,-yHdipole/2,z2dipole); // move
   polyline->SetPoint(16,-xWdipole/2,-yHdipole/2,z1dipole); // go down
   polyline->SetPoint(17,-xWdipole/2,-yHdipole/2,z2dipole); // up again

   polyline->SetLineColor(col);
   return polyline;
}

TPolyLine3D* GetDipoleFront(Color_t col)
{
   TPolyLine3D* polyline = new TPolyLine3D();
   polyline->SetPoint(0,-xWdipole/2,-yHdipole/2,z1dipole);
   polyline->SetPoint(1,+xWdipole/2,-yHdipole/2,z1dipole);
   polyline->SetPoint(2,+xWdipole/2,-yHdipole/2,z2dipole);
   polyline->SetPoint(3,-xWdipole/2,-yHdipole/2,z2dipole);
   polyline->SetPoint(4,-xWdipole/2,-yHdipole/2,z1dipole);
   polyline->SetLineColor(col);
   return polyline;
}

bool skipglitches(TPolyMarker3D* points)
{
	Double_t x,y,z;
	for(int n=0 ; n<points->GetN() ; ++n)
	{
		points->GetPoint(n,x,y,z);
		if(abs(x)>40 || abs(y)>5) return true;
	}
	return false;
}

void WriteGeometry(vector<TPolyMarker3D*>& polm, vector<TPolyLine3D*>& poll, TString process, vector<int>& inacc, vector<TPolyMarker3D*>& clusters, TString suff="")
{
   TCanvas* cnv_pl3d = new TCanvas("cnv_pl3d"+suff,"",500,500);
   TView* view_pl3d = TView::CreateView(1);
   view_pl3d->SetRange(-60,-20,0, +60,+20,zLastLayer+15);
   view_pl3d->ShowAxis();
   
   TCanvas* cnv_pm3d = new TCanvas("cnv_pm3d"+suff,"",500,500);
   TView* view_pm3d = TView::CreateView(1);
   view_pm3d->SetRange(-60,-20,0, +60,+20,zLastLayer+15);
   view_pm3d->ShowAxis();
	
	vector<TPolyLine3D*> staves;
	vector<TPolyLine3D*> fstaves;
	for(unsigned int l=0 ; l<layersz.size() ; ++l)
	{
		double  z  = layersz[l];
		TString io = (layersnames[l].Contains("I")) ? "I" : "O";
		TString pe = (layersnames[l].Contains("P")) ? "P" : "E";
		staves.push_back( GetLayer(pe,io,z,kGreen+3) );
		fstaves.push_back( GetLayerFront(pe,io,z,kGreen+3) );
	}

   TPolyLine3D* dipole  = GetDipole(kGray);
	TPolyLine3D* fdipole = GetDipoleFront(kGray);
	
   cnv_pl3d->cd();
   dipole->Draw();
	for(unsigned int l=0 ; l<staves.size() ; l++) staves[l]->Draw();
	
   cnv_pm3d->cd();
   dipole->Draw();
	for(unsigned int l=0 ; l<staves.size() ; l++) staves[l]->Draw();
	
   for(int i=0 ; i<(int)poll.size() ; ++i)
	{
		/// check acceptance
		if(!inacc[i]) continue;
		/// check for glitches
		if(skipglitches(polm[i])) continue;
		
		cnv_pl3d->cd();
		poll[i]->Draw();
		clusters[i]->Draw();
		
		cnv_pm3d->cd();
		polm[i]->Draw();
	}
   
   TLegend* leg = trkcolleg();
   cnv_pl3d->cd();
	fdipole->Draw();
	for(unsigned int l=0 ; l<fstaves.size() ; l++) fstaves[l]->Draw();
   leg->Draw("same");
	
   cnv_pm3d->cd();
   leg->Draw("same");
	
   cnv_pl3d->SaveAs(storage+"/output/root/"+process+"_tracks_pl3d"+suff+".root");
   cnv_pl3d->SaveAs(storage+"/output/pdf/"+process+"_tracks_pl3d"+suff+".pdf");
   cnv_pm3d->SaveAs(storage+"/output/root/"+process+"_tracks_pm3d"+suff+".root");
   cnv_pm3d->SaveAs(storage+"/output/pdf/"+process+"_tracks_pm3d"+suff+".pdf");
   
   TFile* flines = new TFile(storage+"/data/root/"+process+"_geometry"+suff+".root","RECREATE");
   flines->cd();
   dipole->Write();
   fdipole->Write();
	for(unsigned int l=0 ; l<staves.size() ; l++)  { staves[l]->Write(); fstaves[l]->Write(); }
   leg->Write();
   flines->Close();
}

bool accepttrk(vector<TVector3>& clusters, bool fullacc, double step=0.1, int nMinLayers=4)
{
	/// in acceptance?
   int acc = 0;
   Double_t x,y,z;
	for(unsigned int j=0 ; j<clusters.size() ; ++j)
	{
		x = clusters[j].X();
		y = clusters[j].Y();
		z = clusters[j].Z();
		acc += acceptcls(x,y,z,step);
	}
   return (fullacc) ? (acc>=nMinLayers) : (acc>0);
}


/// for the output tree
int ngen = 0;
int nslv = 0;
int nacc = 0;
vector<double>            wgt;
vector<float>             xvtx;
vector<float>             yvtx;
vector<float>             zvtx;
vector<TLorentzVector>    trkp4;
vector<int>               crg;
vector<vector<int> >      clusters_id;
vector<vector<int> >      clusters_layerid;
vector<vector<int> >      clusters_type;
vector<TPolyMarker3D*>    clusters_xyz;
vector<vector<TVector3> > clusters_r;
vector<TPolyMarker3D*>    trkpts_fullrange;
vector<TPolyMarker3D*>    trkpts;
vector<TPolyLine3D*>      trklin;
vector<int>               acc;
TFile* fOut = 0;
TTree* tOut = 0;
void SetOutTree(TString fOutName, TString side)
{
   // fOut = new TFile(fOutName,"RECREATE");
	fOut->cd();
   tOut = new TTree("dig_"+side,"dig_"+side);
   tOut->Branch("ngen",         &ngen);
   tOut->Branch("nslv",         &nslv);
   tOut->Branch("nacc",         &nacc);
   tOut->Branch("wgt",          &wgt);
   tOut->Branch("xvtx",         &xvtx);
   tOut->Branch("yvtx",         &yvtx);
   tOut->Branch("zvtx",         &zvtx);
   tOut->Branch("crg",          &crg);
   tOut->Branch("trkp4",        &trkp4);
   tOut->Branch("acc",          &acc);
   tOut->Branch("clusters_id",  &clusters_id);
   tOut->Branch("clusters_layerid", &clusters_layerid);
   tOut->Branch("clusters_type",&clusters_type);
   tOut->Branch("clusters_xyz", &clusters_xyz);
   tOut->Branch("clusters_r",   &clusters_r);
   tOut->Branch("trkpts_fullrange", &trkpts_fullrange);
   tOut->Branch("trkpts",       &trkpts);
   tOut->Branch("trklin",       &trklin);
	
	/// remove previous histos if any!
	histos1.clear();
	histos2.clear();
}
void KillOutTree()
{
   fOut->cd();
   tOut->Write();
	for(TMapTSTH1D::iterator it=histos1.begin() ; it!=histos1.end() ; ++it) it->second->Write();
	for(TMapTSTH2D::iterator it=histos2.begin() ; it!=histos2.end() ; ++it) it->second->Write();
   fOut->Write();
   // fOut->Close();
}
void RenameOutTree(TString fOutName, int nFiles)
{
	TString sf = "";
	if(nFiles<10)                      sf = Form("0000%d", nFiles);
	if(nFiles>=10 && nFiles<100)       sf = Form("000%d", nFiles);
	if(nFiles>=100 && nFiles<1000)     sf = Form("00%d", nFiles);
	if(nFiles>=1000 && nFiles<10000)   sf = Form("0%d", nFiles);
	if(nFiles>=10000 && nFiles<100000) sf = Form("%d", nFiles); // assume no more than 99,999 files...
	TString fOutNameNew = fOutName;
	fOutNameNew = fOutNameNew.ReplaceAll(".root","_"+sf+".root");
	gSystem->Exec("mv -f "+fOutName+" "+fOutNameNew);
}

TString FormatEventID(int evnt)
{
	TString sevnt = "";
	if(evnt<10)                        sevnt = Form("000000%d", evnt);
	if(evnt>=10 && evnt<100)           sevnt = Form("00000%d", evnt);
	if(evnt>=100 && evnt<1000)         sevnt = Form("0000%d", evnt);
	if(evnt>=1000 && evnt<10000)       sevnt = Form("000%d", evnt);
	if(evnt>=10000 && evnt<100000)     sevnt = Form("00%d", evnt);
	if(evnt>=100000 && evnt<1000000)   sevnt = Form("0%d", evnt);
	if(evnt>=1000000 && evnt<10000000) sevnt = Form("%d", evnt); // assume no more than 9,999,999 events...
	return sevnt;
}

int toint(TString str)
{
	stringstream strm;
	int x;
	strm << str;
	strm >> x;
	return x;
}


void AddCluster(int slvidx, int index_offset, TString process, TString LYR, double dxMax=5)
{
	if(!det->GetLayer(LYR)->IsITS()) return;
	
	// get the reconstructed propagated to the vertex
	int            layerid = det->GetLayer(LYR)->GetID();
   KMCClusterFwd* cluster = det->GetLayer(LYR)->GetMCCluster();

	if(!cluster)              return;
	if(cluster->GetZLab()==0) return;
	
	Double_t x,y,z;
	x = cluster->GetXLab();
	y = cluster->GetYLab();
	z = cluster->GetZLab();
	// if(x<0) cout << "AddCluster: xyz=("<<x<<","<<y<<","<<z<<")" << endl;
	TVector3 v( x,y,z );
	
	/////////////////////////
	/// !!! dirty fix !!! /// //TODO
	/////////////////////////
	unsigned int ncls = clusters_r[slvidx].size();
	double xprev = (ncls>0) ? clusters_r[slvidx][ncls-1].X() : x;
	if(abs(x-xprev)>dxMax) return;
	if(x*xprev<0)          return;
	
   clusters_xyz[slvidx]->SetNextPoint(x,y,z);
   clusters_r[slvidx].push_back( v );
	clusters_id[slvidx].push_back( layerid*index_offset+slvidx ); // assuming no chance to have >index_offset tracks (and hence clusters) per tracker arm (and hence per layer)
	clusters_layerid[slvidx].push_back( layerid );
	clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
	
	cluster->Reset();
	// det->GetLayer(LYR)->GetMCCluster()->Kill();
	
	// cout << "LYR=" << LYR << ", layerid=" << layerid << endl;
}

void reset_layers_all()
{
	// cout << "Nlayers=" << det->GetLayers()->GetEntries() << endl;
	for(Int_t l=0 ; l<det->GetLayers()->GetEntries() ; l++)
	{
		det->GetLayer(l)->ResetMCClusters();
		det->GetLayer(l)->ResetBgClusters();
		det->GetLayer(l)->ResetMCTracks();
		det->GetLayer(l)->ResetMC();
		det->GetLayer(l)->Reset();
	}
}

void reset_layers_tracks(Int_t skip=-1)
{
	Int_t l0 = (skip>=0) ? skip : 0;
	for(Int_t l=l0 ; l<det->GetLayers()->GetEntries() ; l++)
	{
		det->GetLayer(l)->ResetMCTracks();
	}
}

int main(int argc, char *argv[])
{	
	int argcounter; 
	printf("Program Name Is: %s",argv[0]); 
	if(argc>=2) 
	{ 
		printf("\nNumber Of Arguments Passed: %d",argc); 
		printf("\n----Following Are The Command Line Arguments Passed----"); 
		for(argcounter=0;argcounter<argc;argcounter++) printf("\nargv[%d]: %s",argcounter,argv[argcounter]);
		printf("\n");
	}
	//// minimum requirements
	if(argc<2) { printf("argc<2, exitting now\n"); exit(-1); }
	//// validate inputs
	if(argc==2 and !((TString)argv[1]).Contains("-proc=")) { printf("argc=2 but cannot parse %s\n",argv[1]); exit(-1); }
	if(argc==3 and !((TString)argv[2]).Contains("-path=")) { printf("argc=4 but cannot parse %s\n",argv[3]); exit(-1); }
	if(argc==4 and !((TString)argv[3]).Contains("-evnt=")) { printf("argc=3 but cannot parse %s\n",argv[2]); exit(-1); }
	if(argc==5 and !((TString)argv[4]).Contains("-seed=")) { printf("argc=5 but cannot parse %s\n",argv[4]); exit(-1); }
	//// assign inputs
	TString process = ((TString)argv[1]).ReplaceAll("-proc=",""); // mandatory
	TString path    = ((TString)argv[2]).ReplaceAll("-path=",""); // mandatory
	int     evnt    = (argc>3) ? toint(((TString)argv[3]).ReplaceAll("-evnt=","")) : -1; // job id [optional]
	int     Seed    = (argc>4) ? toint(((TString)argv[4]).ReplaceAll("-seed=","")) : 12345; // seed [optional]
	//// print assigned inputs
	cout << "process=" << process << endl;
	cout << "path=" << path << endl;
	cout << "evnt=" << evnt << endl;
	cout << "Seed=" << Seed << endl;
		
	/// common stuff
	TString eventid = (evnt<0) ? "" : FormatEventID(evnt);
	TString proc = process;
	proc.ReplaceAll("_bkg","");
	
   int outN = 100;
   int nMaxEventsPerFile = 100;
   int nFiles = 1;
	
	int index_offset = (process.Contains("bkg")) ? index_offset_bkg : index_offset_sig; // assuming no chance to have >index_offset tracks (and hence clusters) per tracker arm (and hence per layer)
	
   // output tree
   gInterpreter->GenerateDictionary("vector<TLorentzVector>",       "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>",       "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",         "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",         "vector");
	gInterpreter->GenerateDictionary("vector<vector<TVector3> >",    "vector");
	gInterpreter->GenerateDictionary("vector<vector<TPolyLine3D> >", "vector");
	gSystem->Exec("mkdir -p "+storage+"/data/root/dig");
	TString isflat = (path.Contains("flat")) ? "_flat" : "";
	TString fOutName = storage+"/data/root/dig/dig_"+process+"_"+eventid+isflat+".root";
	fOut = new TFile(fOutName,"RECREATE");
	
	/// stuff which depend on the side
	for(int s=0 ; s<sides.size() ; s++)
	{
		TString side = sides[s];
	
		/// setup the detector for one side at a time!!!
		TString setup = "../setup/setupLUXE_"+process+"_"+side+".txt";
   	gRandom->SetSeed(Seed);  
   	det = new KMCDetectorFwd();
   	det->ReadSetup(setup,setup);
   	det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VertexTelescope
   	det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); // require hit in every layer
   	det->SetMinMSHits(0); // we don't have muon spectrometer
   	det->SetMinTRHits(0); // we don't have muon trigger stations
   	// max number of seeds on each layer to propagate (per muon track)
   	det->SetMaxSeedToPropagate(3000); // relevant only if backgrount is considered
   	// set chi2 cuts
   	det->SetMaxChi2Cl(10.);  // max track to cluster chi2
   	det->SetMaxChi2NDF(3.5); // max total chi2/ndf
   	det->SetMaxChi2Vtx(20e9);  // fiducial cut on chi2 of convergence to vtx
   	// IMPORTANT FOR NON-UNIFORM FIELDS
   	det->SetDefStepAir(1);
   	det->SetMinP2Propagate(0.3); //NA60+
   	det->SetIncludeVertex(kTRUE); // count vertex as an extra measured point
   	det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
   	det->SetApplyBransonPCorrection(-1); // Branson correction, only relevant for setup with MS
		// for reconstruction:
		det->SetErrorScale(1000.);
   	det->Print();
   	// det->BookControlHistos();
	
		///////////////////////////////
		setParametersFromDet(side); ///
		///////////////////////////////

		/// get the particles from a ttree
		// TFile* fIn = new TFile(storage+"/data/root/raw_"+process+".root","READ");
		TFile* fIn = new TFile(path+"/raw_"+process+".root","READ");
		TTree* tIn = (TTree*)fIn->Get("tt");
		int nev = tIn->GetEntries();
		vector<double>* vx    = 0;
		vector<double>* vy    = 0;
		vector<double>* vz    = 0;
		vector<double>* px    = 0;
		vector<double>* py    = 0;
		vector<double>* pz    = 0;
		vector<double>* E     = 0;
		vector<double>* wgt0  = 0;
		vector<int>*    pdgId = 0;
		tIn->SetBranchAddress("vx",&vx);
		tIn->SetBranchAddress("vy",&vy);
		tIn->SetBranchAddress("vz",&vz);
		tIn->SetBranchAddress("px",&px);
		tIn->SetBranchAddress("py",&py);
		tIn->SetBranchAddress("pz",&pz);
		tIn->SetBranchAddress("E",&E);
		tIn->SetBranchAddress("wgt",&wgt0);
		tIn->SetBranchAddress("pdgId",&pdgId);
 
		// output tree
		SetOutTree(fOutName,side);
	
		TString hname = "";
		// hname = "h2_z_vs_x_"+side; histos2.insert( make_pair(hname, new TH2D(hname,";x [cm];z [cm];Tracks",1000,-100,+100, 2000,0,+400)) );
		// hname = "h2_z_vs_y_"+side; histos2.insert( make_pair(hname, new TH2D(hname,";y [cm];z [cm];Tracks",1000,-100,+100, 2000,0,+400)) );
		hname = "h2_y_vs_x_exit_"+side; histos2.insert( make_pair(hname, new TH2D(hname,";x_{exit} [cm];y_{exit} [cm];Tracks",200,-15,+15, 200,-0.2,+0.2)) );
		
		hname = "h2_y_vs_x_L1I_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{1} [cm];y [cm];Tracks",200,-50,+50, 200,-0.5,+0.5)) );
		hname = "h2_y_vs_x_L1O_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{1} [cm];y [cm];Tracks",200,-50,+50, 200,-0.5,+0.5)) );
		hname = "h2_y_vs_x_L4I_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];y [cm];Tracks",200,-50,+50, 200,-0.5,+0.5)) );
		hname = "h2_y_vs_x_L4O_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];y [cm];Tracks",200,-50,+50, 200,-0.5,+0.5)) );

		hname = "h2_E_vs_x_L1I_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{1} [cm];E [GeV];Tracks",200,-50,+50, 170,0,+17)) );
		hname = "h2_E_vs_x_L1O_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{1} [cm];E [GeV];Tracks",200,-50,+50, 170,0,+17)) );
		hname = "h2_E_vs_x_L4I_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];E [GeV];Tracks",200,-50,+50, 170,0,+17)) );
		hname = "h2_E_vs_x_L4O_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];E [GeV];Tracks",200,-50,+50, 170,0,+17)) );
		
		hname = "h2_dx14_vs_x_L4I_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];|x_{4}-x_{1}| [cm];Tracks",200,-50,+50, 100,0,+6)) );
		hname = "h2_dx14_vs_x_L4O_"+side;  histos2.insert( make_pair(hname, new TH2D(hname,";x_{4} [cm];|x_{4}-x_{1}| [cm];Tracks",200,-50,+50, 100,0,+6)) );


	   /// loop on events
	   // for(int iev=0;iev<nev;iev++)
		bool fullloop = (evnt<0);
	   for(int iev=(fullloop)?0:evnt ; (fullloop)?iev<nev:iev==evnt ; iev++)
	   {
			/// reset the layers
			reset_layers_all();
			reset_layers_tracks();
			
	      //// clear
	      ngen = 0;    
	      nslv = 0;
	      nacc = 0;
	 	   for(int i=0;i<(int)clusters_id.size();++i) clusters_id[i].clear();
	 	   for(int i=0;i<(int)clusters_layerid.size();++i) clusters_layerid[i].clear();
	 	   for(int i=0;i<(int)clusters_type.size();++i) clusters_type[i].clear();
	 	   for(int i=0;i<(int)clusters_xyz.size();++i) delete clusters_xyz[i];
	 	   for(int i=0;i<(int)clusters_r.size();++i) clusters_r[i].clear();
	 	   for(int i=0;i<(int)trkpts_fullrange.size();++i) delete trkpts_fullrange[i];
	 	   for(int i=0;i<(int)trkpts.size();++i) delete trkpts[i];
	 	   for(int i=0;i<(int)trklin.size();++i) delete trklin[i];
	      wgt.clear();
	      xvtx.clear();
	      yvtx.clear();
	      zvtx.clear();
	      trkp4.clear();
	      crg.clear();
			clusters_id.clear();
			clusters_layerid.clear();
			clusters_type.clear();
	      clusters_xyz.clear();
	      clusters_r.clear();
	      trkpts_fullrange.clear();
	      trkpts.clear();
	      trklin.clear();
	      acc.clear();
	      
	 	   //// get the next entry
	 	   tIn->GetEntry(iev);
	      if((iev%outN)==0) printf("Done %d out of %d\n",iev,nev);
	      // int nfakeHits = 0;
	 	   TLorentzVector ptmp;
	 	   int ngenall = pdgId->size();
			
			
	      /// start looping on truth particles
	      for(int igen=0 ; igen<ngenall ; igen++)
	      {
	         if(abs(pdgId->at(igen))!=11)
	         {
				   cout << "illegal pdgId: " << pdgId->at(igen) << endl;
	            break;
	         }
				
				//////////////////////////////////////////////////////
				// saving one tree per side!!! ///////////////////////
				if(pdgId->at(igen)<0 && side=="Eside") continue; /////
				if(pdgId->at(igen)>0 && side=="Pside") continue; /////
				//////////////////////////////////////////////////////
				
	         //// all the rest
	         ngen++;
				vector<int> vtmp;
				vector<TVector3> vtmpTVector3d;
	         int q = (pdgId->at(igen)==11) ? -1 : +1;
				ptmp.SetXYZM(px->at(igen), py->at(igen), pz->at(igen), meGeV);
				
				/// reset the layers
				reset_layers_all();
				reset_layers_tracks();
				
	         // prepare the probe
				double vX = (process.Contains("bkg")) ? vx->at(igen) : 0.;
				double vY = (process.Contains("bkg")) ? vy->at(igen) : 0.;
				double vZ = (process.Contains("bkg")) ? vz->at(igen) : 0.;
	         bool slv = det->SolveSingleTrack(ptmp.Pt(),ptmp.Rapidity(),ptmp.Phi(), meGeV, q, vX,vY,vZ, 0,1,99);
	         if(!slv) continue; // reconstruction failed
	         nslv++;
				
	         // get the truth trajectory in the B-field
	         KMCProbeFwd* trutrk = det->GetProbe();
				// double pxyztmp[3];
				// trutrk->GetPXYZ(pxyztmp);
				
				int slvidx = nslv-1;
				wgt.push_back( wgt0->at(igen) );
				xvtx.push_back( vx->at(igen) );
				yvtx.push_back( vy->at(igen) );
				zvtx.push_back( vz->at(igen) );
				crg.push_back( q );
				trkp4.push_back( ptmp );
				if(!process.Contains("bkg")) trkpts_fullrange.push_back( TrackMarker3d(trutrk,0,zLastLayer+15,0.1,trkcol(ptmp.E()),false,true) );
				trkpts.push_back( TrackMarker3d(trutrk,0,zLastLayer+15,0.1,trkcol(ptmp.E())) );
				trklin.push_back( TrackLine3d(trutrk,zLastLayer+15,1,trkcol(ptmp.E())) );
				
				/// fill some histos
				double zfirstinner = (side=="Eside") ? szlayers["EL1I"] : szlayers["PL1I"];
				double zlastinner  = (side=="Eside") ? szlayers["EL4I"] : szlayers["PL4I"];
				double zfirstouter = (side=="Eside") ? szlayers["EL1O"] : szlayers["PL1O"];
				double zlastouter  = (side=="Eside") ? szlayers["EL4O"] : szlayers["PL4O"];
				double xMinI = (side=="Eside") ? xMinEI : xMinPI;
				double xMaxI = (side=="Eside") ? xMaxEI : xMaxPI;
				double xMinO = (side=="Eside") ? xMinEO : xMinPO;
				double xMaxO = (side=="Eside") ? xMaxEO : xMaxPO;
				double x1I=0;
				double x4I=0;
				double x1O=0;
				double x4O=0;
				for(Int_t n=0 ; n<trkpts_fullrange[slvidx]->GetN() ; ++n)
				{
					Double_t xp, yp, zp;
					trkpts_fullrange[slvidx]->GetPoint(n,xp,yp,zp);
					// histos2["h2_z_vs_x"]->Fill(xp,zp);
					// histos2["h2_z_vs_y"]->Fill(yp,zp);
					bool isInnerX = (xp>xMinI && xp<xMaxI);
					bool isOuterX = (xp>xMinO && xp<xMaxO);
					if(zp==z2dipole)                              histos2["h2_y_vs_x_exit_"+side]->Fill(xp,yp);
					
					if(zp==zfirstinner && isInnerX) histos2["h2_y_vs_x_L1I_"+side]->Fill(xp,yp);
					if(zp==zlastinner  && isInnerX) histos2["h2_y_vs_x_L4I_"+side]->Fill(xp,yp);
					if(zp==zfirstouter && isOuterX) histos2["h2_y_vs_x_L1O_"+side]->Fill(xp,yp);
					if(zp==zlastouter  && isOuterX) histos2["h2_y_vs_x_L4O_"+side]->Fill(xp,yp);
					
					if(zp==zfirstinner && isInnerX) histos2["h2_E_vs_x_L1I_"+side]->Fill(xp,trkp4[slvidx].E());
					if(zp==zlastinner  && isInnerX) histos2["h2_E_vs_x_L4I_"+side]->Fill(xp,trkp4[slvidx].E());
					if(zp==zfirstouter && isOuterX) histos2["h2_E_vs_x_L1O_"+side]->Fill(xp,trkp4[slvidx].E());
					if(zp==zlastouter  && isOuterX) histos2["h2_E_vs_x_L4O_"+side]->Fill(xp,trkp4[slvidx].E());
					
					if(zp==zfirstinner && isInnerX) x1I = xp;
					if(zp==zlastinner  && isInnerX) x4I = xp;
					if(zp==zfirstouter && isOuterX) x1O = xp;
					if(zp==zlastouter  && isOuterX) x4O = xp;
				}
				double dx41I = 0;
				double dx410 = 0;
				if(x4I>xMinI && x4I<xMaxI)
				{
					if((x1I>xMinI && x1I<xMaxI)) { dx41I = abs(x4I-x1I); histos2["h2_dx14_vs_x_L4I_"+side]->Fill(x4I,dx41I); }
				}
				if(x4O>xMinO && x4O<xMaxO)
				{
					if     ((x1O>xMinO && x1O<xMaxO)) { dx410 = abs(x4O-x1O); histos2["h2_dx14_vs_x_L4O_"+side]->Fill(x4O,dx410); }
					else if((x1I>xMinI && x1I<xMaxI)) { dx410 = abs(x4O-x1I); histos2["h2_dx14_vs_x_L4O_"+side]->Fill(x4O,dx410); }
				}
				
	         
				clusters_id.push_back( vtmp );
				clusters_layerid.push_back( vtmp );
				clusters_type.push_back( vtmp );
				clusters_xyz.push_back( new TPolyMarker3D() );
				clusters_r.push_back( vtmpTVector3d );
				acc.push_back( 0 );
				
				// get the reconstructed propagated to the vertex
				for(unsigned int k=0 ; k<layersnames.size() ; k++)
				{
					TString LYR = layersnames[k];
					if(!det->GetLayer(LYR)->IsITS())       continue;
					if(crg[slvidx]<0 && LYR.Contains("P")) continue;
					if(crg[slvidx]>0 && LYR.Contains("E")) continue;
					AddCluster(slvidx,index_offset,process,LYR);
				}
				
				int nclusters = clusters_layerid[slvidx].size(); // same as layers hit by the track
				//cout << "slvidx=" << slvidx << ", E=" << trkp4[slvidx].E() << ", Q=" << crg[slvidx] << " --> Nclusters=" << nclusters << endl;
				for(int j=0 ; j<nclusters ; ++j)
				{
				   Double_t x,y,z;
				   // clusters_xyz[slvidx]->GetPoint(j,x,y,z);
				   x = clusters_r[slvidx][j].X();
				   y = clusters_r[slvidx][j].Y();
				   z = clusters_r[slvidx][j].Z();
					int layerid = clusters_layerid[slvidx][j];
					//cout << "point " << j << ", layerid=" << layerid << ", layername=" << det->GetLayer(layerid)->GetName() << ", xyz=(" << x << "," << y << "," << z << ")" << endl;
					// det->GetLayer(layerid)->Print();
				}
				
				
				/// check acceptance
				acc[slvidx] = (accepttrk(clusters_r[slvidx],false));
				if(acc[slvidx]) nacc++;
	      }
			if(iev==0) WriteGeometry(trkpts,trklin,process,acc,clusters_xyz,"_truth");
			if(iev%1==0) cout << "iev=" << iev << " --> ngen=" << ngen << ", nslv=" << nslv << ", nacc=" << nacc << endl;
	      if(nslv!=ngen and !process.Contains("bkg")) cout << "Warning: nslv=" << nslv << ", ngen=" << ngen << " --> problem" << endl;
			
			/// fill the tree
	      fOut->cd();
	      tOut->Fill();
			
			// int nevprocessed = iev+1; // iev starts from 0
			// if(process.Contains("bkg") and fullloop)
			// {
			// 	if(nevprocessed%nMaxEventsPerFile==0)
			// 	{
			// 		cout << "Killing file #" << nFiles << " and redefining new one" << endl;
			//    	KillOutTree();
			//
			// 		cout << "Renaming file #" << nFiles << " to include its index" << endl;
			// 		RenameOutTree(fOutName,nFiles);
			// 		nFiles++; // must propagate!
			//
			// 		// reset the file and tree as usual
			// 		cout << "Resetting file #" << nFiles << " and tree" << endl;
			// 		SetOutTree(fOutName,side);
			// 	}
			// }
				
	   }
	   KillOutTree();
		if(process.Contains("bkg") and eventid=="") RenameOutTree(fOutName,nFiles);
	}
	return 0;
}
