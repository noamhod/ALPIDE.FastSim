#if !defined(__CINT__) || defined(__MAKECINT__)
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
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
#include "TChain.h"
#include "TTreeStream.h"
#include "TList.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TObject.h"
#include "TStopwatch.h"
#include <iterator> 
#include <map>
#include <sstream>
#endif

using namespace std;

struct Cluster {
	int type;
	int lyrid;
	TString lyrnm;
	int clsid;
	int npixl;
	int shape;
	double charge;
	TVector3 r;
};

typedef map<int, int>                TMapii;
typedef map<TString, int >           TMapTSi;
typedef map<TString, float >         TMapTSf;
typedef map<TString, double >        TMapTSd;
typedef map<TString, vector<int> >   TMapTSvi;
typedef map<TString, vector<float> > TMapTSvf;
typedef map<int,TString>             TMapiTS;
typedef map<TString, TH1D* >         TMapTSTH1D;
typedef map<TString, vector<Cluster> > TMapTSvCls; // formatted as 

TString storage =  gSystem->ExpandPathName("$STORAGEDIR");

int nMinHits = 4;

double vX=0,vY=0,vZ=0; // event vertex
KMCDetectorFwd* det = 0;

double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;
double meGeV2 = meGeV*meGeV;
double cm2m = 0.01;

// for matching these must be the same as in digitisation
int index_offset_bkg = 100000;
int index_offset_sig = 10000000;

//// temp stuff:
double bestdistance = -1;
int bestmatchtrui = -1;
int bestmatchreci = -1;

//// uncertainties
float dxAlignmentXFEL = 0.005; //0.005; // cm
float dyAlignmentXFEL = 0.005; //0.005; // cm
float XvariationSign = +1.;
float YvariationSign = +1.;
float dxAlignmentInTray = 0.000; // cm
float dyAlignmentInTray = 0.000; // cm
bool doMisalignmentX = false;
bool doMisalignmentY = false;

//// seed energies
double EseedMin = 1.0; // GeV
double EseedMaxBPPP = 16.0; // GeV
double EseedMaxTRIDENT = 5.0; // GeV

// double yDipoleExitMin = -0.05; ## cm --> TODO: need tuning
// double yDipoleExitMax = +0.05; ## cm --> TODO: need tuning
// double xAbsMargins = 0.025; # cm --> TODO: need tuning
// double yAbsMargins = 0.025 if(proc=="glaser") else 0.1 # cm --> TODO: need tuning

vector<TString> sides{"Eside","Pside"};
// vector<TString> iolyr{"Inner","Outer"};
// vector<TString> coord{"x","y","z"};
// vector<TString> attri{"id","type"};
vector<TString> layersnames;
vector<double>  layersz;
vector<double> zlayer;
TMapiTS layers;
TMapTSi szlayers;
TMapTSi silayers;
TMapiTS islayers;

TMapTSvCls cached_clusters; /// formatted per side per layer e.g. as: cached_clusters[side+"_"+layerid][i].r.X() or cached_clusters[side+"_"+layerid][i].clsid
// TMapTSvi   cached_clusters_att; /// attributes (type and id)
TMapii     cached_clusters_all_ids; /// 
TMapii     cached_clusters_id2lyr; /// 


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

bool accept(double x, double y)
{
	bool failx = (x<xPsideL || (x>xPsideR && x<xEsideL) || x>xEsideR);
	bool faily = (y>yUp || y<yDn);
	if(failx || faily) return false;
	return true;
}

void SetLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
}

Color_t trkcol(double E)
{
	if     (E>=14)           return kBlack;
	else if(E<14. and E>=12) return kRed;
	else if(E<12. and E>=10) return 95;
	else if(E<10. and E>=8.) return 91;
	else if(E<8.  and E>=7.) return 80;
	else if(E<7.  and E>=6.) return 71;
	else if(E<6.  and E>=5.) return 65;
	else if(E<5.  and E>=4.) return 60;
	else if(E<4.  and E>=3.) return 53;
	else if(E<3.  and E>=2.) return 51;
	else                     return 6;
	return kRed;
}

TLegend* trkcolleg()
{
	TLegend* leg = new TLegend(0.12,0.30,0.50,0.60);
	leg->SetFillStyle(4000); // will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	TLine* l1 = new TLine();  l1->SetLineColor(trkcol(14)); leg->AddEntry(l1,  "#it{E}[GeV]#geq14",    "l");
	TLine* l2 = new TLine();  l2->SetLineColor(trkcol(12)); leg->AddEntry(l2,  "12#leq#it{E}[GeV]<14", "l");
	TLine* l3 = new TLine();  l3->SetLineColor(trkcol(10)); leg->AddEntry(l3,  "10#leq#it{E}[GeV]<12", "l");
	TLine* l4 = new TLine();  l4->SetLineColor(trkcol(8));  leg->AddEntry(l4,  "8#leq#it{E}[GeV]<10",  "l");
	TLine* l5 = new TLine();  l5->SetLineColor(trkcol(7));  leg->AddEntry(l5,  "7#leq#it{E}[GeV]<8",   "l");
	TLine* l6 = new TLine();  l6->SetLineColor(trkcol(6));  leg->AddEntry(l6,  "6#leq#it{E}[GeV]<7",   "l");
	TLine* l7 = new TLine();  l7->SetLineColor(trkcol(5));  leg->AddEntry(l7,  "5#leq#it{E}[GeV]<6",   "l");
	TLine* l8 = new TLine();  l8->SetLineColor(trkcol(4));  leg->AddEntry(l8,  "4#leq#it{E}[GeV]<5",   "l");
	TLine* l9 = new TLine();  l9->SetLineColor(trkcol(3));  leg->AddEntry(l9, "3#leq#it{E}[GeV]<4",   "l");
	TLine* l10 = new TLine(); l10->SetLineColor(trkcol(2)); leg->AddEntry(l10, "2#leq#it{E}[GeV]<3",   "l");
	TLine* l11 = new TLine(); l11->SetLineColor(trkcol(1)); leg->AddEntry(l11, "#it{E}[GeV]<2",        "l");
	return leg;
}

TPolyLine3D* TrackLine3d(const KMCProbeFwd* source, Double_t zMax, Double_t step=1, Color_t col=kBlack)
{
	double xyz[3];
	source->GetXYZ(xyz);
	double zCurr = xyz[2]; //source->GetZ();
	int nZ = (zMax - zCurr)/step + 1;
	if (nZ<2) {
		printf("bad limits\n");
		return 0;
	}
	KMCProbeFwd tmp(*source);
	double xp[nZ],yp[nZ],zp[nZ];
	xp[0] = xyz[0];
	yp[0] = xyz[1];
	zp[0] = xyz[2];
	int nz = 0;
	for (int iz=1;iz<nZ;iz++) {
		if(!det->PropagateToZBxByBz(&tmp, TMath::Min(tmp.GetZ()+step, zMax), step)) break; //propagation may fail..
		tmp.GetXYZ(xyz);
		xp[iz] = xyz[0];
		yp[iz] = xyz[1];
		zp[iz] = xyz[2];
		nz++;
	}
	TPolyLine3D *polyline = new TPolyLine3D(nz+1);
	polyline->SetLineColor(col);
	for (int i=0;i<nz+1;i++) {
		polyline->SetPoint(i,xp[i],yp[i],zp[i]);
	}
	return polyline;
}

// bool islayer(double z)
// {
// 	for(int j=0 ; j<(int)zlayer.size() ; ++j)
// 	{
// 		double dz = abs(zlayer.[j]-z);
// 		if(dz<1.e-6) return true;
// 	}
// 	return false;
// }
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



TPolyMarker3D* TrackMarker3d(const KMCProbeFwd* source, double zmin, double zmax, double zstep, Color_t col=kBlack)
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
	for(int iz=1;iz<nZ;iz++) {
		if(!det->PropagateToZBxByBz(&tmp, tmp.GetZ()+zstep, zstep)) break; //propagation may fail...
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
		if(!islayer(zp[i])) continue;
		polymarker->SetPoint(n,xp[i],yp[i],zp[i]);
		n++;
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


bool accepttrk(vector<TPolyMarker3D*>& polm, int itrk)
{
	/// in acceptance?
	int nlayers = 4;
	int acctrk = 0;
	for (int i=0 ; i<polm[itrk]->GetN() ; i++)
	{
		Double_t xr,yr,zr;
		polm[itrk]->GetPoint(i,xr,yr,zr);
		if(zr<zLastLayer+15) continue; //// count only the active layers
		int inacclayer = accept(xr,yr);
		acctrk += inacclayer;
	}
	return (acctrk==nlayers);
}

bool foundinvec(int x, vector<int>& v)
{
	vector<int>::iterator it = find(v.begin(),v.end(),x);
	return (it!=v.end());
}

int getvecindex(int x, vector<int>& v)
{
	vector<int>::iterator it = find(v.begin(),v.end(),x);
	return (it!=v.end()) ? distance(v.begin(),it) : -1;
}

int CheckMatchingByID(int id1, int id2, int id3, int id4)
{	
	if(id1<index_offset_sig) return -1;
	int nNeg = (id1<0)+(id2<0)+(id3<0)+(id4<0);
 	if(nNeg>1) {cout << "too many negative clusters: " << id1 << "," << id2 << "," << id3 << "," << id4 << endl; return -1;}
	int i1 = id1-1*index_offset_sig;
	int i2 = id2-2*index_offset_sig;
	int i3 = id3-3*index_offset_sig;
	int i4 = id4-4*index_offset_sig;
	vector<int> posids;
	vector<int> posixs;
	if(id1>=0) { posids.push_back(id1); posixs.push_back(i1); }
	if(id2>=0) { posids.push_back(id2); posixs.push_back(i2); }
	if(id3>=0) { posids.push_back(id3); posixs.push_back(i3); }
	if(id4>=0) { posids.push_back(id4); posixs.push_back(i4); }
	if(posids.size()==3 and (posixs[0]!=posixs[1] or posixs[0]!=posixs[2] or posixs[1]!=posixs[2])) {cout << "from 3, idi!=idj: " << id1 << "," << id2 << "," << id3 << "," << id4 << endl; return -1;}
	if(posids.size()==4 and (posixs[0]!=posixs[1] or posixs[0]!=posixs[2] or posixs[0]!=posixs[3] or posixs[1]!=posixs[2] or posixs[1]!=posixs[3] or posixs[2]!=posixs[3])) {cout << "from 4, idi!=idj: " << id1 << "," << id2 << "," << id3 << "," << id4 << endl; return -1;}
	return posixs[0];
}


void prepare_cahced_clusters()
{
	/// clear first
	for(TMapTSvCls::iterator it=cached_clusters.begin() ; it!=cached_clusters.end() ; ++it) it->second.clear();
	cached_clusters.clear();
	/// then rebuild
	vector<Cluster> vc;
	vector<int> vi;
	for(TMapiTS::iterator it=layers.begin() ; it!=layers.begin() ; ++it)
	{
		// int     layerid   = it->first;
		TString layername = it->second;
		cached_clusters.insert( make_pair(layername,vc) );
	}
}


void clear_cached_clusters()
{
	cached_clusters_id2lyr.clear();
	for(TMapTSvCls::iterator it=cached_clusters.begin() ; it!=cached_clusters.end() ; ++it) it->second.clear();
}


int cache_clusters(vector<vector<TVector3> >* clusters_r, vector<vector<int> >* clusters_type, vector<vector<int> >* clusters_id, vector<vector<int> >* clusters_layerid, TString side, int nMaxToCache=1410065407)
{
	int ncached = 0;
	int ntrks = (int)clusters_r->size();
	int ntrkmax = (nMaxToCache>0 && nMaxToCache<ntrks) ? nMaxToCache : ntrks;
	for(int i=0 ; i<ntrkmax ; i++)
	{
		// if(acc) // if the vector is provided (for background only)
		// {
		// 	if(!acc->at(i)) continue; // check acceptnce
		// }
		for(Int_t j=0 ; j<clusters_r->at(i).size() ; ++j)
		{
			double x = clusters_r->at(i)[j].X();
			double y = clusters_r->at(i)[j].Y();
			double z = clusters_r->at(i)[j].Z();
			if(x<0 and side=="Pside") continue;
			if(x>0 and side=="Eside") continue;
			int lyrid = clusters_layerid->at(i)[j];
			TString lyrname = layers[lyrid];
			int clstype = clusters_type->at(i)[j];
			int clsid = clusters_id->at(i)[j];
			
			x = (doMisalignmentX) ? x+XvariationSign*dxAlignmentXFEL : x;
			y = (doMisalignmentY) ? y+YvariationSign*dyAlignmentXFEL : y;
			
			Cluster cls;
			cls.type  = clstype;
			cls.lyrid = lyrid;
			cls.lyrnm = lyrname;
			cls.clsid = clsid;
			cls.npixl = -1;
			cls.shape = -1;
			cls.charge = -1;
			cls.r.SetXYZ(x,y,z);
			
			cached_clusters[lyrname].push_back(cls);
			cached_clusters_id2lyr.insert(make_pair(clsid,lyrid));
			ncached++;
		}
	}
	return ncached;
}


void reset_layers_all()
{
	for(Int_t l=0 ; l<det->GetLayers()->GetEntries() ; l++)
	{
		det->GetLayer(l)->ResetBgClusters();
		det->GetLayer(l)->ResetMCTracks();
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

void embed_cluster(int iLayer, float x, float y, float z, int id)
{
	/// set the clusters of the seed
	double clxyzTrk[3];
	double clxyzLab[3];
	clxyzLab[0]=x;
	clxyzLab[1]=y;
	clxyzLab[2]=z;
	KMCProbeFwd::Lab2Trk(clxyzLab, clxyzTrk);
	det->GetLayer(iLayer)->AddBgCluster(clxyzTrk[0], clxyzTrk[1], clxyzTrk[2], id);
}


void add_all_clusters(TString side)
{
	for(TMapiTS::iterator it=layers.begin() ; it!=layers.end() ; ++it)
	{
		int     lid = it->first;
		TString slr = it->second;
		for(unsigned int j=0 ; j<cached_clusters[slr].size() ; ++j)
		{
			float x = cached_clusters[slr][j].r.X();
			float y = cached_clusters[slr][j].r.Y();
			float z = cached_clusters[slr][j].r.Z();
			int cid = cached_clusters[slr][j].clsid;
			embed_cluster(lid,x,y,z,cid);	
		}
	}
	/// must sort clusters
	for(TMapiTS::iterator it=layers.begin() ; it!=layers.end() ; ++it)
	{
		int lid = it->first;
		det->GetLayer(lid)->GetMCCluster()->Kill();
		det->GetLayer(lid)->SortBGClusters(); /// sort!!!
		/// after sorting, need to map the cluster ids to their indices!!!
		for(int n=0 ; n<det->GetLayer(lid)->GetNBgClusters() ; ++n)
		{
			int id = det->GetLayer(lid)->GetBgCluster(n)->GetTrID();
			cached_clusters_all_ids.insert( make_pair(id,n) );
		}
	}
}


void print_all_clusters(TString side, bool doprint = true)
{	
	if(!doprint) return;
	for(TMapiTS::iterator it=layers.begin() ; it!=layers.end() ; ++it)
	{
		int     ilr = it->first;
		TString slr = it->second;
		for(int c=0 ; c<det->GetLayer(ilr)->GetNBgClusters() ; c++)
		{
			KMCClusterFwd* cluster = det->GetLayer(ilr)->GetBgCluster(c);
			int  id = cluster->GetTrID();
			float x = cluster->GetXLab();
			float y = cluster->GetYLab();
			float z = cluster->GetZLab();
			cout << "side=" << side << ", layer=" << ilr << ", id=" << id << " --> r={" << x << ", " << y << ", " << z << "}" << endl;
		}
		cout << endl;
	}
	for(TMapii::iterator it=cached_clusters_all_ids.begin() ; it!=cached_clusters_all_ids.end() ; it++)
	{
		cout << "id=" << it->first << " --> index=" << it->second << endl;
	}
}

int fill_output_clusters(TString side, vector<vector<TVector3> >& r, vector<int>& ctype, vector<int>& cid)
{
	int nclusters = 0;
	for(TMapiTS::iterator it=layers.begin() ; it!=layers.end() ; ++it)
	{
		vector<TVector3> v3tmp;
		r.push_back(v3tmp);
		int     ilr = it->first;
		TString slr = it->second;
		for(int c=0 ; c<det->GetLayer(ilr)->GetNBgClusters() ; c++)
		{
			KMCClusterFwd* cluster = det->GetLayer(ilr)->GetBgCluster(c);
			int  id  = cluster->GetTrID();
			float x  = cluster->GetXLab();
			float y  = cluster->GetYLab();
			float z  = cluster->GetZLab();
	
			/// TODO!!!
			// int idx = getvecindex(id, cached_clusters_att["id_"+islayers[ilr]+"_"+side]);
			// int typ = (idx>=0) ? cached_clusters_att["type_"+islayers[l]+"_"+side][ idx ] : -3;
			// cout << "id=" << id << ", idx=" << idx << ", typ=" << typ << endl;
			// if(idx<0) cout << "WARNING: cannot find in=" << id << " in cached clusters vector" << endl;
			
			int thisindex = r.size()-1;
			
			TVector3 point(x,y,z);
			r[thisindex].push_back(point);
			// ctype.push_back(typ);
			ctype.push_back(-999); //TODO!!!
			// cid.push_back(id);
			cid.push_back(-999); //TODO!!!
			nclusters++;
		}
	}
	return nclusters;
}

TVector2 rUnit2(TVector2 r1, TVector2 r2)
{
	TVector2 r = (r2-r1).Unit();
	return r;
}

float xofz(float* r1,float* r2, float z)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dz==0)
	{
		cout << "ERROR in xofz: dz=0" << endl;
		exit(-1);
	}
	float a = dx/dz;
	float b = r1[0]-a*r1[2];
	float x = a*z+b;
	// cout << "in xofz: x=" << x << ", dz=" << dz << endl;
	return x;
}
	
float yofz(float* r1,float* r2, float z)
{
	float dz = r2[2]-r1[2];
	float dy = r2[1]-r1[1];
	if(dz==0)
	{
		cout << "ERROR in yofz: dz=0" << endl;
		exit(-1);
	}
	float a = dy/dz;
	float b = r1[1]-a*r1[2];
	float y = a*z+b;
	return y;
}

float zofx(float* r1,float* r2, float x)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dx==0)
	{
		cout << "ERROR in zofx: dx=0" << endl;
		exit(-1);
	}
	float a = dz/dx;
	float b = r1[2]-a*r1[0];
	float z = a*x+b;
	return z;
}


///TODO!!!
bool check_clusters(unsigned int i1, unsigned int i4, TString side)
{
	TString slyr1 = (side=="Eside") ? "EL1I" : "PL1I";
	TString slyr2 = (side=="Eside") ? "EL2I" : "PL2I";
	TString slyr3 = (side=="Eside") ? "EL3I" : "PL3I";
	TString slyr4 = (side=="Eside") ? "EL4I" : "PL4I";
	
	float yAbsMargins = 0.03; //0.02; // cm (a "road" of 200 microns around the line between r4 and r1)
	float xAbsMargins = 0.03; //0.02; // cm (a "road" of 200 microns around the line between r4 and r1)
	float r1min[3]; r1min[0]=cached_clusters[slyr1][i1].r.X()-xAbsMargins; r1min[1]=cached_clusters[slyr1][i1].r.Y()-yAbsMargins; r1min[2]=cached_clusters[slyr1][i1].r.Z();
	float r1max[3]; r1max[0]=cached_clusters[slyr1][i1].r.X()+xAbsMargins; r1max[1]=cached_clusters[slyr1][i1].r.Y()+yAbsMargins; r1max[2]=cached_clusters[slyr1][i1].r.Z();
	float r4min[3]; r4min[0]=cached_clusters[slyr4][i4].r.X()-xAbsMargins; r4min[1]=cached_clusters[slyr4][i4].r.Y()-yAbsMargins; r4min[2]=cached_clusters[slyr4][i4].r.Z();
	float r4max[3]; r4max[0]=cached_clusters[slyr4][i4].r.X()+xAbsMargins; r4max[1]=cached_clusters[slyr4][i4].r.Y()+yAbsMargins; r4max[2]=cached_clusters[slyr4][i4].r.Z();

	/// check possible clusters in layer 2
	float y2min = yofz(r1min,r4min,(float)szlayers[slyr2]);
	float y2max = yofz(r1max,r4max,(float)szlayers[slyr2]);
	float x2min = xofz(r1min,r4min,(float)szlayers[slyr2]);
	float x2max = xofz(r1max,r4max,(float)szlayers[slyr2]);
	bool accept2 = false;
	for(unsigned int i2=0 ; i2<cached_clusters[slyr2].size() ; ++i2)
	{
		bool acceptyz = ( cached_clusters[slyr2][i2].r.Y()>=y2min and cached_clusters[slyr2][i2].r.Y()<=y2max );
		if(!acceptyz) continue;
		bool acceptxz = ( cached_clusters[slyr2][i2].r.X()>=x2min and cached_clusters[slyr2][i2].r.X()<=x2max );
		if(!acceptxz) continue;
		accept2 = true;
		break;
	}
	if(!accept2) return false;

	float y3min = yofz(r1min,r4min,(float)szlayers[slyr3]);
	float y3max = yofz(r1max,r4max,(float)szlayers[slyr3]);
	float x3min = xofz(r1min,r4min,(float)szlayers[slyr3]);
	float x3max = xofz(r1max,r4max,(float)szlayers[slyr3]);
	bool accept3 = false;
	for(unsigned int i3=0 ; i3<cached_clusters[slyr3].size() ; ++i3)
	{
		bool acceptyz = ( cached_clusters[slyr3][i3].r.Y()>=y3min and cached_clusters[slyr3][i3].r.Y()<=y3max );
		if(!acceptyz) continue;
		bool acceptxz = ( cached_clusters[slyr3][i3].r.X()>=x3min and cached_clusters[slyr3][i3].r.X()<=x3max );
		if(!acceptxz) continue;
		accept3 = true;
		break;
	}
	if(!accept3) return false;

	return true;
}

bool adaptive_dx14_vs_x4(double x4, double x1, TF1* fDxvsX, double width0, int tirals=3)
{
	double dx = abs(x4-x1);
	if(dx>5) return false;
	double width = width0;
	double step  = width0/4.;
	for(int i=1 ; i<=tirals ; ++i)
	{
		bool passUp = (dx<(fDxvsX->Eval(x4)+width));
		bool passDn = (dx>(fDxvsX->Eval(x4)-width));
		// cout << i << ": width=" << width << ", passUp=" << passUp << ", passDn=" << passDn << endl;
		if(passUp && passDn) return true;
		width = width0 + i*step; // increase the width by 20% each time
	}
	return false;
}



bool makeseed_nonuniformB(TString process, float* r1, float* r4, unsigned int i1, unsigned int i4, TString side, TLorentzVector& p, TF1* fEvsX, TF1* fDxvsX)
{
	if(abs(r1[0])>=abs(r4[0]))            return false; // |x1| must be smaller than |x4|
	if(r1[0]*r4[0]<0)                     return false; // not on the same side...
	if(r1[2]==r4[2])                      return false; // trivial, make sure z is not the same
	float yDipoleExitAbsMax = (process=="glaser") ? 0.12 : 0.7; // cm
	float xDipoleExitAbsMin = (process=="glaser") ? 2.  : 4. ; // cm
	float xDipoleExitAbsMax = (process=="glaser") ? 15. : 30.; // cm
	float yDipoleExit = yofz(r1,r4,zDipoleExit);
	float xDipoleExit = xofz(r1,r4,zDipoleExit);
	if(abs(yDipoleExit)>yDipoleExitAbsMax) return false; // the track should point to |y|<yDipoleExitAbsMax at the dipole exit
	if(abs(xDipoleExit)<xDipoleExitAbsMin) return false; // the track should point to |x|>xDipoleExitAbsMin at the dipole exit
	if(abs(xDipoleExit)>xDipoleExitAbsMax) return false; // the track should point to |x|<xDipoleExitAbsMax at the dipole exit
	// if(abs(r4[0]-r1[0])>(fDxvsX->Eval(r4[0])*(1+0.1))) return false; // new cut!!
	// if(abs(r4[0]-r1[0])<(fDxvsX->Eval(r4[0])*(1-0.1))) return false; // new cut!!
	if(abs(r4[0]-r1[0])>5.5 || abs(r4[0]-r1[0])<0.5) return false; // new cut!!
	if(abs(r4[0]-r1[0])>(fDxvsX->Eval(r4[0])+0.1))   return false; // new cut!!
	if(abs(r4[0]-r1[0])<(fDxvsX->Eval(r4[0])-0.1))   return false; // new cut!!
	// if(!adaptive_dx14_vs_x4(r4[0],r1[0],fDxvsX,0.2,4)) return false;
	// if(!check_clusters(i1,i4,side))        return false; // minimum one cluster at layer 2 and one at layer 3 ///TODO!!!

	TRandom rnd;
	rnd.SetSeed();
	double posneg = rnd.Uniform(-1,+1);
	double pxgaus = rnd.Gaus(7.2e-4,5.0e-4);

	double xExit = abs(xofz(r1,r4,zDipoleExit)); // in cm!
	double P = fEvsX->Eval(r4[0]); // in GeV

	TVector2 v1(r1[2],r1[1]);
	TVector2 v4(r4[2],r4[1]);
	TVector2 u = rUnit2(v1,v4);
	double uz = u.X();
	double uy = u.Y();
	double px = (posneg>=0) ? pxgaus : -pxgaus;
	double py = P*uy;
	double pz = P*uz;
	p.SetPxPyPzE(px,py,pz,TMath::Sqrt(px*px + py*py + pz*pz + meGeV2));
	// if(i4==0 and side=="Eside") cout << "px=" << px << ", py=" << py << ", pz=" << pz << endl;
	// cout << "side=" << side << ", px=" << px << ", py=" << py << ", pz=" << pz << endl;
	float EseedMax = (process=="glaser") ? EseedMaxBPPP : EseedMaxTRIDENT;	// GeV
	if(p.E()<EseedMin or p.E()>EseedMax) return false;

	return true;
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



int main(int argc, char *argv[])
{	
	int argcounter; 
	printf("Program Name Is: %s",argv[0]);
	if(argc>=2) 
	{ 
		//gInterpreter->GenerateDictionary("vector<vector<TVector3> >",    "vector");
		printf("\nNumber Of Arguments Passed: %d",argc); 
		printf("\n----Following Are The Command Line Arguments Passed----"); 
		for(argcounter=0;argcounter<argc;argcounter++) printf("\nargv[%d]: %s",argcounter,argv[argcounter]);
		printf("\n");
	}
	//// minimum requirements
	if(argc<2) { printf("argc<2, exitting now\n"); exit(-1); }
	//// validate inputs
	if(argc==2 and !((TString)argv[1]).Contains("-proc=")) { printf("argc=2 but cannot parse %s\n",argv[1]); exit(-1); }
	if(argc==3 and !((TString)argv[2]).Contains("-dobg=")) { printf("argc=3 but cannot parse %s\n",argv[2]); exit(-1); }
	if(argc==4 and !((TString)argv[3]).Contains("-evnt=")) { printf("argc=4 but cannot parse %s\n",argv[3]); exit(-1); }
	if(argc==5 and !((TString)argv[4]).Contains("-seed=")) { printf("argc=5 but cannot parse %s\n",argv[4]); exit(-1); }
	if(argc==6 and !((TString)argv[5]).Contains("-ntrk=")) { printf("argc=6 but cannot parse %s\n",argv[5]); exit(-1); }
	//// assign inputs
	TString process = ((TString)argv[1]).ReplaceAll("-proc=",""); // mandatory
	int     dobg     = (argc>2) ? toint(((TString)argv[2]).ReplaceAll("-dobg=","")) : 0; // job id [optional]
	int     evnt     = (argc>2) ? toint(((TString)argv[3]).ReplaceAll("-evnt=","")) : -1; // job id [optional]
	int     Seed     = (argc>3) ? toint(((TString)argv[4]).ReplaceAll("-seed=","")) : 12345; // seed [optional]
	int     nsigtrks = (argc>4) ? toint(((TString)argv[5]).ReplaceAll("-ntrk=","")) : -1; // job id [optional]
	//// print assigned inputs
	cout << "process=" << process << endl;
	cout << "dobkg?=" << dobg << endl;
	cout << "evnt=" << evnt << endl;
	cout << "nsigtrks=" << nsigtrks << endl;
	cout << "Seed=" << Seed << endl;
	
	// TString proc = process;
	TString eventid = (evnt<0) ? "" : FormatEventID(evnt);
	TStopwatch stopwatch;
	
	
	/// get the B-field vs xExit functions to read off the 
	// TFile* fEvsx = new TFile(storage+"/output/root/test_bfield_fit.root","READ");
	// TF1* fEvsX_pos = (TF1*)fEvsx->Get("fEvsX_pos");
	// TF1* fEvsX_ele = (TF1*)fEvsx->Get("fEvsX_ele");
	
	TFile* fFits = new TFile(storage+"/output/root/test_bfield_fit2.root","READ");
	

	TF1* fEvsX_L1I_Eside = (TF1*)fFits->Get("h2_E_vs_x_L1I_Eside");
	TF1* fEvsX_L1I_Pside = (TF1*)fFits->Get("h2_E_vs_x_L1I_Pside");
	TF1* fEvsX_L4I_Eside = (TF1*)fFits->Get("h2_E_vs_x_L4I_Eside");
	TF1* fEvsX_L4I_Pside = (TF1*)fFits->Get("h2_E_vs_x_L4I_Pside");
	TF1* fDx14vsX_L4I_Eside = (TF1*)fFits->Get("h2_dx14_vs_x_L4I_Eside");
	TF1* fDx14vsX_L4I_Pside = (TF1*)fFits->Get("h2_dx14_vs_x_L4I_Pside");
	cout << "setup fits from files" << endl;

	int outN = (process=="elaser") ? 10 : 10;

	if(process=="elaser")
	{
		// resetToTridentGeometry();
		cout << "Doing only Pside!" << endl;
		sides.clear();
		sides.push_back("Pside"); /// do not reconstruct the Eside
	}
	
	// output tree
	cout << "Setting the output tree" << endl;
	gInterpreter->GenerateDictionary("vector<TLorentzVector>",    "vector");
	gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>",    "vector");
	gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",      "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",      "vector");
	gInterpreter->GenerateDictionary("vector<vector<TVector3> >", "vector");
	gSystem->Exec("mkdir -p "+storage+"/data/root/rec");
	
	
	/// monitoring histograms
	Int_t nlogebins = 30;
	Double_t logemin = 1.;
	Double_t logemax = 17.5;
	Double_t logebins[nlogebins+1];
	SetLogBins(nlogebins,logemin,logemax,logebins);
	TMapTSTH1D histos;
	TString hname = "";
	hname = "h_dErel_sed_gen_Eside"  ; histos.insert( make_pair(hname, new TH1D(hname,"Seed vs Gen;(E_{seed}-E_{gen})/E_{gen};Tracks",150,-0.03,+0.03)) );
	hname = "h_dErel_sed_gen_Pside"  ; histos.insert( make_pair(hname, new TH1D(hname,"Seed vs Gen;(E_{seed}-E_{gen})/E_{gen};Tracks",150,-0.03,+0.03)) );
	hname = "h_dErel_rec_gen_Eside"  ; histos.insert( make_pair(hname, new TH1D(hname,"Rec vs Gen;(E_{rec}-E_{gen})/E_{gen};Tracks",150,-0.03,+0.03)) );
	hname = "h_dErel_rec_gen_Pside"  ; histos.insert( make_pair(hname, new TH1D(hname,"Rec vs Gen;(E_{rec}-E_{gen})/E_{gen};Tracks",150,-0.03,+0.03)) );
	hname = "h_chi2_Eside"           ; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_chi2_Pside"           ; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_chi2_matched_Eside"   ; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_chi2_matched_Pside"   ; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_chi2_nonmatched_Eside"; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_chi2_nonmatched_Pside"; histos.insert( make_pair(hname, new TH1D(hname,";#chi^2;Tracks",150,0,15)) );
	hname = "h_E_tru_all_Eside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)) );
	hname = "h_E_tru_all_Pside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{all} [GeV];Tracks",68,0,17)) );
	hname = "h_E_tru_sed_mat_Eside"  ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{mat} [GeV];Tracks",68,0,17)) );
	hname = "h_E_tru_sed_mat_Pside"  ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{mat} [GeV];Tracks",68,0,17)) );
	hname = "h_E_tru_rec_mat_Eside"  ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{mat} [GeV];Tracks",68,0,17)) );
	hname = "h_E_tru_rec_mat_Pside"  ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru}^{mat} [GeV];Tracks",68,0,17)) );
	hname = "h_E_eff_sed_Eside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru} [GeV];Tracks",68,0,17)) );
	hname = "h_E_eff_sed_Pside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru} [GeV];Tracks",68,0,17)) );
	hname = "h_E_eff_rec_Eside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru} [GeV];Tracks",68,0,17)) );
	hname = "h_E_eff_rec_Pside"      ; histos.insert( make_pair(hname, new TH1D(hname,";#it{E}_{tru} [GeV];Tracks",68,0,17)) );
	
	
	
	
	/////////////////////
	/// loop on the sides
	for(unsigned int s=0 ; s<sides.size() ; ++s)
	{
		TString side = sides[s];
		
		TFile* fOut = new TFile(storage+"/data/root/rec/rec_"+process+"_"+eventid+"_"+side+".root","RECREATE");
		TTree* tOut = new TTree("reco","reco");
		/// all clusters output branches
		vector<TPolyMarker3D*>     all_clusters_xyz;
		vector<vector<TVector3> >  all_clusters_r;
		vector<int>                all_clusters_type;
		vector<int>                all_clusters_id;
		tOut->Branch("all_clusters_xyz",  &all_clusters_xyz);
		tOut->Branch("all_clusters_r",    &all_clusters_r);
		tOut->Branch("all_clusters_type", &all_clusters_type);
		tOut->Branch("all_clusters_id",   &all_clusters_id);
		/// truth output branches
		vector<int>            true_acc;
		vector<float>          true_wgt;
		vector<float>          true_x;
		vector<float>          true_y;
		vector<float>          true_z;
		vector<float>          true_q;
		vector<TLorentzVector> true_p;
		vector<TPolyMarker3D*> true_trckmar;
		vector<TPolyLine3D*>   true_trcklin;
		vector<vector<int> >   true_rec_imatch;
		vector<vector<int> >   true_clusters_id;
		tOut->Branch("true_acc",         &true_acc);
		tOut->Branch("true_wgt",         &true_wgt);
		tOut->Branch("true_x",           &true_x);
		tOut->Branch("true_y",           &true_y);
		tOut->Branch("true_z",           &true_z);
		tOut->Branch("true_q",           &true_q);
		tOut->Branch("true_p",           &true_p);
		tOut->Branch("true_trckmar",     &true_trckmar);
		tOut->Branch("true_trcklin",     &true_trcklin);
		tOut->Branch("true_rec_imatch",  &true_rec_imatch);
		tOut->Branch("true_clusters_id", &true_clusters_id);
		/// background tracks output branches
		vector<int>              bkgr_acc;
		vector<float>            bkgr_wgt;
		vector<float>            bkgr_x;
		vector<float>            bkgr_y;
		vector<float>            bkgr_z;
		vector<float>            bkgr_q;
		vector<TLorentzVector>   bkgr_p;
		vector<TPolyMarker3D*>   bkgr_trckmar;
		vector<TPolyLine3D*>     bkgr_trcklin;
		vector<vector<int> >     bkgr_clusters_id;
		tOut->Branch("bkgr_acc",         &bkgr_acc);
		tOut->Branch("bkgr_wgt",         &bkgr_wgt);
		tOut->Branch("bkgr_x",           &bkgr_x);
		tOut->Branch("bkgr_y",           &bkgr_y);
		tOut->Branch("bkgr_z",           &bkgr_z);
		tOut->Branch("bkgr_q",           &bkgr_q);
		tOut->Branch("bkgr_p",           &bkgr_p);
		tOut->Branch("bkgr_trckmar",     &bkgr_trckmar);
		tOut->Branch("bkgr_trcklin",     &bkgr_trcklin);
		tOut->Branch("bkgr_clusters_id", &bkgr_clusters_id);
		/// seeds output branches
		vector<int>            seed_type;
		vector<vector<int> >   seed_clusters_id;
		vector<float>          seed_q;
		vector<TLorentzVector> seed_p;
		tOut->Branch("seed_type",        &seed_type);
		tOut->Branch("seed_clusters_id", &seed_clusters_id);
		tOut->Branch("seed_q",           &seed_q);
		tOut->Branch("seed_p",           &seed_p);
		/// reconstructed clusters output branches
		vector<float>             reco_q;
		vector<TLorentzVector>    reco_p;
		vector<float>             reco_x;
		vector<float>             reco_y;
		vector<float>             reco_z;
		vector<vector<TVector3> > reco_trck_cls_r;
		vector<TPolyMarker3D*>   reco_trckmar;
		vector<TPolyLine3D*>     reco_trcklin;
		vector<float>            reco_chi2dof;
		vector<int>              reco_ismtchd;
		vector<int>              reco_ixmtchd;
		vector<int>              reco_idmtchd;
		vector<vector<int> >     reco_clusters_id;
		vector<double>           reco_Tgl;
		vector<double>           reco_Snp; // the slope in X direction: probe->GetTrack()->GetSnp()
		vector<double>           reco_alpha;
		vector<double>           reco_signedinvpT; // new: the curvature (q/Pyz): probe->GetTrack()->GetSigned1Pt()
		vector<double>           reco_sigmaY2;
		vector<double>           reco_sigmaZY;
		vector<double>           reco_sigmaZ2;
		vector<double>           reco_sigmaSnpY;
		vector<double>           reco_sigmaSnpZ;
		vector<double>           reco_sigmaSnp2; // probe->GetTrack()->GetSigmaSnp2()
		vector<double>           reco_sigmaTglY;
		vector<double>           reco_sigmaTglZ;
		vector<double>           reco_sigmaTglSnp;
		vector<double>           reco_sigmaTgl2;
		vector<double>           reco_sigma1PtY;
		vector<double>           reco_sigma1PtZ;
		vector<double>           reco_sigma1PtSnp;
		vector<double>           reco_sigma1PtTgl;
		vector<double>           reco_sigma1Pt2;
		vector<double>           reco_invpT;
		vector<double>           reco_signedpT;
		tOut->Branch("reco_q",           &reco_q);
		tOut->Branch("reco_p",           &reco_p);
		tOut->Branch("reco_x",           &reco_x);
		tOut->Branch("reco_y",           &reco_y);
		tOut->Branch("reco_z",           &reco_z);
		tOut->Branch("reco_trck_cls_r",  &reco_trck_cls_r);
		tOut->Branch("reco_trckmar",     &reco_trckmar);
		tOut->Branch("reco_trcklin",     &reco_trcklin);
		tOut->Branch("reco_chi2dof",     &reco_chi2dof);
		tOut->Branch("reco_ismtchd",     &reco_ismtchd);
		tOut->Branch("reco_ixmtchd",     &reco_ixmtchd);
		tOut->Branch("reco_idmtchd",     &reco_idmtchd);
		tOut->Branch("reco_clusters_id", &reco_clusters_id);
		tOut->Branch("reco_Tgl",         &reco_Tgl        );
		tOut->Branch("reco_Snp",         &reco_Snp        );
		tOut->Branch("reco_alpha",       &reco_alpha      );
		tOut->Branch("reco_signedinvpT", &reco_signedinvpT);
		tOut->Branch("reco_sigmaY2",     &reco_sigmaY2    );
		tOut->Branch("reco_sigmaZY",     &reco_sigmaZY    );
		tOut->Branch("reco_sigmaZ2",     &reco_sigmaZ2    );
		tOut->Branch("reco_sigmaSnpY",   &reco_sigmaSnpY  );
		tOut->Branch("reco_sigmaSnpZ",   &reco_sigmaSnpZ  );
		tOut->Branch("reco_sigmaSnp2",   &reco_sigmaSnp2  );
		tOut->Branch("reco_sigmaTglY",   &reco_sigmaTglY  );
		tOut->Branch("reco_sigmaTglZ",   &reco_sigmaTglZ  );
		tOut->Branch("reco_sigmaTglSnp", &reco_sigmaTglSnp);
		tOut->Branch("reco_sigmaTgl2",   &reco_sigmaTgl2  );
		tOut->Branch("reco_sigma1PtY",   &reco_sigma1PtY  );
		tOut->Branch("reco_sigma1PtZ",   &reco_sigma1PtZ  );
		tOut->Branch("reco_sigma1PtSnp", &reco_sigma1PtSnp);
		tOut->Branch("reco_sigma1PtTgl", &reco_sigma1PtTgl);
		tOut->Branch("reco_sigma1Pt2",   &reco_sigma1Pt2  );
		tOut->Branch("reco_invpT",       &reco_invpT      );
		tOut->Branch("reco_signedpT",    &reco_signedpT   );
		
	
		TString setup = "../setup/setupLUXE_"+process+"_"+side+".txt";
		det = new KMCDetectorFwd();
		det->ReadSetup(setup,setup);
		det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VertexTelescope
		// det->SetMinITSHits( det->GetNumberOfActiveLayersITS() ); // require hit in every layer
		det->SetMinITSHits( nMinHits ); // require hit in at least 4 layers //TODO!!!
		det->SetMinMSHits(0); // we don't have muon spectrometer
		det->SetMinTRHits(0); // we don't have muon trigger stations
		// max number of seeds on each layer to propagate (per muon track)
		det->SetMaxSeedToPropagate(3000); // relevant only if background is considered
		// set chi2 cuts
		// det->SetMaxChi2Cl(10.);  // max track to cluster chi2
		// det->SetMaxChi2Cl(10.);  // max track to cluster chi2
		det->SetMaxChi2Cl(10.);  // max track to cluster chi2
		// det->SetMaxChi2NDF(3.5); // max total chi2/ndf
		// det->SetMaxChi2NDF((process=="elaser")?15.:5.); // max total chi2/ndf
		// det->SetMaxChi2NDF((process=="elaser")?15.:5.); // max total chi2/ndf
		det->SetMaxChi2NDF((process=="elaser")?15.:5.); // max total chi2/ndf
		det->SetMaxChi2Vtx(20e9);  // fiducial cut on chi2 of convergence to vtx
		// det->SetMaxChi2Vtx(1e3);  // fiducial cut on chi2 of convergence to vtx
		// det->SetMaxChi2Vtx(500);  // fiducial cut on chi2 of convergence to vtx
		det->SetDefStepAir(1); // IMPORTANT FOR NON-UNIFORM FIELDS
		det->SetMinP2Propagate(0.3); //NA60+
		det->SetIncludeVertex(kTRUE); // count vertex as an extra measured point
		det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
		det->SetApplyBransonPCorrection(-1); // Branson correction, only relevant for setup with MS
		// for reconstruction:
		// det->SetErrorScale(500.);
		// det->SetErrorScale( (process=="elaser")?500.:200. );
		det->SetErrorScale( (process=="elaser")?500.:300. ); // was 400 earlier, can be also anywhere up to 1000
		det->Print();
		// det->BookControlHistos();
		
		///////////////////////////////
		setParametersFromDet(side); ///
		///////////////////////////////
		
		
		/// get the input signal clusters
		int                        sig_ngen          = 0;
		int                        sig_nslv          = 0;
		int                        sig_nacc          = 0;
		vector<double>*            sig_wgt           = 0;
		vector<int>*               sig_crg           = 0;
		vector<float>*             sig_xvtx          = 0;
		vector<float>*             sig_yvtx          = 0;
		vector<float>*             sig_zvtx          = 0;
		vector<TLorentzVector>*    sig_trkp4         = 0;
		vector<int>*               sig_acc           = 0;
		vector<vector<int> >*      sig_clusters_layerid = 0;
		vector<vector<int> >*      sig_clusters_id   = 0;
		vector<vector<int> >*      sig_clusters_type = 0;
		vector<TPolyMarker3D*>*    sig_clusters_xyz  = 0;
		vector<vector<TVector3> >* sig_clusters_r    = 0;
		vector<TPolyMarker3D*>*    sig_trkpts        = 0;
		vector<TPolyLine3D*>*      sig_trklin        = 0;
		cout << "Getting signal clusters from tree" << endl;
		TFile* fSig = new TFile(storage+"/data/root/dig/dig_"+process+"_"+eventid+".root","READ");
		TTree* tSig = (TTree*)fSig->Get("dig_"+side);
		tSig->SetBranchAddress("ngen",         &sig_ngen);
		tSig->SetBranchAddress("nslv",         &sig_nslv);
		tSig->SetBranchAddress("nacc",         &sig_nacc);
		tSig->SetBranchAddress("wgt",          &sig_wgt);
		tSig->SetBranchAddress("crg",          &sig_crg);
	   tSig->SetBranchAddress("xvtx",         &sig_xvtx);
	   tSig->SetBranchAddress("yvtx",         &sig_yvtx);
	   tSig->SetBranchAddress("zvtx",         &sig_zvtx);
		tSig->SetBranchAddress("trkp4",        &sig_trkp4);
		tSig->SetBranchAddress("acc",          &sig_acc);
		tSig->SetBranchAddress("clusters_layerid",  &sig_clusters_layerid);
		tSig->SetBranchAddress("clusters_id",  &sig_clusters_id);
		tSig->SetBranchAddress("clusters_type",&sig_clusters_type);
		tSig->SetBranchAddress("clusters_xyz", &sig_clusters_xyz);
		tSig->SetBranchAddress("clusters_r",   &sig_clusters_r);
		tSig->SetBranchAddress("trkpts",       &sig_trkpts);
		tSig->SetBranchAddress("trklin",       &sig_trklin);
		
		
		/// get the input background clusters
		int                        bkg_ngen          = 0;
		int                        bkg_nslv          = 0;
		int                        bkg_nacc          = 0;
		vector<double>*            bkg_wgt           = 0;
		vector<int>*               bkg_crg           = 0;
		vector<float>*             bkg_xvtx          = 0;
		vector<float>*             bkg_yvtx          = 0;
		vector<float>*             bkg_zvtx          = 0;
		vector<TLorentzVector>*    bkg_trkp4         = 0;
		vector<int>*               bkg_acc           = 0;
		vector<vector<int> >*      bkg_clusters_layerid = 0;
		vector<vector<int> >*      bkg_clusters_id   = 0;
		vector<vector<int> >*      bkg_clusters_type = 0;
		vector<TPolyMarker3D*>*    bkg_clusters_xyz  = 0;
		vector<vector<TVector3> >* bkg_clusters_r    = 0;
		vector<TPolyMarker3D*>*    bkg_trkpts        = 0;
		vector<TPolyLine3D*>*      bkg_trklin        = 0;
		TFile* fBkg = 0;
		TTree* tBkg = 0;
		if(dobg)
		{
			/// get the background clusters
			cout << "Getting background clusters from tree" << endl;
			// TChain* tBkg = new TChain("dig");
			// tBkg->Add(storage+"/data/root/dig_"+process+"_bkg_0*.root");
			// cout << "---- TChain content ----" << endl;
			// tBkg->ls();
			// cout << "------------------------" << endl;
			fBkg = new TFile(storage+"/data/root/dig/dig_"+process+"_bkg_"+eventid+".root","READ");
			tBkg = (TTree*)fBkg->Get("dig_"+side);
			tBkg->SetBranchAddress("ngen",         &bkg_ngen);
			tBkg->SetBranchAddress("nslv",         &bkg_nslv);
			tBkg->SetBranchAddress("nacc",         &bkg_nacc);
			tBkg->SetBranchAddress("wgt",          &bkg_wgt);
			tBkg->SetBranchAddress("crg",          &bkg_crg);
	   	tBkg->SetBranchAddress("xvtx",         &bkg_xvtx);
	   	tBkg->SetBranchAddress("yvtx",         &bkg_yvtx);
	   	tBkg->SetBranchAddress("zvtx",         &bkg_zvtx);
			tBkg->SetBranchAddress("trkp4",        &bkg_trkp4);
			tBkg->SetBranchAddress("acc",          &bkg_acc);
			tBkg->SetBranchAddress("clusters_layerid",  &bkg_clusters_layerid);
			tBkg->SetBranchAddress("clusters_id",  &bkg_clusters_id);
			tBkg->SetBranchAddress("clusters_type",&bkg_clusters_type);
			tBkg->SetBranchAddress("clusters_xyz", &bkg_clusters_xyz);
			tBkg->SetBranchAddress("clusters_r",   &bkg_clusters_r);
			tBkg->SetBranchAddress("trkpts",       &bkg_trkpts);
			tBkg->SetBranchAddress("trklin",       &bkg_trklin);
		}
	 
		/// prepare the dictionaries
		prepare_cahced_clusters();
	
		/// for timing
		Double_t av_cputime  = 0;
		Double_t av_realtime = 0;
	 
		/// loop on events
		Int_t nsigevents = tSig->GetEntries();
		Int_t nbkgevents = (dobg) ? tBkg->GetEntries() : -1;
		if(dobg && nbkgevents<nsigevents)
		{
			cout << "ERROR: nbkgevents<nsigevents" << endl;
			cout << "       nsigevents=" << nsigevents << endl;
			cout << "       nbkgevents=" << nbkgevents << endl;
			exit(-1);
		}
		cout << "Starting loop over signal events with nsigevents=" << nsigevents << endl;
		for(int iev=0 ; iev<tSig->GetEntries() ; iev++)
		// for(int iev=0 ; iev<tSig->GetEntries() ; iev++)
		{
			stopwatch.Start();
			
			//////////////////////////////
			//// get the next input entry
			tSig->GetEntry(iev); /// signal
			if(dobg) tBkg->GetEntry(iev); /// background
			
			////////////////////////////////////////////
			/// clear output vectors: digitized clusters
			for(unsigned int x=0 ; x<all_clusters_xyz.size() ; ++x)
			{
				delete all_clusters_xyz[x];
				all_clusters_r[x].clear();
			}
			all_clusters_xyz.clear();
			all_clusters_r.clear();
			all_clusters_type.clear();
			all_clusters_id.clear();
			/// clear output vectors: truth signal physics
			for(unsigned int x=0 ; x<true_rec_imatch.size() ; ++x) true_rec_imatch[x].clear();
			true_rec_imatch.clear();
			for(unsigned int x=0 ; x<true_clusters_id.size() ; ++x) true_clusters_id[x].clear();
			true_clusters_id.clear();
			true_acc.clear();
			true_wgt.clear();
			true_x.clear();
			true_y.clear();
			true_z.clear();
			true_q.clear();
			true_p.clear();
			true_trckmar.clear();
			true_trcklin.clear();
			/// clear output vectors: truth background physics
			bkgr_acc.clear();
			bkgr_wgt.clear();
			bkgr_x.clear();
			bkgr_y.clear();
			bkgr_z.clear();
			bkgr_q.clear();
			bkgr_p.clear();
			bkgr_trckmar.clear();
			bkgr_trcklin.clear();
			for(unsigned int x=0 ; x<bkgr_clusters_id.size() ; ++x) bkgr_clusters_id[x].clear();
			bkgr_clusters_id.clear();
			/// clear output vectors: seeds
			seed_type.clear();
			for(unsigned int x=0 ; x<seed_clusters_id.size() ; ++x) seed_clusters_id[x].clear();
			seed_clusters_id.clear();
			seed_q.clear();
			seed_p.clear();
			/// clear output vectors: reconstruction
			reco_q.clear();
			reco_p.clear();
			reco_x.clear();
			reco_y.clear();
			reco_z.clear();
			for(unsigned int x=0 ; x<reco_trck_cls_r.size() ; ++x) reco_trck_cls_r[x].clear();
			for(unsigned int x=0 ; x<reco_trckmar.size() ; ++x) delete reco_trckmar[x];
			for(unsigned int x=0 ; x<reco_trcklin.size() ; ++x) delete reco_trcklin[x];
			reco_trck_cls_r.clear();
			reco_trckmar.clear();
			reco_trcklin.clear();
			reco_chi2dof.clear();
			reco_ismtchd.clear();
			reco_ixmtchd.clear();
			reco_idmtchd.clear();
			for(unsigned int x=0 ; x<reco_clusters_id.size() ; ++x) reco_clusters_id[x].clear();
			reco_clusters_id.clear();
			reco_Tgl.clear();         
			reco_Snp.clear();
			reco_alpha.clear();
			reco_signedinvpT.clear(); 
			reco_sigmaY2.clear();
			reco_sigmaZY.clear();
			reco_sigmaZ2.clear();
			reco_sigmaSnpY.clear();
			reco_sigmaSnpZ.clear();
			reco_sigmaSnp2.clear();
			reco_sigmaTglY.clear();
			reco_sigmaTglZ.clear();
			reco_sigmaTglSnp.clear();
			reco_sigmaTgl2.clear();
			reco_sigma1PtY.clear();
			reco_sigma1PtZ.clear();
			reco_sigma1PtSnp.clear();
			reco_sigma1PtTgl.clear();
			reco_sigma1Pt2.clear();
			reco_invpT.clear();
			reco_signedpT.clear();
			
			//// clear cached clusters
			clear_cached_clusters(); /// clear for both sides
			
			/// rest all the layers of the detector (including inactive if any)
			reset_layers_all(); // reset both sides
			
			/// fill truth signal tracks:
			vector<int> vitmp;
			for(unsigned int t=0 ; t<sig_crg->size() ; ++t)
			{
				if(side=="Eside" && sig_crg->at(t)>0) continue;
				if(side=="Pside" && sig_crg->at(t)<0) continue;
				
				vector<int> vtruid;
				for(int k=0 ; k<sig_clusters_id->at(t).size() ; ++k) vtruid.push_back( sig_clusters_id->at(t)[k] );
				true_clusters_id.push_back( vtruid );
				true_acc.push_back( sig_acc->at(t) );
				true_wgt.push_back( sig_wgt->at(t) );
				true_x.push_back( sig_xvtx->at(t) );
				true_y.push_back( sig_yvtx->at(t) );
				true_z.push_back( sig_zvtx->at(t) );
				true_q.push_back( sig_crg->at(t) );
				true_p.push_back( sig_trkp4->at(t) );
				true_trckmar.push_back( sig_trkpts->at(t) );
				true_trcklin.push_back( sig_trklin->at(t) );
				true_rec_imatch.push_back( vitmp );
			}
			
			/// fill truth background tracks:
			int nbtrks = (dobg) ? (int)bkg_crg->size() : -1;
			int nbmax  = nbtrks; //(nMaxBkgTrks>0 && nMaxBkgTrks<nbtrks) ? nMaxBkgTrks : nbtrks;
			if(dobg)
			{
				for(int b=0 ; b<nbmax ; ++b)
				{
					if(!bkg_acc->at(b)) continue; // ignore tracks out of acceptance!
					
					vector<int> vbkgid;
					for(int k=0 ; k<bkg_clusters_id->at(b).size() ; ++k) vbkgid.push_back( bkg_clusters_id->at(b)[k] );
					bkgr_clusters_id.push_back( vbkgid );
					bkgr_acc.push_back( bkg_acc->at(b) );
					bkgr_wgt.push_back( bkg_wgt->at(b) );
					bkgr_x.push_back( bkg_xvtx->at(b) );
					bkgr_y.push_back( bkg_yvtx->at(b) );
					bkgr_z.push_back( bkg_zvtx->at(b) );
					bkgr_q.push_back( bkg_crg->at(b) );
					bkgr_p.push_back( bkg_trkp4->at(b) );
				   bkgr_trckmar.push_back( bkg_trkpts->at(b) );
				   bkgr_trcklin.push_back( bkg_trklin->at(b) );
				}
			}
			
			/// clear this side's indices
			cached_clusters_all_ids.clear();
			
			/// set the charge
			float crg = (side=="Eside") ? -1 : +1;
		   
			/// globals (per side)
			unsigned int n_truth = 0;
			unsigned int n_seeds = 0;
			unsigned int n_sedmt = 0;
			unsigned int n_solve = 0;
			unsigned int n_recos = 0;
			unsigned int n_match = 0;
			unsigned int n_trumt = 0;
			
			/// count truth per side
			for(unsigned int t=0 ; t<true_q.size() ; ++t)
			{
				if(side=="Eside" and true_q[t]>0) continue;
				if(side=="Pside" and true_q[t]<0) continue;
				n_truth++;
			}
		   
			/// make a pool of all signal clusters
			int ncached_signal_clusters = cache_clusters(sig_clusters_r,sig_clusters_type,sig_clusters_id,sig_clusters_layerid,side,nsigtrks);

			/// make a pool of all background and noise clusters
			int ncached_background_clusters = (dobg) ? cache_clusters(bkg_clusters_r,bkg_clusters_type,bkg_clusters_id,bkg_clusters_layerid,side): -1;
			
			/// rest all the layers of the detector (including inactive if any)
			reset_layers_all(); // reset both sides 
			
			/// add all clusters to the detector
			add_all_clusters(side);
			print_all_clusters(side,false);
			
			/// write out all clusters when these are sorted
			int all_clusters = fill_output_clusters(side,all_clusters_r,all_clusters_type,all_clusters_id);
			// cout << "ncached_signal_clusters=" << ncached_signal_clusters << ", ncached_background_clusters=" << ncached_background_clusters << ", all_clusters=" << all_clusters << endl;
			
			/// offset for signal id's !!!
			// int sigoffset = 100000; // should be multiplied by the layer number
			int index_offset = (process.Contains("bkg")) ? index_offset_bkg : index_offset_sig; // assuming no chance to have >index_offset tracks (and hence clusters) per tracker arm (and hence per layer)
			
			
			/// run over all clusters of layer 4 in the pool --> these are the seeds for the KalmanFilter fit
			TString slyr4 = (side=="Eside") ? "EL4I" : "PL4I";
			TString slyr1 = (side=="Eside") ? "EL1I" : "PL1I";
			int     ilyr4 = silayers[slyr4];
			int     ilyr1 = silayers[slyr1];
			// cout << "ilyr1=" << ilyr1 << ", ilyr4=" << ilyr4 << endl;
			
			/// loop on seeds
			for(unsigned int i4=0 ; i4<cached_clusters[slyr4].size() ; ++i4)
			// for(unsigned int i4=0 ; i4<1 ; ++i4)
			{
				// reset all tracks from all layers
				reset_layers_tracks();
				
				vector<TLorentzVector> pseeds;
				for(unsigned int i1=0 ; i1<cached_clusters[slyr1].size() ; ++i1)
				{	
					// reset all tracks from all layers but layer 0
					reset_layers_tracks(0);
					
					/// find the momentum of the seed
					TLorentzVector pseed;
					float r1[3]; r1[0]=cached_clusters[slyr1][i1].r.X(); r1[1]=cached_clusters[slyr1][i1].r.Y(); r1[2]=cached_clusters[slyr1][i1].r.Z();
					float r4[3]; r4[0]=cached_clusters[slyr4][i4].r.X(); r4[1]=cached_clusters[slyr4][i4].r.Y(); r4[2]=cached_clusters[slyr4][i4].r.Z();
					
					bool seed = makeseed_nonuniformB(process,r1,r4,i1,i4,side,pseed, (side=="Pside")?fEvsX_L4I_Pside:fEvsX_L4I_Eside, (side=="Pside")?fDx14vsX_L4I_Pside:fDx14vsX_L4I_Eside );
					if(!seed) continue; // cannot make a meaningful seed
					pseeds.push_back(pseed);
					bool issig  = ( cached_clusters[slyr1][i1].type==1 && cached_clusters[slyr4][i4].type==1 );
					bool sameid = ( (cached_clusters[slyr1][i1].clsid-ilyr1*index_offset_sig)==(cached_clusters[slyr4][i4].clsid-ilyr4*index_offset_sig) );
					seed_type.push_back( issig and sameid );
					vector<int> vidseed{cached_clusters[slyr1][i1].clsid,-1,-1,cached_clusters[slyr4][i4].clsid};
					seed_clusters_id.push_back(vidseed);
					seed_q.push_back(crg);
					seed_p.push_back(pseed);
					n_seeds++;
				} // end of loop on clusters in layer 1
				if(n_seeds<1) continue;
				cout << "nseeds=" << n_seeds << " for i4=" << i4 << " out of " << cached_clusters[slyr4].size() << " clusters in layer4 (with " << n_recos << " recos)" << endl;
				for(int d=0 ; d<seed_type.size() ; ++d) {if(seed_type[d]) cout << "at least one true seed" << endl; break;}
				
				bool doPrint = false;
				if(doPrint) cout << "\n\n\n########################################## calling SolveSingleTrackViaKalmanMC_Noam_multiseed for i4=" << i4 << " ######################################" << endl;
				// prepare the probe from the seed and do the KF fit
				
				bool solved = det->SolveSingleTrackViaKalmanMC_Noam_multiseed(pseeds,meGeV,crg,99,doPrint);
				if(!solved) continue; // reconstruction failed
				n_solve++;
				
				// get the reconstructed propagated to the vertex 
				KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack();
				if(!trw)            continue; // track was not reconstructed
				if(trw->IsKilled()) continue; // track was killed
				n_recos++;
				
				if(doPrint) {cout << "Track fit succeeded, associated clusters are:" << endl; trw->Print("clid");}
				
				/// save clusters of the track
				vector<TVector3> v3tmp;
				reco_trck_cls_r.push_back( v3tmp );
				
				/// the rec index - important to set it here!
				unsigned int irec = reco_trck_cls_r.size()-1;
				
				/// get the clusters of the winner tracK
				vector<int> win_cls_id;
				vector<int> win_cls_inx;
				TMapii win_cls_id2lr;
				const int* probeclusters = trw->GetClID();
				int nprobeclusters = sizeof(probeclusters);
				// cout << "nprobeclusters=" << nprobeclusters << endl;
				// trw->Print("clid etp");
				trw->Print("clid");
				for(int l=0 ; l<=nprobeclusters ; ++l)
				{
					int cid = probeclusters[l];
					if(cid<0) continue;
					int cix = cached_clusters_all_ids[cid];
					int ilr = cached_clusters_id2lyr[cid]; /// essential to find the correct layer!
					TString lname = layers[ilr];
					win_cls_id.push_back( cid ); // provide active layer ID, not the physical ones (most are passive)
					win_cls_inx.push_back( cix );
					win_cls_id2lr.insert(make_pair(cid,ilr));
					if(doPrint) cout << "going to kill cluster id=" << cid << " in ix=" << cix << ", from " << det->GetLayer(ilr)->GetNBgClusters() << " on layer " << ilr << endl;
					// cout << "cid=" << cid << " in layer=" << ilr << ": " << lname << endl;
					double xwin = det->GetLayer(ilr)->GetBgCluster(cix)->GetXLab();
					double ywin = det->GetLayer(ilr)->GetBgCluster(cix)->GetYLab();
					double zwin = det->GetLayer(ilr)->GetBgCluster(cix)->GetZLab();
					reco_trck_cls_r[irec].push_back(TVector3(xwin,ywin,zwin)); // fill before killing!
					det->GetLayer(ilr)->GetBgCluster(cix)->Kill();
				}
				/// save the clusters' id of the winner track 
				reco_clusters_id.push_back( win_cls_id );
				
				/// reco kinematics etc
				TLorentzVector prec;
				double pxyz[3];
				double xyz[3];
				trw->GetPXYZ(pxyz);
				trw->GetXYZ(xyz);
				prec.SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
				float chi2dof = trw->GetNormChi2();
				reco_chi2dof.push_back( chi2dof );
				reco_q.push_back( crg );
				reco_p.push_back( prec );
				reco_x.push_back( xyz[0] );
				reco_y.push_back( xyz[1] );
				reco_z.push_back( xyz[2] );
				reco_trckmar.push_back( TrackMarker3d(trw,0,zLastLayer+1,0.1,trkcol(prec.E())) );
				reco_trcklin.push_back( TrackLine3d(trw,zLastLayer+1,1,trkcol(prec.E())) );


				/// rec-tru matching
				int ismatched =  0;
				int ixmatched = -1;
				int idmatched = -1;
				int imatch = -1;
				for(TMapii::iterator it=win_cls_id2lr.begin() ; it!=win_cls_id2lr.end() ; ++it)
				{
					int cid = it->first;
					int lid = it->second;
					int ix = cid-lid*index_offset_sig;
					if(it==win_cls_id2lr.begin()) imatch = ix;
					if(ix!=imatch)
					{
						imatch = -1;
						break;
					}
					// cout << "lid=" << lid << " cid=" << cid << endl;
				}
				if(imatch>=0)
				{
					cout << "found match index: " << imatch << ", Etru=" << true_p[imatch].E() << ", Erec=" << prec.E() << endl;
					ismatched = 1;
					ixmatched = imatch;
					idmatched = true_clusters_id[imatch][0];
					true_rec_imatch[imatch].push_back( irec );
					n_match++;
				}
				else { ismatched =  0; ixmatched = -1; idmatched = -1; }
				reco_ismtchd.push_back( ismatched );
				reco_ixmtchd.push_back( ixmatched );
				reco_idmtchd.push_back( idmatched );
				cout << "n_recos=" << n_recos << ", n_match=" << n_match << endl;
				
				
				/// more kinematics
				TrackPar* trk = trw->GetTrack();
				reco_Tgl.push_back( trk->GetTgl() );
				reco_Snp.push_back( trk->GetSnp() );
				reco_alpha.push_back( trk->GetAlpha() );
				reco_signedinvpT.push_back( trk->GetSigned1Pt() );
				reco_sigmaY2.push_back( trk->GetSigmaY2() );
				reco_sigmaZY.push_back( trk->GetSigmaZY() );
				reco_sigmaZ2.push_back( trk->GetSigmaZ2() );
				reco_sigmaSnpY.push_back( trk->GetSigmaSnpY() );
				reco_sigmaSnpZ.push_back( trk->GetSigmaSnpZ() );
				reco_sigmaSnp2.push_back( trk->GetSigmaSnp2() );
				reco_sigmaTglY.push_back( trk->GetSigmaTglY() );
				reco_sigmaTglZ.push_back( trk->GetSigmaTglZ() );
				reco_sigmaTglSnp.push_back( trk->GetSigmaTglSnp() );
				reco_sigmaTgl2.push_back( trk->GetSigmaTgl2() );
				reco_sigma1PtY.push_back( trk->GetSigma1PtY() );
				reco_sigma1PtZ.push_back( trk->GetSigma1PtZ() );
				reco_sigma1PtSnp.push_back( trk->GetSigma1PtSnp() );
				reco_sigma1PtTgl.push_back( trk->GetSigma1PtTgl() );
				reco_sigma1Pt2.push_back( trk->GetSigma1Pt2() );
				reco_invpT.push_back( trk->OneOverPt() );
				reco_signedpT.push_back( trk->GetSignedPt() );
				
				pseeds.clear(); /// this is maybe redundant
			} // end of loop on clusters in layer 4
			
			///////////////////////////////////////////////
			/// post-processing per side histos to fill ///
			///////////////////////////////////////////////
			
			/// seed matching
			for(unsigned int t=0 ; t<true_q.size() ; ++t)
			{
				if(side=="Eside" and true_q[t]>0) continue;
				if(side=="Pside" and true_q[t]<0) continue;
				if(true_rec_imatch[t].size()>0) n_trumt++;
				// else
				// {
				// 	cout << "This truth track is not matched: Etru[" << t << "]=" << true_p[t].E() << " GeV with bestdistance=" << bestdistance << ", bestmatchtrui=" << bestmatchtrui << ", bestmatchreci=" << bestmatchreci << endl;
				// }
				histos["h_E_tru_all_"+side]->Fill( true_p[t].E() );
				int truid1 = true_clusters_id[t][0];
				int truid4 = true_clusters_id[t][3];
				// cout << "truid1=" << truid1 << ", truid4=" << truid4 << endl;
				for(unsigned int s=0 ; s<seed_p.size() ; ++s)
				{
					if(seed_type[s]!=1)                continue; // has to be signal track
					// cout << "seed_clusters_id["<<s<<"][0]=" << seed_clusters_id[s][0] << ", seed_clusters_id["<<s<<"][3]=" << seed_clusters_id[s][3] << endl;
					if(seed_clusters_id[s][0]!=truid1) continue; // match cluster id of layer 1
					if(seed_clusters_id[s][3]!=truid4) continue; // match cluster id of layer 4 
					histos["h_dErel_sed_gen_"+side]->Fill((seed_p[s].E()-true_p[t].E())/true_p[t].E());
					histos["h_E_tru_sed_mat_"+side]->Fill(seed_p[s].E());
					n_sedmt++;
					break;
				}
			}
			
			/// TODO: add a vector for all truth tracks, to have an inner vector of all matched reco tracks.
			/// TODO: then need to check if the truth track has more than 1 reco track and take the better one when filling.
			vector<int> ixtrumatched;
			for(unsigned int k=0 ; k<reco_ismtchd.size() ; ++k)
			{
				if(side=="Eside" and reco_q[k]>0) continue;
				if(side=="Pside" and reco_q[k]<0) continue;
				
				histos["h_chi2_"+side]->Fill( reco_chi2dof[k] ); // fill regardless of matching
				
				// /// TODO: now I skip if more than one tru track matched (later implement something to take the best one)
				// if(true_rec_imatch[reco_ixmtchd[k]].size()>1) continue;
				
				if(reco_ismtchd[k]==1 and reco_ixmtchd[k]>=0 and !foundinvec(reco_ixmtchd[k],ixtrumatched))
				{	
					ixtrumatched.push_back( reco_ixmtchd[k] ); /// fill and check in next iterations to avoid repetition

					histos["h_chi2_matched_"+side]->Fill( reco_chi2dof[k] );
					
					// bool accept = (reco_p[k].E()>1. and reco_p[k].E()<17.5);
					// if(!accept) continue;
					histos["h_E_tru_rec_mat_"+side]->Fill( sig_trkp4->at(reco_ixmtchd[k]).E() );
					histos["h_dErel_rec_gen_"+side]->Fill( (reco_p[k].E()-sig_trkp4->at(reco_ixmtchd[k]).E())/sig_trkp4->at(reco_ixmtchd[k]).E() );
				}
				else
				{
					histos["h_chi2_nonmatched_"+side]->Fill( reco_chi2dof[k] );
				}
			}
			
			
			/// summarize
			int mateff = (int)((float)n_trumt/(float)n_truth*100.);
			cout << "Event #" << iev << ", "<< side << ": n_truth=" << n_truth
				<< ", n_seeds=" << n_seeds
					<< ", n_sedmt=" << n_sedmt
						<< ", n_solve=" << n_solve
							<< ", n_recos=" << n_recos
								<< ", n_match=" << n_match
									<< ", n_trumt=" << n_trumt
										<< ", eff(rec,mat)=" << mateff << "%"<< endl;
		
			fOut->cd();
			tOut->Fill();
			
			stopwatch.Stop();
			Double_t cputime  = stopwatch.CpuTime();
			Double_t realtime = stopwatch.RealTime();
			av_cputime  += cputime;
			av_realtime += realtime;
			cout << "Event #" << iev << ": CPU time=" << cputime << ", Real time=" << realtime << endl;
			if((iev%outN)==0) printf("Done %d out of %d --> CPUav=%g, REAL=%g\n",iev,nsigevents,av_cputime/(iev+1),av_realtime/(iev+1));
		} // end of loop on events
	
		histos["h_E_eff_sed_Eside"]->Divide(histos["h_E_tru_sed_mat_Eside"],histos["h_E_tru_all_Eside"]);
		histos["h_E_eff_sed_Pside"]->Divide(histos["h_E_tru_sed_mat_Pside"],histos["h_E_tru_all_Pside"]);
		histos["h_E_eff_rec_Eside"]->Divide(histos["h_E_tru_rec_mat_Eside"],histos["h_E_tru_all_Eside"]);
		histos["h_E_eff_rec_Pside"]->Divide(histos["h_E_tru_rec_mat_Pside"],histos["h_E_tru_all_Pside"]);
		fOut->cd();
		tOut->Write();
		for(TMapTSTH1D::iterator it=histos.begin() ; it!=histos.end() ; ++it) it->second->Write();
		fOut->Write();
		fOut->Close();
	} // end of loop on sides
	
	return 0;
}
