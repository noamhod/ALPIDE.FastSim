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
#include <filesystem>
#include <map>
#include <sstream>
#include <limits> 
#endif


/*
1.   cluster id - DONE
2.   change the geometry (setup file, and code) - DONE
3.   deal with the gaps between the chips (for old digitization code) - NOT URGENT
4.   change the output of old digitization code same as clusters tree of Allpix - NOT URGENT
5-.  fDx14vsXMap -- change key names like the DyvsY map - DONE
6.   PropagatetoBxByBz in the track class very finely in the B-field in KF (out of the B-field this can very coarse) - NOT URGENT
7.   check energy loss in the material of the window/air - NOT URGENT
8+.  do the matching properly (from itru and all the other clusters in the winner track) (only for itru>=0) - DONE
8a+. calculate the matched counting efficiency (all reconstructed tracks) - DONE
8b+. calculate the selected matched counting efficiency (only selected tracks) - DONE
9-.  change all histograms to h_x to h_x_rec (only those that contain kinematic stuff) - DONE
10-. add more histograms for selected: h_sel_x (same kinematic block) - DONE
11-. plot signal only, background only and combined in different colours in all the variables (in 9 above) to identify potential cuts - DONE
12+. change the hardcoded path for the allpix digitized file - DONE
13+. formulate a set of cuts and fill h_sel_x - DONE
14. ask Sasha to generate flat signal
*/

using namespace std;

TString storage = gSystem->ExpandPathName("$STORAGEDIR");

bool debug = false;
bool debug2 = false;

struct Cluster
{
	int type;
	int lyridKF;
	TString lyrnmKF;
	long clsid;
	int npixl;
	int npixlx;
	int npixly;
	int shape;
	double charge;
	TVector3 r;
	int lyr;
	int chip;
	int cellx;
	int celly;
	vector<int> trksid;
	vector<int> trkstype;
	vector<TLorentzVector> trksp;
	int issig;
};



typedef std::numeric_limits< double>dbl;
typedef map<int, int> TMapii;
typedef map<TString, int> TMapTSi;
typedef map<TString, float> TMapTSf;
typedef map<TString, double> TMapTSd;
typedef map<TString, vector<int>> TMapTSvi;
typedef map<TString, vector<float>> TMapTSvf;
typedef map<int, TString> TMapiTS;
typedef map<TString, TH1D *> TMapTSTH1D;
typedef map<TString, TH2D *> TMapTSTH2D;
typedef map<TString, vector<Cluster>> TMapTSvCls;		  // formatted as
typedef map<TString, map<int, vector<int>>> TMapTSMapivi; // this is needed for the lookup table
typedef map<TString, TAxis *> TMapTSAxis;
typedef map<TString, map<int, int>> TMapTSMapii;
typedef map<TString, TF1 *> TMapTSTF1;
typedef map<TString, TString> TMapTSTS;
typedef map<int, TMapTSTS> TMapiTMapTSTS;




enum seedcuts{
	FAIL_SLOPE=-1,
	FAIL_SIDE=-2,
	FAIL_SAMEZ=-3,
	FAIL_DIPEXITYABS=-4,
	FAIL_DIPEXITXLOW=-5,
	FAIL_DIPEXITXHIGH=-6,
	FAIL_DXSIZE=-7,
	FAIL_RWXHIGH=-8,
	FAIL_RWXLOW=-9,
	FAIL_RWYHIGH=-10,
	FAIL_RWYLOW=-11,
	FAIL_NHITS3=-12,
	FAIL_NHITS2=-13,
	FAIL_NHITS1=-14,
	FAIL_E=-15,
	PASS=1
};

TMapiTS seedcutnames = {{FAIL_SLOPE,"slope"}, {FAIL_SIDE,"side"}, {FAIL_SAMEZ,"same z"}, {FAIL_DIPEXITYABS,"dipole exit |y|"},
								{FAIL_DIPEXITXLOW,"dipole exit |xmin|"}, {FAIL_DIPEXITXHIGH,"dipole exit |xmax|"},
								{FAIL_DXSIZE,"dx14 size"}, {FAIL_RWXHIGH,"road width x+dx"}, {FAIL_RWXLOW,"road width x-dx"},
								{FAIL_RWYHIGH,"road width y+dy"}, {FAIL_RWYLOW,"road width y-dy"},
								{FAIL_NHITS3,"insufficient hits in L3"}, {FAIL_NHITS2,"insufficient hits in L2"}, {FAIL_NHITS1,"insufficient hits in L1"},
								{FAIL_E,"energy"},{PASS,"pass"}};


int nMinHits = 4;

double vX = 0, vY = 0, vZ = 0; // event vertex
KMCDetectorFwd *det = 0;

double meMeV  = 0.5109989461; // MeV
double meGeV  = meMeV / 1000.;
double meGeV2 = meGeV * meGeV;
double cm2m   = 0.01;
double um2cm  = 0.0001;

// for matching these must be the same as in digitisation
int index_offset_bkg = 100000;
int index_offset_sig = 10000000;

//// temp stuff:
double bestdistance = -1;
int bestmatchtrui = -1;
int bestmatchreci = -1;

//// uncertainties
float dxAlignmentXFEL = 0.005; // 0.005; // cm
float dyAlignmentXFEL = 0.005; // 0.005; // cm
float XvariationSign = +1.;
float YvariationSign = +1.;
float dxAlignmentInTray = 0.000; // cm
float dyAlignmentInTray = 0.000; // cm
bool doMisalignmentX = false;
bool doMisalignmentY = false;

//// seed energies
double EseedMinGLaser = 0.5;  // GeV
double EseedMinELaser = 0.5;  // GeV
double EseedMaxGLaser = 14.0; // GeV
double EseedMaxELaser = 12.0; // GeV

vector<TString> sides{"Eside", "Pside"};
// vector<TString> sides{"Pside", "Eside"};
vector<TString> layersnames;
vector<double> layersz;
vector<double> zlayer;
TMapiTS layers;
TMapTSi szlayers;
TMapTSi silayers;
TMapiTS islayers;

TMapTSvCls cached_clusters;		/// formatted per side per layer e.g. as: cached_clusters[side+"_"+layerid][i].r.X() or cached_clusters[side+"_"+layerid][i].clsid
TMapii cached_clusters_all_ids; ///
TMapii cached_clusters_id2lyr;	///
TMapTSMapii cached_clusters_id2ix;

/// lookup table for cluster id's bin
TMapTSMapivi lookupTable;
TMapTSAxis axisMap;
int nAxisBins = 10000;

/// road/cone width
double rwxL1 = 0.030; // in cm, road width in x along the possible cluster where we embed the clusters
double rwxL2 = 0.025; // in cm, road width in x along the possible cluster where we embed the clusters
double rwxL3 = 0.020; // in cm, road width in x along the possible cluster where we embed the clus2s
double rwyL1 = 0.030; // in cm, road width in y along the possible cluster where we embed the clusters
double rwyL2 = 0.025; // in cm, road width in y along the possible cluster where we embed the clusters
double rwyL3 = 0.020; // in cm, road width in y along the possible cluster where we embed the clusters
TF2* fRWxPI = 0;
TF2* fRWxPO = 0;
TF2* fRWxEI = 0;
TF2* fRWxEO = 0;
TF2* fRWy   = 0;


//// staves geometry
double Hstave = -999;		
double Lstave = -999;			
double xPsideL = -999;
double xPsideR = -999;
double xEsideL = -999;
double xEsideR = -999;
double yUp = -999;
double yDn = -999;

//// dipole geometry
double xWdipole = -999;
double yHdipole = -999;
double z1dipole = -999;
double z2dipole = -999;
double zDipoleExit = -999;

// dipole field
double B = -999;
double LB = -999; // meters

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

TMapTSi staveNumbers;


double zLastLayer = -999;
double zFirstLayer = -999;
TString LastLayer = "";
TString FirstLayer = "";

int toint(TString str)
{
	stringstream strm;
	int x;
	strm << str;
	strm >> x;
	return x;
}
int tofloat(TString str)
{
	stringstream strm;
	float x;
	strm << str;
	strm >> x;
	return x;
}
string tostring(int n)
{
	stringstream strm;
	string str;
	strm << n;
	strm >> str;
	return str;
}

void setParametersFromDet(TString side, TString proc)
{
	cout << "====================================" << endl;
	cout << "============DEFINITIONS=============" << endl;
	cout << "====================================" << endl;
	KMCLayerFwd *layer_outer = (side=="Eside") ? det->GetLayer("EL1O") : det->GetLayer("PL1O");
	KMCLayerFwd *layer_inner = (side=="Eside") ? det->GetLayer("EL1I") : det->GetLayer("PL1I");

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

	yUp = layer_outer->GetYMax();
	yDn = layer_outer->GetYMin();
	cout << "yUp=" << yUp << ", yDn=" << yDn << endl;

	//// get the Bfield from the setup
	TVirtualMagField *fld = TGeoGlobalMagField::Instance()->GetField();
	MagField *fldm = (MagField *)fld;
	const double BfieldKG = fldm->GetBVals(0, 1);			   // region 0, the B field is only in the y direction (0,B,0), hence the index 1
	const double *BfieldXminObj = fldm->GetXMin();			   // region 0
	const double *BfieldXmaxObj = fldm->GetXMax();			   // region 0
	const double *BfieldYminObj = fldm->GetYMin();			   // region 0
	const double *BfieldYmaxObj = fldm->GetYMax();			   // region 0
	const double *BfieldZminObj = fldm->GetZMin();			   // region 0
	const double *BfieldZmaxObj = fldm->GetZMax();			   // region 0
	const std::string BFunction = fldm->GetFunctionForm(0, 1); /// 1 is for y component of the field
	double BfieldValTesla = BfieldKG / 10;					   /// the B field is only in the y direction (0,B,0), hence the index 1
	double BfieldXmin = BfieldXminObj[0];					   // region 0
	double BfieldXmax = BfieldXmaxObj[0];					   // region 0
	double BfieldYmin = BfieldYminObj[0];					   // region 0
	double BfieldYmax = BfieldYmaxObj[0];					   // region 0
	double BfieldZmin = BfieldZminObj[0];					   // region 0
	double BfieldZmax = BfieldZmaxObj[0];					   // region 0
	cout << "BfieldValTesla=" << BfieldValTesla << ", BfieldXmin=" << BfieldXmin << ", BfieldXmax=" << BfieldXmax << ", BfieldYmin=" << BfieldYmin << ", BfieldYmax=" << BfieldYmax << ", BfieldZmin=" << BfieldZmin << ", BfieldZmax=" << BfieldZmax << endl;
	cout << "Bfunction: " << BFunction << std::endl;
	xWdipole = BfieldXmax-BfieldXmin;
	yHdipole = BfieldYmax-BfieldYmin;
	z1dipole = BfieldZmin;
	z2dipole = BfieldZmax;
	zDipoleExit = z2dipole;
	B = BfieldValTesla;
	LB = z2dipole-z1dipole;
	cout << "xWdipole=" << xWdipole << ", yHdipole=" << yHdipole << ", z1dipole=" << z1dipole << ", z2dipole=" << z2dipole << endl;

	zEL1I = (side=="Eside") ? det->GetLayer("EL1I")->GetZ()  : -999;
	iEL1I = (side=="Eside") ? det->GetLayer("EL1I")->GetID() : -999;
	zEL1O = (side=="Eside") ? det->GetLayer("EL1O")->GetZ()  : -999;
	iEL1O = (side=="Eside") ? det->GetLayer("EL1O")->GetID() : -999;
	zPL1I = (side=="Pside") ? det->GetLayer("PL1I")->GetZ()  : -999;
	iPL1I = (side=="Pside") ? det->GetLayer("PL1I")->GetID() : -999;
	zPL1O = (side=="Pside") ? det->GetLayer("PL1O")->GetZ()  : -999;
	iPL1O = (side=="Pside") ? det->GetLayer("PL1O")->GetID() : -999;
	cout << "zEL1I=" << zEL1I << ", zEL1O=" << zEL1O << ", zPL1I=" << zPL1I << ", zPL1O=" << zPL1O << endl;

	zEL2I = (side=="Eside") ? det->GetLayer("EL2I")->GetZ()  : -999;
	iEL2I = (side=="Eside") ? det->GetLayer("EL2I")->GetID() : -999;
	zEL2O = (side=="Eside") ? det->GetLayer("EL2O")->GetZ()  : -999;
	iEL2O = (side=="Eside") ? det->GetLayer("EL2O")->GetID() : -999;
	zPL2I = (side=="Pside") ? det->GetLayer("PL2I")->GetZ()  : -999;
	iPL2I = (side=="Pside") ? det->GetLayer("PL2I")->GetID() : -999;
	zPL2O = (side=="Pside") ? det->GetLayer("PL2O")->GetZ()  : -999;
	iPL2O = (side=="Pside") ? det->GetLayer("PL2O")->GetID() : -999;
	cout << "zEL2I=" << zEL2I << ", zEL2O=" << zEL2O << ", zPL2I=" << zPL2I << ", zPL2O=" << zPL2O << endl;

	zEL3I = (side=="Eside") ? det->GetLayer("EL3I")->GetZ()  : -999;
	iEL3I = (side=="Eside") ? det->GetLayer("EL3I")->GetID() : -999;
	zEL3O = (side=="Eside") ? det->GetLayer("EL3O")->GetZ()  : -999;
	iEL3O = (side=="Eside") ? det->GetLayer("EL3O")->GetID() : -999;
	zPL3I = (side=="Pside") ? det->GetLayer("PL3I")->GetZ()  : -999;
	iPL3I = (side=="Pside") ? det->GetLayer("PL3I")->GetID() : -999;
	zPL3O = (side=="Pside") ? det->GetLayer("PL3O")->GetZ()  : -999;
	iPL3O = (side=="Pside") ? det->GetLayer("PL3O")->GetID() : -999;
	cout << "zEL3I=" << zEL3I << ", zEL3O=" << zEL3O << ", zPL3I=" << zPL3I << ", zPL3O=" << zPL3O << endl;

	zEL4I = (side=="Eside") ? det->GetLayer("EL4I")->GetZ()  : -999;
	iEL4I = (side=="Eside") ? det->GetLayer("EL4I")->GetID() : -999;
	zEL4O = (side=="Eside") ? det->GetLayer("EL4O")->GetZ()  : -999;
	iEL4O = (side=="Eside") ? det->GetLayer("EL4O")->GetID() : -999;
	zPL4I = (side=="Pside") ? det->GetLayer("PL4I")->GetZ()  : -999;
	iPL4I = (side=="Pside") ? det->GetLayer("PL4I")->GetID() : -999;
	zPL4O = (side=="Pside") ? det->GetLayer("PL4O")->GetZ()  : -999;
	iPL4O = (side=="Pside") ? det->GetLayer("PL4O")->GetID() : -999;
	cout << "zEL4I=" << zEL4I << ", zEL4O=" << zEL4O << ", zPL4I=" << zPL4I << ", zPL4O=" << zPL4O << endl;

	if(side=="Eside")
	{
		layersnames = {"EL1I", "EL1O", "EL2I", "EL2O", "EL3I", "EL3O", "EL4I", "EL4O"};
		zlayer = {0, z1dipole, z2dipole, zEL1I, zEL1O, zEL2I, zEL2O, zEL3I, zEL3O, zEL4I, zEL4O};
		layersz = {zEL1I, zEL1O, zEL2I, zEL2O, zEL3I, zEL3O, zEL4I, zEL4O};
		layers = {{iEL1I, "EL1I"}, {iEL1O, "EL1O"}, {iEL2I, "EL2I"}, {iEL2O, "EL2O"}, {iEL3I, "EL3I"}, {iEL3O, "EL3O"}, {iEL4I, "EL4I"}, {iEL4O, "EL4O"}};
		szlayers = {{"EL1I", zEL1I}, {"EL1O", zEL1O}, {"EL2I", zEL2I}, {"EL2O", zEL2O}, {"EL3I", zEL3I}, {"EL3O", zEL3O}, {"EL4I", zEL4I}, {"EL4O", zEL4O}};
		silayers = {{"EL1I", iEL1I}, {"EL1O", iEL1O}, {"EL2I", iEL2I}, {"EL2O", iEL2O}, {"EL3I", iEL3I}, {"EL3O", iEL3O}, {"EL4I", iEL4I}, {"EL4O", iEL4O}};
		islayers = {{iEL1I, "EL1I"}, {iEL1O, "EL1O"}, {iEL2I, "EL2I"}, {iEL2O, "EL2O"}, {iEL3I, "EL3I"}, {iEL3O, "EL3O"}, {iEL4I, "EL4I"}, {iEL4O, "EL4O"}};
	}
	if(side=="Pside")
	{
		layersnames = {"PL1I", "PL1O", "PL2I", "PL2O", "PL3I", "PL3O", "PL4I", "PL4O"};
		zlayer = {0, z1dipole, z2dipole, zPL1I, zPL1O, zPL2I, zPL2O, zPL3I, zPL3O, zPL4I, zPL4O};
		layersz = {zPL1I, zPL1O, zPL2I, zPL2O, zPL3I, zPL3O, zPL4I, zPL4O};
		layers = {{iPL1I, "PL1I"}, {iPL1O, "PL1O"}, {iPL2I, "PL2I"}, {iPL2O, "PL2O"}, {iPL3I, "PL3I"}, {iPL3O, "PL3O"}, {iPL4I, "PL4I"}, {iPL4O, "PL4O"}};
		szlayers = {{"PL1I", zPL1I}, {"PL1O", zPL1O}, {"PL2I", zPL2I}, {"PL2O", zPL2O}, {"PL3I", zPL3I}, {"PL3O", zPL3O}, {"PL4I", zPL4I}, {"PL4O", zPL4O}};
		silayers = {{"PL1I", iPL1I}, {"PL1O", iPL1O}, {"PL2I", iPL2I}, {"PL2O", iPL2O}, {"PL3I", iPL3I}, {"PL3O", iPL3O}, {"PL4I", iPL4I}, {"PL4O", iPL4O}};
		islayers = {{iPL1I, "PL1I"}, {iPL1O, "PL1O"}, {iPL2I, "PL2I"}, {iPL2O, "PL2O"}, {iPL3I, "PL3I"}, {iPL3O, "PL3O"}, {iPL4I, "PL4I"}, {iPL4O, "PL4O"}};
	}

	zLastLayer = (side=="Eside") ? zEL4I : zPL4I;
	zFirstLayer = (side=="Eside") ? zEL1O : zPL1O;
	LastLayer = (side=="Eside") ? "EL4I" : "PL4I";
	FirstLayer = (side=="Eside") ? "EL1O" : "PL1O";
	cout << "zLastLayer=" << zLastLayer << " (" << LastLayer << "), zFirstLayer=" << zFirstLayer << " (" << FirstLayer << ")" << endl;
	
	if(proc=="elaser" && side=="Pside") staveNumbers.insert(make_pair("iMin",0));
	if(proc=="elaser" && side=="Pside") staveNumbers.insert(make_pair("nMax",8));
	if(proc=="glaser" && side=="Pside") staveNumbers.insert(make_pair("iMin",0));
	if(proc=="glaser" && side=="Pside") staveNumbers.insert(make_pair("nMax",8));
	if(proc=="glaser" && side=="Eside") staveNumbers.insert(make_pair("iMin",8));
	if(proc=="glaser" && side=="Eside") staveNumbers.insert(make_pair("nMax",16));
	
	/// initialize the lookup table
	for (size_t i = 0; i<layersnames.size(); i++)
	{
		TString lname = layersnames.at(i);
		if(lname.Contains("P") && lname.Contains("I"))      axisMap.insert(make_pair(lname, new TAxis(nAxisBins, xMinPI, xMaxPI)));
		else if(lname.Contains("P") && lname.Contains("O")) axisMap.insert(make_pair(lname, new TAxis(nAxisBins, xMinPO, xMaxPO)));
		else if(lname.Contains("E") && lname.Contains("I")) axisMap.insert(make_pair(lname, new TAxis(nAxisBins, xMinEI, xMaxEI)));
		else if(lname.Contains("E") && lname.Contains("O")) axisMap.insert(make_pair(lname, new TAxis(nAxisBins, xMinEO, xMaxEO)));

		// If x is underflow or overflow, attempt to extend the axis if TAxis::kCanExtend is true. Otherwise, return 0 or fNbins+1.
		axisMap[lname]->SetCanExtend(0);
		vector<int> temp1;
		for (int b = 1; b<axisMap[lname]->GetNbins()+1; ++b)
		{
			lookupTable[lname].insert(make_pair(b, temp1));
		}
	}
	
	/// initialise the road width cuts
	if(side=="Pside")
	{
		fRWxPI = new TF2("fRWxPI","(0.06-0.01*x)/2 * (1+(abs(y-"+(TString)tostring(xMinPI)+")/("+(TString)tostring(Lstave)+"))^0.1)",1.,3.,xMinPI,xMaxPI); /// y is x4
		fRWxPO = new TF2("fRWxPO","(0.06-0.01*x)/2 * (1+(abs(y-"+(TString)tostring(xMinPO)+")/("+(TString)tostring(Lstave)+"))^0.1)",1.,3.,xMinPO,xMaxPO); /// y is x4
	}
	if(side=="Eside")
	{
		fRWxEI = new TF2("fRWxEI","(0.06-0.01*x)/2 * (1+(abs(y-"+(TString)tostring(xMaxEI)+")/("+(TString)tostring(Lstave)+"))^0.1)",1.,3.,xMinEI,xMaxEI); /// y is x4
		fRWxEO = new TF2("fRWxEO","(0.06-0.01*x)/2 * (1+(abs(y-"+(TString)tostring(xMaxEO)+")/("+(TString)tostring(Lstave)+"))^0.1)",1.,3.,xMinEO,xMaxEO); /// y is x4
	}
	fRWy = new TF2("fRWy","(0.06-0.01*x)/2 * (1+(abs(y)/("+(TString)tostring(Hstave)+"/2))^0.1)",1.,3.,yDn,yUp); /// y is y4
	
	cout << "++++++++++++++++++++++++++++++++++++" << endl;
}

unsigned int EncodeClusterId(unsigned int layer, unsigned int chip, unsigned int cellx, unsigned int celly)
{
	// layer should go into the first 4 bits [layer value ranges from 0 to 15]
	// chip should go into the 4 bits after that [chip value ranges from 0 to 8]
	// cellx should go into the 10 bits after that [cellx value ranges from 0 to 1023]
	// celly should go into the 9 bits after that [celly value ranges from 0 to 511]
	return (static_cast<unsigned int>(layer)<<23) + (static_cast<unsigned int>(chip)<<19) + (static_cast<unsigned int>(cellx)<<9) + (static_cast<unsigned int>(celly)<<0);
}
unsigned int DecodeLayer(const unsigned int clsid)
{
	// layer is the first 4 bits
	int check1 = 15;
	return (clsid>>23)&check1;
}
unsigned int DecodeChip(const unsigned int clsid)
{
	// chip is the 4 bits after that
	unsigned int check2 = 15;
	return (clsid>>19)&check2;
}
unsigned int DecodeCellX(const unsigned int clsid)
{
	// cellx is the 10 bits after that
	unsigned int check3 = 1023;
	return (clsid>>9)&check3;
}
unsigned int DecodeCellY(const unsigned int clsid)
{
	// celly is the last 9 bits
	unsigned int check4 = 511;
	return (clsid>>0)&check4;
}
unsigned int mapFullSim2KFLayer(TString side, int lyrid)
{  
	unsigned int KFLyr = -999;
	if(side=="Pside")
	{
		if(lyrid==0) KFLyr = silayers["PL1I"]; // 4;
		if(lyrid==1) KFLyr = silayers["PL1O"]; // 2;
		if(lyrid==2) KFLyr = silayers["PL2I"]; // 8;
		if(lyrid==3) KFLyr = silayers["PL2O"]; // 6;
		if(lyrid==4) KFLyr = silayers["PL3I"]; // 13;
		if(lyrid==5) KFLyr = silayers["PL3O"]; // 11;
		if(lyrid==6) KFLyr = silayers["PL4I"]; // 17;
		if(lyrid==7) KFLyr = silayers["PL4O"]; // 15;
	}
	else
	{
		if(lyrid==8)  KFLyr = silayers["EL1I"]; // 4;
		if(lyrid==9)  KFLyr = silayers["EL1O"]; // 2;
		if(lyrid==10) KFLyr = silayers["EL2I"]; // 8;
		if(lyrid==11) KFLyr = silayers["EL2O"]; // 6;
		if(lyrid==12) KFLyr = silayers["EL3I"]; // 13;
		if(lyrid==13) KFLyr = silayers["EL3O"]; // 11;
		if(lyrid==14) KFLyr = silayers["EL4I"]; // 17;
		if(lyrid==15) KFLyr = silayers["EL4O"]; // 15;
	}
	return KFLyr;
}


vector<TString> splitString(string s1)
{
    stringstream test(s1);
    string segment;
    vector<TString> seglist;
    while(getline(test, segment, '_'))
    {
        seglist.push_back(segment);
    }
    return seglist;
}
int getfiles(string path, TMapiTMapTSTS& fmap, TString side)
{
    int nevents = 0;
    for (const auto & entry : filesystem::directory_iterator(path))
    {
        string file = entry.path();
        // cout << file << endl;
        string barename = file.substr(file.find_last_of("/")+1);
        vector<TString> words = splitString(barename);
        int event = toint(words[words.size()-1].ReplaceAll(".root","").ReplaceAll("Event", ""));
        TString stave = words[words.size()-2];
        // cout << "event=" << event << ", stave=" << stave << endl;
		auto it = fmap.find(event);
		if(it==fmap.end())
		{
			TMapTSTS tempm;
			tempm[stave] = file;
			fmap.insert(make_pair(event,tempm));
		} 
		else fmap[event].insert(make_pair(stave, file));
    }
    for(TMapiTMapTSTS::iterator it1=fmap.begin(); it1!=fmap.end(); ++it1)
    {
        for(TMapTSTS::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2)
        {
            cout << "event: " << it1->first << " stave: " << it2->first << " filename: " << it2->second << endl;
        }
    }
	nevents = fmap.size();
	return nevents;
}


void setLogBins(Int_t nbins, Double_t min, Double_t max, Double_t* xpoints)
{
	Double_t logmin  = log10(min);
	Double_t logmax  = log10(max);
	Double_t logbinwidth = (Double_t)( (logmax-logmin)/(Double_t)nbins );
	xpoints[0] = min;
	for(Int_t i=1 ; i<=nbins ; i++) xpoints[i] = TMath::Power( 10,(logmin + i*logbinwidth) );
}


bool accept(double x, double y)
{
	bool failx = (x<xPsideL || (x>xPsideR && x<xEsideL) || x>xEsideR);
	bool faily = (y>yUp || y<yDn);
	if(failx || faily)
	return false;
	return true;
}

TVector2 rUnit2(TVector2 r1, TVector2 r2)
{
	TVector2 r = (r2-r1).Unit();
	return r;
}

float xofz(float *r1, float *r2, float z)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dz==0)
	{
		cout << "ERROR in xofz: dz=0" << endl;
		exit(-1);
	}
	float a = dx / dz;
	float b = r1[0]-a * r1[2];
	float x = a * z+b;
	// cout << "in xofz: x=" << x << ", dz=" << dz << endl;
	return x;
}

float xofz(float x1, float x2, float z1, float z2, float z)
{
	float r1[3] = {x1, 0, z1};
	float r2[3] = {x2, 0, z2};
	return xofz(r1, r2, z);
}

float yofz(float *r1, float *r2, float z)
{
	float dz = r2[2]-r1[2];
	float dy = r2[1]-r1[1];
	if(dz==0)
	{
		cout << "ERROR in yofz: dz=0" << endl;
		exit(-1);
	}
	float a = dy / dz;
	float b = r1[1]-a * r1[2];
	float y = a * z+b;
	return y;
}

float yofz(float y1, float y2, float z1, float z2, float z)
{
	float r1[3] = {0, y1, z1};
	float r2[3] = {0, y2, z2};
	return yofz(r1, r2, z);
}

float zofx(float *r1, float *r2, float x)
{
	float dz = r2[2]-r1[2];
	float dx = r2[0]-r1[0];
	if(dx==0)
	{
		cout << "ERROR in zofx: dx=0" << endl;
		exit(-1);
	}
	float a = dz / dx;
	float b = r1[2]-a * r1[0];
	float z = a * x+b;
	return z;
}

float zofx(float x1, float x2, float z1, float z2, float x)
{
	float r1[3] = {x1, 0, z1};
	float r2[3] = {x2, 0, z2};
	return zofx(r1, r2, x);
}

Color_t trkcol(double E)
{
	if(E >= 14)                return kBlack;
	else if(E<14. and E >= 12) return kRed;
	else if(E<12. and E >= 10) return 95;
	else if(E<10. and E >= 8.) return 91;
	else if(E<8. and E >= 7.)  return 80;
	else if(E<7. and E >= 6.)  return 71;
	else if(E<6. and E >= 5.)  return 65;
	else if(E<5. and E >= 4.)  return 60;
	else if(E<4. and E >= 3.)  return 53;
	else if(E<3. and E >= 2.)  return 51;
	else return 6;
	return kRed;
}

TLegend *trkcolleg()
{
	TLegend *leg = new TLegend(0.12, 0.30, 0.50, 0.60);
	leg->SetFillStyle(4000); // will be transparent
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->SetBorderSize(0);
	TLine *l1 = new TLine();
	l1->SetLineColor(trkcol(14));
	leg->AddEntry(l1, "#it{E}[GeV]#geq14", "l");
	TLine *l2 = new TLine();
	l2->SetLineColor(trkcol(12));
	leg->AddEntry(l2, "12#leq#it{E}[GeV]<14", "l");
	TLine *l3 = new TLine();
	l3->SetLineColor(trkcol(10));
	leg->AddEntry(l3, "10#leq#it{E}[GeV]<12", "l");
	TLine *l4 = new TLine();
	l4->SetLineColor(trkcol(8));
	leg->AddEntry(l4, "8#leq#it{E}[GeV]<10", "l");
	TLine *l5 = new TLine();
	l5->SetLineColor(trkcol(7));
	leg->AddEntry(l5, "7#leq#it{E}[GeV]<8", "l");
	TLine *l6 = new TLine();
	l6->SetLineColor(trkcol(6));
	leg->AddEntry(l6, "6#leq#it{E}[GeV]<7", "l");
	TLine *l7 = new TLine();
	l7->SetLineColor(trkcol(5));
	leg->AddEntry(l7, "5#leq#it{E}[GeV]<6", "l");
	TLine *l8 = new TLine();
	l8->SetLineColor(trkcol(4));
	leg->AddEntry(l8, "4#leq#it{E}[GeV]<5", "l");
	TLine *l9 = new TLine();
	l9->SetLineColor(trkcol(3));
	leg->AddEntry(l9, "3#leq#it{E}[GeV]<4", "l");
	TLine *l10 = new TLine();
	l10->SetLineColor(trkcol(2));
	leg->AddEntry(l10, "2#leq#it{E}[GeV]<3", "l");
	TLine *l11 = new TLine();
	l11->SetLineColor(trkcol(1));
	leg->AddEntry(l11, "#it{E}[GeV]<2", "l");
	return leg;
}

TPolyLine3D *TrackLine3d(const KMCProbeFwd *source, Double_t zMax, Double_t step = 1, Color_t col = kBlack)
{
	double xyz[3];
	source->GetXYZ(xyz);
	double zCurr = xyz[2]; // source->GetZ();
	int nZ = (zMax-zCurr) / step+1;
	if(nZ<2)
	{
		printf("bad limits\n");
		return 0;
	}
	KMCProbeFwd tmp(*source);
	double xp[nZ], yp[nZ], zp[nZ];
	xp[0] = xyz[0];
	yp[0] = xyz[1];
	zp[0] = xyz[2];
	int nz = 0;
	for (int iz = 1; iz<nZ; iz++)
	{
		if(!det->PropagateToZBxByBz(&tmp, TMath::Min(tmp.GetZ()+step, zMax), step)) break; // propagation may fail..
		tmp.GetXYZ(xyz);
		xp[iz] = xyz[0];
		yp[iz] = xyz[1];
		zp[iz] = xyz[2];
		nz++;
	}
	TPolyLine3D *polyline = new TPolyLine3D(nz+1);
	polyline->SetLineColor(col);
	for (int i = 0; i<nz+1; i++)
	{
		polyline->SetPoint(i, xp[i], yp[i], zp[i]);
	}
	return polyline;
}

bool islayer(double z, int layerindex = -1, double stepsize = 1)
{
	if(layerindex >= 0)
	{
		double dz = abs(zlayer[layerindex]-z);
		if(dz<stepsize) return true;
	}
	else
	{
		for (int j = 0; j<(int)zlayer.size(); ++j)
		{
			double dz = abs(zlayer[j]-z);
			if(dz<stepsize / 2.)
			{
				return true;
			}
		}
	}
	return false;
}

TPolyMarker3D *TrackMarker3d(const KMCProbeFwd *source, double zmin, double zmax, double zstep, Color_t col = kBlack)
{
	KMCProbeFwd tmp(*source);
	int nZ = (int)(zmax-zmin) / zstep;
	double xp[nZ], yp[nZ], zp[nZ];
	double xyz[3];
	tmp.GetXYZ(xyz);
	xp[0] = xyz[0];
	yp[0] = xyz[1];
	zp[0] = xyz[2];
	int nz = 0;
	for (int iz = 1; iz<nZ; iz++)
	{
		if(!det->PropagateToZBxByBz(&tmp, tmp.GetZ()+zstep, zstep)) break; // propagation may fail...
		tmp.GetXYZ(xyz);
		xp[iz] = xyz[0];
		yp[iz] = xyz[1];
		zp[iz] = xyz[2];
		nz++;
	}
	TPolyMarker3D *polymarker = new TPolyMarker3D(zlayer.size());
	polymarker->SetMarkerColor(col);
	int n = 0;
	for (int i = 0; i<nz+1; i++)
	{
		if(!islayer(zp[i])) continue;
		polymarker->SetPoint(n, xp[i], yp[i], zp[i]);
		n++;
	}
	return polymarker;
}

TPolyLine3D *GetLayer(TString side, TString io, double z, Color_t col)
{
	Int_t n = 5;

	Double_t xIE[] = {xMinEI, xMinEI, xMaxEI, xMaxEI, xMinEI};
	Double_t xIP[] = {xMinPI, xMinPI, xMaxPI, xMaxPI, xMinPI};

	Double_t xOE[] = {xMinEO, xMinEO, xMaxEO, xMaxEO, xMinEO};
	Double_t xOP[] = {xMinPO, xMinPO, xMaxPO, xMaxPO, xMinPO};

	Double_t y[] = {yDn, yUp, yUp, yDn, yDn};

	Double_t zC[] = {z, z, z, z, z};

	TPolyLine3D *polyline = 0;
	if(side=="P" && io=="I") polyline = new TPolyLine3D(n, xIP, y, zC);
	if(side=="E" && io=="I") polyline = new TPolyLine3D(n, xIE, y, zC);
	if(side=="P" && io=="O") polyline = new TPolyLine3D(n, xOP, y, zC);
	if(side=="E" && io=="O") polyline = new TPolyLine3D(n, xOE, y, zC);
	polyline->SetLineColor(col);
	return polyline;
}

TPolyLine3D *GetLayerFront(TString side, TString io, double z, Color_t col)
{
	Int_t n = 4;

	Double_t xIE[] = {xMinEI, xMinEI, xMaxEI, xMaxEI};
	Double_t xIP[] = {xMinPI, xMinPI, xMaxPI, xMaxPI};

	Double_t xOE[] = {xMinEO, xMinEO, xMaxEO, xMaxEO};
	Double_t xOP[] = {xMinPO, xMinPO, xMaxPO, xMaxPO};

	Double_t y[] = {yUp, yDn, yDn, yUp};

	Double_t zC[] = {z, z, z, z};

	TPolyLine3D *polyline = 0;
	if(side=="P" && io=="I") polyline = new TPolyLine3D(n, xIP, y, zC);
	if(side=="E" && io=="I") polyline = new TPolyLine3D(n, xIE, y, zC);
	if(side=="P" && io=="O") polyline = new TPolyLine3D(n, xOP, y, zC);
	if(side=="E" && io=="O") polyline = new TPolyLine3D(n, xOE, y, zC);
	polyline->SetLineColor(col);
	return polyline;
}

TPolyLine3D *GetDipole(Color_t col)
{
	TPolyLine3D *polyline = new TPolyLine3D();
	polyline->SetPoint(0, -xWdipole / 2, -yHdipole / 2, z1dipole);
	polyline->SetPoint(1, -xWdipole / 2, +yHdipole / 2, z1dipole);
	polyline->SetPoint(2, +xWdipole / 2, +yHdipole / 2, z1dipole);
	polyline->SetPoint(3, +xWdipole / 2, -yHdipole / 2, z1dipole);
	polyline->SetPoint(4, -xWdipole / 2, -yHdipole / 2, z1dipole);

	polyline->SetPoint(5, -xWdipole / 2, -yHdipole / 2, z2dipole); // go up

	polyline->SetPoint(6, -xWdipole / 2, +yHdipole / 2, z2dipole); // move
	polyline->SetPoint(7, -xWdipole / 2, +yHdipole / 2, z1dipole); // go down
	polyline->SetPoint(8, -xWdipole / 2, +yHdipole / 2, z2dipole); // up again

	polyline->SetPoint(9, +xWdipole / 2, +yHdipole / 2, z2dipole);	// move
	polyline->SetPoint(10, +xWdipole / 2, +yHdipole / 2, z1dipole); // go down
	polyline->SetPoint(11, +xWdipole / 2, +yHdipole / 2, z2dipole); // up again

	polyline->SetPoint(12, +xWdipole / 2, -yHdipole / 2, z2dipole); // move
	polyline->SetPoint(13, +xWdipole / 2, -yHdipole / 2, z1dipole); // go down
	polyline->SetPoint(14, +xWdipole / 2, -yHdipole / 2, z2dipole); // up again

	polyline->SetPoint(15, -xWdipole / 2, -yHdipole / 2, z2dipole); // move
	polyline->SetPoint(16, -xWdipole / 2, -yHdipole / 2, z1dipole); // go down
	polyline->SetPoint(17, -xWdipole / 2, -yHdipole / 2, z2dipole); // up again

	polyline->SetLineColor(col);
	return polyline;
}

TPolyLine3D *GetDipoleFront(Color_t col)
{
	TPolyLine3D *polyline = new TPolyLine3D();
	polyline->SetPoint(0, -xWdipole / 2, -yHdipole / 2, z1dipole);
	polyline->SetPoint(1, +xWdipole / 2, -yHdipole / 2, z1dipole);
	polyline->SetPoint(2, +xWdipole / 2, -yHdipole / 2, z2dipole);
	polyline->SetPoint(3, -xWdipole / 2, -yHdipole / 2, z2dipole);
	polyline->SetPoint(4, -xWdipole / 2, -yHdipole / 2, z1dipole);
	polyline->SetLineColor(col);
	return polyline;
}

bool skipglitches(TPolyMarker3D *points)
{
	Double_t x, y, z;
	for (int n = 0; n<points->GetN(); ++n)
	{
		points->GetPoint(n, x, y, z);
		if(abs(x)>40 || abs(y)>5) return true;
	}
	return false;
}

void WriteGeometry(vector<TPolyMarker3D *> &polm, vector<TPolyLine3D *> &poll, TString process, vector<int> &inacc, vector<TPolyMarker3D *> &clusters, TString suff = "")
{
	TCanvas *cnv_pl3d = new TCanvas("cnv_pl3d"+suff, "", 500, 500);
	TView *view_pl3d = TView::CreateView(1);
	view_pl3d->SetRange(-60, -20, 0, +60, +20, zLastLayer+15);
	view_pl3d->ShowAxis();

	TCanvas *cnv_pm3d = new TCanvas("cnv_pm3d"+suff, "", 500, 500);
	TView *view_pm3d = TView::CreateView(1);
	view_pm3d->SetRange(-60, -20, 0, +60, +20, zLastLayer+15);
	view_pm3d->ShowAxis();

	vector<TPolyLine3D *> staves;
	vector<TPolyLine3D *> fstaves;
	for (unsigned int l = 0; l<layersz.size(); ++l)
	{
		double z = layersz[l];
		TString io = (layersnames[l].Contains("I")) ? "I" : "O";
		TString pe = (layersnames[l].Contains("P")) ? "P" : "E";
		staves.push_back(GetLayer(pe, io, z, kGreen+3));
		fstaves.push_back(GetLayerFront(pe, io, z, kGreen+3));
	}

	TPolyLine3D *dipole = GetDipole(kGray);
	TPolyLine3D *fdipole = GetDipoleFront(kGray);

	cnv_pl3d->cd();
	dipole->Draw();
	for (unsigned int l = 0; l<staves.size(); l++) staves[l]->Draw();

	cnv_pm3d->cd();
	dipole->Draw();
	for (unsigned int l = 0; l<staves.size(); l++) staves[l]->Draw();

	for (int i = 0; i<(int)poll.size(); ++i)
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

	TLegend *leg = trkcolleg();
	cnv_pl3d->cd();
	fdipole->Draw();
	for (unsigned int l = 0; l<fstaves.size(); l++) fstaves[l]->Draw();
	leg->Draw("same");

	cnv_pm3d->cd();
	leg->Draw("same");

	cnv_pl3d->SaveAs(storage+"/output/root/"+process+"_tracks_pl3d"+suff+".root");
	cnv_pl3d->SaveAs(storage+"/output/pdf/"+process+"_tracks_pl3d"+suff+".pdf");
	cnv_pm3d->SaveAs(storage+"/output/root/"+process+"_tracks_pm3d"+suff+".root");
	cnv_pm3d->SaveAs(storage+"/output/pdf/"+process+"_tracks_pm3d"+suff+".pdf");

	
	TFile *flines = new TFile(storage+"/data/root/"+process+"_geometry"+suff+".root", "RECREATE");
	flines->cd();
	dipole->Write();
	fdipole->Write();
	for (unsigned int l = 0; l<staves.size(); l++)
	{
		staves[l]->Write();
		fstaves[l]->Write();
	}
	leg->Write();
	flines->Close();
}

bool accepttrk(vector<TPolyMarker3D *> &polm, int itrk)
{
	/// in acceptance?
	int nlayers = 4;
	int acctrk = 0;
	for (int i = 0; i<polm[itrk]->GetN(); i++)
	{
		Double_t xr, yr, zr;
		polm[itrk]->GetPoint(i, xr, yr, zr);
		if(zr<zLastLayer+15) continue; //// count only the active layers
		int inacclayer = accept(xr, yr);
		acctrk += inacclayer;
	}
	return (acctrk==nlayers);
}

bool foundinvec(int x, vector<int> &v)
{
	vector<int>::iterator it = find(v.begin(), v.end(), x);
	return (it != v.end());
}

int getvecindex(int x, vector<int> &v)
{
	vector<int>::iterator it = find(v.begin(), v.end(), x);
	return (it != v.end()) ? distance(v.begin(), it) : -1;
}

// removed const
bool matchTrkId(Cluster& lhs, Cluster& rhs)
{
	for(size_t i=0; i < lhs.trksid.size(); i++)
	{
		if(foundinvec(lhs.trksid[i], rhs.trksid)) return true;
	}
	return false; 
}

void prepare_cached_clusters()
{
	/// clear first
	for (TMapTSvCls::iterator it = cached_clusters.begin(); it != cached_clusters.end(); ++it) it->second.clear();
	cached_clusters.clear();
	
	/// then rebuild
	TMapii mapii;
	vector<Cluster> vc;
	for (TMapiTS::iterator it = layers.begin(); it != layers.begin(); ++it)
	{
		TString layername = it->second;
		cached_clusters.insert(make_pair(layername, vc));
		cached_clusters_id2ix.insert(make_pair(layername, mapii));
	}
}

void clear_cached_clusters()
{
	cached_clusters_id2lyr.clear();
	for (TMapTSvCls::iterator it = cached_clusters.begin(); it != cached_clusters.end(); ++it) it->second.clear();
	cached_clusters.clear();
	for (TMapTSMapii::iterator it = cached_clusters_id2ix.begin(); it != cached_clusters_id2ix.end(); ++it) it->second.clear();
	cached_clusters_id2ix.clear();
}
 
void cache_cluster(TVector3*               cls_r,
						 int                     cls_id,
						 int                     cls_isSig,
						 int                     cls_size,
						 int                     cls_sizex,
						 int                     cls_sizey,
						 double                  cls_charge,
						 vector<int>*            cls_trkids,
						 vector<int>*            cls_type,
						 vector<TLorentzVector>* cls_trksp,
						 TString side, TMapTSTH1D& histos1, TMapTSTH2D& histos2)
{
	double x = cls_r->X()/10.; // mm to cm
	double y = cls_r->Y()/10.; // mm to cm
	double z = cls_r->Z()/10.; // mm to cm
	int lyrid_FS = DecodeLayer(cls_id);
	// int chip_FS  = DecodeChip(cls_id);
	// int cellx_FS = DecodeCellX(cls_id);
	// int celly_FS = DecodeCellY(cls_id);
	// cout << "lyrid_FS: " << lyrid_FS << endl;
	int lyrid_KF = mapFullSim2KFLayer(side, lyrid_FS); 
	TString lyrname_KF = layers[lyrid_KF];
	// cout << "1. (x,y,z) " << x << ", " << y << ", " << z << ") lyrid_FS: " << lyrid_FS << " lyrid_KF: " << lyrid_KF << " lyrname_KF: " << lyrname_KF << endl; 

	x = (doMisalignmentX) ? x+XvariationSign * dxAlignmentXFEL : x;
	y = (doMisalignmentY) ? y+YvariationSign * dyAlignmentXFEL : y;

	Cluster cls;
	
	cls.lyridKF = lyrid_KF;
	cls.lyrnmKF = lyrname_KF;
	cls.clsid = cls_id;
	cls.npixl = cls_size;
	cls.npixlx = cls_sizex;
	cls.npixly = cls_sizey;
	cls.shape = -1;
	cls.charge = cls_charge;
	cls.r.SetXYZ(x, y, z);
	cls.issig = cls_isSig;
	int typeSummary = -999;
	for(size_t c=0; c<cls_type->size(); c++)
	{
		if(cls_type->at(c)==1) typeSummary=1;
		if(cls_type->at(c)==2 && typeSummary!=1) typeSummary=2;
		if(cls_type->at(c)==0 && typeSummary!=1 && typeSummary!=2) typeSummary=0;
		cls.trkstype.push_back(cls_type->at(c));
		cls.trksp.push_back(cls_trksp->at(c));
		cls.trksid.push_back(cls_trkids->at(c));
	}
	cls.type = typeSummary;
	// cout << "2. (x,y,z) = (" << x <<", " << y << ", " << z << ") lyrid_FS: " << lyrid_FS << " lyrname_KF: " << lyrname_KF << " clsid: " << cluster_id << " type: " << cluster_isSig << endl; 
	// cout << " chip_FS: " << chip_FS << ", cellX_FS: " << cellx_FS << ", cellY_FS: " << celly_FS << endl;
	
	/// fill occupancy plots
	if(lyrname_KF.Contains("L1I")) histos2["h_all_occ_L1I_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L2I")) histos2["h_all_occ_L2I_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L3I")) histos2["h_all_occ_L3I_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L4I")) histos2["h_all_occ_L4I_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L1O")) histos2["h_all_occ_L1O_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L2O")) histos2["h_all_occ_L2O_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L3O")) histos2["h_all_occ_L3O_"+side]->Fill(x,y);
	if(lyrname_KF.Contains("L4O")) histos2["h_all_occ_L4O_"+side]->Fill(x,y);
	
	/// fill cluster size plots
	if(lyrname_KF.Contains("L1I")) { histos1["h_all_csize_L1I_"+side]->Fill(cls_size); histos1["h_all_csizex_L1I_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L1I_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L2I")) { histos1["h_all_csize_L2I_"+side]->Fill(cls_size); histos1["h_all_csizex_L2I_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L2I_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L3I")) { histos1["h_all_csize_L3I_"+side]->Fill(cls_size); histos1["h_all_csizex_L3I_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L3I_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L4I")) { histos1["h_all_csize_L4I_"+side]->Fill(cls_size); histos1["h_all_csizex_L4I_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L4I_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L1O")) { histos1["h_all_csize_L1O_"+side]->Fill(cls_size); histos1["h_all_csizex_L1O_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L1O_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L2O")) { histos1["h_all_csize_L2O_"+side]->Fill(cls_size); histos1["h_all_csizex_L2O_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L2O_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L3O")) { histos1["h_all_csize_L3O_"+side]->Fill(cls_size); histos1["h_all_csizex_L3O_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L3O_"+side]->Fill(cls_sizey); }
	if(lyrname_KF.Contains("L4O")) { histos1["h_all_csize_L4O_"+side]->Fill(cls_size); histos1["h_all_csizex_L4O_"+side]->Fill(cls_sizex); histos1["h_all_csizey_L4O_"+side]->Fill(cls_sizey); }
	
	cached_clusters[lyrname_KF].push_back(cls);

	cached_clusters_id2lyr.insert(make_pair(cls_id, lyrid_KF));

	int index = cached_clusters[lyrname_KF].size()-1;
	cached_clusters_id2ix[lyrname_KF].insert(make_pair(cls_id, index));

	int bin = axisMap[lyrname_KF]->FindBin(x);

	lookupTable[lyrname_KF][bin].push_back(cls_id);
}

void reset_layers_all()
{
	for (Int_t l = 0; l<det->GetLayers()->GetEntries(); l++)
	{
		det->GetLayer(l)->ResetBgClusters();
		det->GetLayer(l)->ResetMCTracks();
		det->GetLayer(l)->Reset();
	}
}

void reset_layers_tracks(Int_t skip = -1)
{
	Int_t l0 = (skip >= 0) ? skip : 0;
	for (Int_t l = l0; l<det->GetLayers()->GetEntries(); l++)
	{
		det->GetLayer(l)->ResetMCTracks();
	}
}

/// clear the lookup table after each event
void clear_lookup_table()
{
	for (size_t i = 0; i<layersnames.size(); i++)
	{
		TString lname = layersnames.at(i);
		{
			for (int b = 1; b<axisMap[lname]->GetNbins()+1; ++b)
				lookupTable[lname][b].clear();
		}
	}
}

void remove_from_lookup_table(vector<Cluster> wincls)
{
   for(size_t i=0; i<wincls.size(); ++i)
   {
	   int clsid = wincls[i].clsid;
	   TString lyrnameKF = layers[wincls[i].lyridKF];
	   //int clsix = cached_clusters_id2ix[lyrnameKF][clsid];
	   double x = wincls[i].r.X();
	   int bin = axisMap[lyrnameKF]->FindBin(x);
	   int clsix = getvecindex(clsid, lookupTable[lyrnameKF][bin]);
	   if(clsix<0) continue;
	   lookupTable[lyrnameKF][bin].erase(lookupTable[lyrnameKF][bin].begin()+clsix);
   }
}

void embed_cluster(Cluster &cls)
{
	/// set the clusters of the seed
	double clxyzTrk[3];
	double clxyzLab[3];
	clxyzLab[0] = cls.r.X();
	clxyzLab[1] = cls.r.Y();
	clxyzLab[2] = cls.r.Z();
	KMCProbeFwd::Lab2Trk(clxyzLab, clxyzTrk);
	det->GetLayer(cls.lyridKF)->AddBgCluster(clxyzTrk[0], clxyzTrk[1], clxyzTrk[2], cls.clsid);
}

int embed_selective(TString side, TString lyrnumber, TString io, double xPivot, double x4, double y4, double z4, TMapTSTF1& fDy14vsYMap, TF1* fEvsXL4, vector<int> &embedded_clusters, double rwscale1=1, double rwscale2=1, double rwscale3=1)
{
	TString slyr  = side.ReplaceAll("side", "")+"L"+lyrnumber+io;
	TString slyr4 = side.ReplaceAll("side", "")+"L4"+io;
	double E = fEvsXL4->Eval(z4); // in GeV

	double rwx1 = (E>3) ? rwxL1*rwscale1 : rwxL1*rwscale1*2;
	double rwy1 = (E>3) ? rwyL1*rwscale1 : rwyL1*rwscale1*2;
	double rwx2 = (E>3) ? rwxL2*rwscale2 : rwxL2*rwscale2*2;
	double rwy2 = (E>3) ? rwyL2*rwscale2 : rwyL2*rwscale2*2;
	double rwx3 = (E>3) ? rwxL3*rwscale3 : rwxL3*rwscale3*2;
	double rwy3 = (E>3) ? rwyL3*rwscale3 : rwyL3*rwscale3*2;
	double rwxM = -1;
	double rwyM = -1;
	if(slyr.Contains("1")) { rwxM = rwx1; rwyM = rwy1; }
	if(slyr.Contains("2")) { rwxM = rwx2; rwyM = rwy2; }
	if(slyr.Contains("3")) { rwxM = rwx3; rwyM = rwy3; }
	int binUp = axisMap[slyr]->FindBin(xPivot+rwxM);
	int binDown = axisMap[slyr]->FindBin(xPivot-rwxM);

	int nembedded = 0;
	for (int bin = binDown; bin <= binUp; ++bin)
	{
		for (size_t k = 0; k<lookupTable[slyr][bin].size(); k++)
		{
			int clsid = lookupTable[slyr][bin][k];
			if(foundinvec(clsid, embedded_clusters)) continue;
			int index = cached_clusters_id2ix[slyr][clsid];
			double dycut = fDy14vsYMap[slyr4]->Eval(y4);
			double dy14 = y4-cached_clusters[slyr][index].r.Y();
			if(dy14>(dycut+rwyM)) continue;
			if(dy14<(dycut-rwyM)) continue;
			if(debug) cout << "embedding: slyr " << slyr << " index " << index << endl;
			embed_cluster(cached_clusters[slyr][index]);
			embedded_clusters.push_back(clsid);
			nembedded++;
		}
	}
	return nembedded;
}

/// turning on only those clusters which are in a specific bin along the rw
void add_all_clusters(TString side, TString slyr, int i4, TMapTSTF1& fDx14vsXMap, TMapTSTF1& fDy14vsYMap, TF1* fEvsXL4, vector<int>& embedded_clusters, double rwscale1=1, double rwscale2=1, double rwscale3=1)
{
	/// first embed the fourth layer pivot cluter
	// int ilyr = silayers[slyr];
	// int truix4 = cached_clusters[slyr][i4].clsid-ilyr*index_offset_sig;
	double x4 = cached_clusters[slyr][i4].r.X();
	double y4 = cached_clusters[slyr][i4].r.Y();
	double z4 = cached_clusters[slyr][i4].r.Z();
	if(!foundinvec(cached_clusters[slyr][i4].clsid, embedded_clusters))
	{
		embed_cluster(cached_clusters[slyr][i4]);
		embedded_clusters.push_back(cached_clusters[slyr][i4].clsid);
	}


	TString sname = (side=="Pside") ? "P" : "E";
	/// II: if x4 and x1 are in the inner layer
	double dxabs1I  = fDx14vsXMap[sname+"L4I"]->Eval(x4);
	double x1IPivot = (side=="Pside") ? (x4-dxabs1I) : (x4+dxabs1I);
	double z1I      = (side=="Pside") ? zPL1I : zEL1I; 


	/// OO: if x4 and x1 are in the outer layer
	double dxabs1O  = fDx14vsXMap[sname+"L4O"]->Eval(x4);
	double x1OPivot = (side=="Pside") ? (x4-dxabs1O) : (x4+dxabs1O);
	double z1O      = (side=="Pside") ? zPL1O : zEL1O; 


	/// OI: if x4 is in the outer layer and x1 is in the inner layer
	double dxabs1X  = fDx14vsXMap[sname+"L4X"]->Eval(x4);
	double x1XPivot = (side=="Pside") ? (x4-dxabs1X) : (x4+dxabs1X);
	double z1X      = (side=="Pside") ? zPL1I : zEL1I; 


	/// II: find the x along layer 2 and 3
	double x2IPivotII = xofz(x1IPivot, x4, z1I, z4, (side=="Pside") ? zPL2I : zEL2I);
	double x2OPivotII = xofz(x1IPivot, x4, z1I, z4, (side=="Pside") ? zPL2O : zEL2O);
	double x3IPivotII = xofz(x1IPivot, x4, z1I, z4, (side=="Pside") ? zPL3I : zEL3I);
	double x3OPivotII = xofz(x1IPivot, x4, z1I, z4, (side=="Pside") ? zPL3O : zEL3O);


	/// OO: find the x along layer 2 and 3
	double x2IPivotOO = xofz(x1OPivot, x4, z1O, z4, (side=="Pside") ? zPL2I : zEL2I);
	double x2OPivotOO = xofz(x1OPivot, x4, z1O, z4, (side=="Pside") ? zPL2O : zEL2O);
	double x3IPivotOO = xofz(x1OPivot, x4, z1O, z4, (side=="Pside") ? zPL3I : zEL3I);
	double x3OPivotOO = xofz(x1OPivot, x4, z1O, z4, (side=="Pside") ? zPL3O : zEL3O);


	/// OI: find the x along layer 2 and 3
	double x2IPivotOI = xofz(x1XPivot, x4, z1X, z4, (side=="Pside") ? zPL2I : zEL2I);
	double x2OPivotOI = xofz(x1XPivot, x4, z1X, z4, (side=="Pside") ? zPL2O : zEL2O);
	double x3IPivotOI = xofz(x1XPivot, x4, z1X, z4, (side=="Pside") ? zPL3I : zEL3I);
	double x3OPivotOI = xofz(x1XPivot, x4, z1X, z4, (side=="Pside") ? zPL3O : zEL3O);

	
	/// II: find the bins in layer 1 where the x values lie
	int n1I_II = embed_selective(side, "1", "I", x1IPivot,   x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2I_II = embed_selective(side, "2", "I", x2IPivotII, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2O_II = embed_selective(side, "2", "O", x2OPivotII, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3I_II = embed_selective(side, "3", "I", x3IPivotII, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3O_II = embed_selective(side, "3", "O", x3OPivotII, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);

	/// OO: find the bins in layer 1 where the x values lie
	int n1I_OO = embed_selective(side, "1", "O", x1OPivot,   x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2I_OO = embed_selective(side, "2", "I", x2IPivotOO, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2O_OO = embed_selective(side, "2", "O", x2OPivotOO, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3I_OO = embed_selective(side, "3", "I", x3IPivotOO, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3O_OO = embed_selective(side, "3", "O", x3OPivotOO, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);

	/// OI: find the bins in layer 1 where the x values lie
	int n1I_OI = embed_selective(side, "1", "I", x1IPivot,   x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2I_OI = embed_selective(side, "2", "I", x2IPivotOI, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n2O_OI = embed_selective(side, "2", "O", x2OPivotOI, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3I_OI = embed_selective(side, "3", "I", x3IPivotOI, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);
	int n3O_OI = embed_selective(side, "3", "O", x3OPivotOI, x4,y4,z4, fDy14vsYMap,fEvsXL4, embedded_clusters, rwscale1,rwscale2,rwscale3);


	if(0) cout << "n1I_II=" << n1I_II << endl;
	if(0) cout << "n2I_II=" << n2I_II << endl;
	if(0) cout << "n2O_II=" << n2O_II << endl;
	if(0) cout << "n3I_II=" << n3I_II << endl;
	if(0) cout << "n3O_II=" << n3O_II << endl;
	
	if(0) cout << "n1I_OO=" << n1I_OO << endl;
	if(0) cout << "n2I_OO=" << n2I_OO << endl;
	if(0) cout << "n2O_OO=" << n2O_OO << endl;
	if(0) cout << "n3I_OO=" << n3I_OO << endl;
	if(0) cout << "n3O_OO=" << n3O_OO << endl;

	if(0) cout << "n1I_OI=" << n1I_OI << endl;
	if(0) cout << "n2I_OI=" << n2I_OI << endl;
	if(0) cout << "n2O_OI=" << n2O_OI << endl;
	if(0) cout << "n3I_OI=" << n3I_OI << endl;
	if(0) cout << "n3O_OI=" << n3O_OI << endl;

	/// must sort clusters
	for (TMapiTS::iterator it = layers.begin(); it != layers.end(); ++it)
	{
		int lid = it->first;
		det->GetLayer(lid)->GetMCCluster()->Kill();
		det->GetLayer(lid)->SortBGClusters(); /// sort!!!
		/// after sorting, need to map the cluster ids to their indices!!!
		for (int n = 0; n<det->GetLayer(lid)->GetNBgClusters(); ++n)
		{
			int id = det->GetLayer(lid)->GetBgCluster(n)->GetTrID();
			cached_clusters_all_ids.insert(make_pair(id, n));
		}
	}
}

void print_all_clusters(TString side, bool doprint = true)
{
	if(!doprint) return;
	for (TMapiTS::iterator it = layers.begin(); it != layers.end(); ++it)
	{
		int ilr = it->first;
		TString slr = it->second;
		for (int c = 0; c<det->GetLayer(ilr)->GetNBgClusters(); c++)
		{
			KMCClusterFwd *cluster = det->GetLayer(ilr)->GetBgCluster(c);
			int id = cluster->GetTrID();
			float x = cluster->GetXLab();
			float y = cluster->GetYLab();
			float z = cluster->GetZLab();
			cout << "side=" << side << ", layer=" << ilr << ", id=" << id << " --> r={" << x << ", " << y << ", " << z << "}" << endl;
		}
		cout << endl;
	}
	// for (TMapii::iterator it = cached_clusters_all_ids.begin(); it != cached_clusters_all_ids.end(); it++)
	// {
	// 	cout << "id=" << it->first << " --> index=" << it->second << endl;
	// }
}



int makeseed_nonuniformB(TString process, float *r1, float *r4, TString side, TLorentzVector &p, TF1* fEvsXL4, TF1 *fDxvsX, TF1 *fDyvsY, bool doPrint=false)
{
	if(abs(r1[0])>=abs(r4[0])) return FAIL_SLOPE; // |x1| must be smaller than |x4|
	if(r1[0]*r4[0]<0)          return FAIL_SIDE; // not on the same side...
	if(r1[2]==r4[2])           return FAIL_SAMEZ; // trivial, make sure z is not the same
	float yDipoleExitAbsMax = (process=="glaser") ? 0.6  : 0.75; // cm
	float xDipoleExitAbsMin = (process=="glaser") ? 1.8  : 1.8; // cm
	float xDipoleExitAbsMax = (process=="glaser") ? 30.  : 30.;  // cm
	float yDipoleExit = yofz(r1, r4, zDipoleExit);
	float xDipoleExit = xofz(r1, r4, zDipoleExit);
	if(abs(yDipoleExit)>yDipoleExitAbsMax) return FAIL_DIPEXITYABS; // the track should point to |y|<yDipoleExitAbsMax at the dipole exit
	if(abs(xDipoleExit)<xDipoleExitAbsMin) return FAIL_DIPEXITXLOW; // the track should point to |x|>xDipoleExitAbsMin at the dipole exit
	if(abs(xDipoleExit)>xDipoleExitAbsMax) return FAIL_DIPEXITXHIGH; // the track should point to |x|<xDipoleExitAbsMax at the dipole exit
	float absdx41max = (process=="glaser") ? 8.0 : 8.0; //7.6; // cm, similar for elaser and glaser-derived from flat signal
	float absdx41min = (process=="glaser") ? 0.7 : 0.7; //0.8; // cm, similar for elaser and glaser-derived from flat signal
	double E = fEvsXL4->Eval(r4[0]); // in GeV
	double rwx = (E>3) ? rwxL1 : rwxL1*2;
	double rwy = (E>3) ? rwyL1 : rwyL1*2;
	if(abs(r4[0]-r1[0])>absdx41max || abs(r4[0]-r1[0])<absdx41min) return FAIL_DXSIZE; // new cut!!
	if(abs(r4[0]-r1[0])>(fDxvsX->Eval(r4[0])+rwx))                 return FAIL_RWXHIGH; // new cut!!
	if(abs(r4[0]-r1[0])<(fDxvsX->Eval(r4[0])-rwx))                 return FAIL_RWXLOW; // new cut!!
	if((r4[1]-r1[1])>(fDyvsY->Eval(r4[1])+rwy))                    return FAIL_RWYHIGH; // new cut!!
	if((r4[1]-r1[1])<(fDyvsY->Eval(r4[1])-rwy))                    return FAIL_RWYLOW; // new cut!!

	double P = sqrt(E*E-meGeV2);
	TVector2 v1(r1[2], r1[1]);
	TVector2 v4(r4[2], r4[1]);
	TVector2 u = rUnit2(v1, v4);
	double uz = u.X();
	double uy = u.Y();
	double px = 0;
	double py = P * uy;
	double pz = P * uz;
	
	p.SetPxPyPzE(px, py, pz, E);
	float EseedMin = (process=="glaser") ? EseedMinGLaser : EseedMinELaser; // GeV
	float EseedMax = (process=="glaser") ? EseedMaxGLaser : EseedMaxELaser; // GeV
	if(p.E()<EseedMin or p.E()>EseedMax) return FAIL_E;

	return PASS;
}

TString FormatEventID(int evnt)
{
	TString sevnt = "";
	if(evnt<10)                          sevnt = Form("000000%d", evnt);
	if(evnt >= 10 && evnt<100)           sevnt = Form("00000%d", evnt);
	if(evnt >= 100 && evnt<1000)         sevnt = Form("0000%d", evnt);
	if(evnt >= 1000 && evnt<10000)       sevnt = Form("000%d", evnt);
	if(evnt >= 10000 && evnt<100000)     sevnt = Form("00%d", evnt);
	if(evnt >= 100000 && evnt<1000000)   sevnt = Form("0%d", evnt);
	if(evnt >= 1000000 && evnt<10000000) sevnt = Form("%d", evnt); // assume no more than 9,999,999 events...
	return sevnt;
}




vector<int> nclusters_on_layer(vector<int>& embedded_clusters, TString side)
{
	vector<int> ncls_on_layer = {0,0,0,0};
 	for(int k=0 ; k<embedded_clusters.size() ; ++k)
	{
		unsigned int ilayer_FS = DecodeLayer(embedded_clusters[k]);
		unsigned int ilayer_KF = mapFullSim2KFLayer(side, ilayer_FS);
		if(islayers[ilayer_KF].Contains("1")) ncls_on_layer[0]++;
		if(islayers[ilayer_KF].Contains("2")) ncls_on_layer[1]++;
		if(islayers[ilayer_KF].Contains("3")) ncls_on_layer[2]++;
		if(islayers[ilayer_KF].Contains("4")) ncls_on_layer[3]++;
	}
	return ncls_on_layer;
}


int nmatched(int truid4, vector<Cluster>& wincls)
{
	int nmat = 0;
	for(size_t c=0; c<wincls.size(); ++c)
	{
		for(size_t i=0; i<wincls[c].trksid.size(); ++i)
		{
			if(wincls[c].trksid[i]==truid4) { nmat++; break; }
		}
	}
	return nmat;
}




int main(int argc, char *argv[])
{
	int argcounter;
	printf("Program Name Is: %s", argv[0]);
	if(argc >= 2)
	{
		printf("\nNumber Of Arguments Passed: %d", argc);
		printf("\n----Following Are The Command Line Arguments Passed----");
		for (argcounter = 0; argcounter<argc; argcounter++) printf("\nargv[%d]: %s", argcounter, argv[argcounter]);
		printf("\n");
	}
	//// minimum requirements
	if(argc<2) { printf("argc<2, exitting now\n"); exit(-1); }
	//// validate inputs
	if(argc==2 and !((TString)argv[1]).Contains("-proc=")) { printf("argc=2 but cannot parse %s\n", argv[1]); exit(-1); }
	if(argc==3 and !((TString)argv[2]).Contains("-smpl=")) { printf("argc=3 but cannot parse %s\n", argv[2]); exit(-1); }
	if(argc==4 and !((TString)argv[3]).Contains("-evnt=")) { printf("argc=4 but cannot parse %s\n", argv[3]); exit(-1); }
	if(argc==5 and !((TString)argv[4]).Contains("-seed=")) { printf("argc=5 but cannot parse %s\n", argv[4]); exit(-1); }
	if(argc==6 and !((TString)argv[5]).Contains("-ntrk=")) { printf("argc=6 but cannot parse %s\n", argv[5]); exit(-1); }
	//// assign inputs
	TString process  = ((TString)argv[1]).ReplaceAll("-proc=", "");	// mandatory
	TString smplname = ((TString)argv[2]).ReplaceAll("-smpl=", "");	// signal name
	int evnt     = (argc>3) ? toint(((TString)argv[3]).ReplaceAll("-evnt=", "")) : -1;	 // job id [optional]
	int Seed     = (argc>4) ? toint(((TString)argv[4]).ReplaceAll("-seed=", "")) : 12345;	 // seed [optional]
	int nsigtrks = (argc>5) ? toint(((TString)argv[5]).ReplaceAll("-ntrk=", "")) : -1; // job id [optional]
	
	
	// string path  = "/Volumes/Study/Weizmann_PostDoc/AllPix2Study/AllPixProceessedOutput/Signal_e0ppw_3.0"; /// this needs to be taken as an argument
	// string path  = "/Volumes/Study/Weizmann_PostDoc/AllPix2Study/AllPixProceessedOutput/Signal_e0ppw_3.0_Background"; /// this needs to be taken as an argument
	//string path  = "/Volumes/Study/Weizmann_PostDoc/AllPix2Study/AllPixProceessedOutput/ELaserBackground"; /// this needs to be taken as an argument
	//// print assigned inputs
	cout << "process=" << process << endl;
	cout << "smplname=" << smplname << endl;
	cout << "evnt=" << evnt << endl;
	cout << "nsigtrks=" << nsigtrks << endl;
	cout << "Seed=" << Seed << endl;
	
	TString digdir = storage+"/data/root/dig/"+process+"/"+smplname;
	
	// making the dir if it is not existing
	TString recdir = storage+"/data/root/rec/"+process+"/"+smplname;
	cout << "making dir (if not exists already): " << recdir << endl;
	gSystem->Exec("mkdir -p "+recdir);

	// TString proc = process;
	TString eventid = (evnt<0) ? "" : FormatEventID(evnt);
	TStopwatch stopwatch;
	TStopwatch stopwatch1;

	/// get the B-field vs xExit functions to read off the
	TString fFitsName = storage+"/output/root/inputs_for_reco_"+process+"_flat.root";
	TFile *fFits = new TFile(fFitsName, "READ");
	
	// TF1 *fEvsX_L1I_Eside = (TF1 *)fFits->Get("h2_E_vs_x_L1I_Eside");
	// TF1 *fEvsX_L1I_Pside = (TF1 *)fFits->Get("h2_E_vs_x_L1I_Pside");
	TF1* fEvsX_L4I_Eside = (TF1*)fFits->Get("h2_E_vs_x_L4I_Eside");
	TF1* fEvsX_L4I_Pside = (TF1*)fFits->Get("h2_E_vs_x_L4I_Pside");

	// TF1 *fEvsX_L1O_Eside = (TF1 *)fFits->Get("h2_E_vs_x_L1O_Eside");
	// TF1 *fEvsX_L1O_Pside = (TF1 *)fFits->Get("h2_E_vs_x_L1O_Pside");
	TF1* fEvsX_L4O_Eside = (TF1*)fFits->Get("h2_E_vs_x_L4O_Eside");
	TF1* fEvsX_L4O_Pside = (TF1*)fFits->Get("h2_E_vs_x_L4O_Pside");

	TF1 *fDx14vsX_L4I_Eside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4I_Eside");
	TF1 *fDx14vsX_L4I_Pside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4I_Pside");

	TF1 *fDx14vsX_L4O_Eside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4O_Eside");
	TF1 *fDx14vsX_L4O_Pside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4O_Pside");

	TF1 *fDx14vsX_L4X_Eside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4X_Eside");
	TF1 *fDx14vsX_L4X_Pside = (TF1 *)fFits->Get("h2_dx14_vs_x_L4X_Pside");
	
	
	TF1 *fDy14vsY_L4I_Eside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4I_Eside");
	TF1 *fDy14vsY_L4I_Pside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4I_Pside");
	
	TF1 *fDy14vsY_L4O_Eside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4O_Eside");
	TF1 *fDy14vsY_L4O_Pside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4O_Pside");
	
	TF1 *fDy14vsY_L4X_Eside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4X_Eside");
	TF1 *fDy14vsY_L4X_Pside = (TF1 *)fFits->Get("h2_dy14_vs_y_L4X_Pside");
	
  
	TMapTSTF1 fDx14vsXMap = {{"EL4I", fDx14vsX_L4I_Eside}, {"PL4I", fDx14vsX_L4I_Pside},
									 {"EL4O", fDx14vsX_L4O_Eside}, {"PL4O", fDx14vsX_L4O_Pside},
									 {"EL4X", fDx14vsX_L4X_Eside}, {"PL4X", fDx14vsX_L4X_Pside}};
	TMapTSTF1 fDy14vsYMap = {{"EL4I", fDy14vsY_L4I_Eside}, {"PL4I", fDy14vsY_L4I_Pside},
									 {"EL4O", fDy14vsY_L4O_Eside}, {"PL4O", fDy14vsY_L4O_Pside},
									 {"EL4X", fDy14vsY_L4X_Eside}, {"PL4X", fDy14vsY_L4X_Pside}};
					

	cout << "setup fits from files" << endl;

	int outN = (process=="elaser") ? 10 : 10;

	if(process=="elaser")
	{
		cout << "Doing only Pside!" << endl;
		sides.clear();
		sides.push_back("Pside"); /// do not reconstruct the Eside
	}

	// output tree
	cout << "Setting the output tree" << endl;
	gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyLine3D*>", "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
	gInterpreter->GenerateDictionary("vector<vector<TVector3> >", "vector");
	
	/// binning
	Int_t   nlogEbins0 = 10;
	Double_t  logEmin0 = 0.5;
	Double_t  logEmax0 = 16.5;
	Double_t logEbins0[nlogEbins0+1];
	setLogBins(nlogEbins0,logEmin0,logEmax0,logEbins0);
	Int_t   nlogEbins1 = 15;
	Double_t  logEmin1 = 0.5;
	Double_t  logEmax1 = 16.5;
	Double_t logEbins1[nlogEbins1+1];
	setLogBins(nlogEbins1,logEmin1,logEmax1,logEbins1);
	Int_t   nlogEbins2 = 40;
	Double_t  logEmin2 = 0.5;
	Double_t  logEmax2 = 16.5;
	Double_t logEbins2[nlogEbins2+1];
	setLogBins(nlogEbins2,logEmin2,logEmax2,logEbins2);
	Int_t   nlogEbins3 = 80;
	Double_t  logEmin3 = 0.5;
	Double_t  logEmax3 = 16.5;
	Double_t logEbins3[nlogEbins3+1];
	setLogBins(nlogEbins3,logEmin3,logEmax3,logEbins3);
	

	/////////////////////
	/// loop on the sides
	for (unsigned int s = 0; s<sides.size(); ++s)
	{
		TString side = sides[s];
		cout << "starting: " << side << endl;

		TFile *fOut = new TFile(recdir+"/rec_"+process+"_"+eventid+"_"+side+".root", "RECREATE");
		
		/// temporary variables
		vector<vector<int>> true_rec_imatch;
		vector<vector<int>> true_clusters_id;
		vector<float> reco_q;
		vector<TLorentzVector> reco_p;
		vector<float> reco_dErel;
		vector<float> reco_dpzrel;
		vector<float> reco_x;
		vector<float> reco_y;
		vector<float> reco_z;
		vector<vector<TVector3>> reco_trck_cls_r;
		vector<TPolyMarker3D *> reco_trckmar;
		vector<TPolyLine3D *> reco_trcklin;
		vector<float> reco_chi2dof;
		vector<int> reco_ismtchd;
		vector<int> reco_ixmtchd;
		vector<int> reco_idmtchd;
		vector<vector<int>> reco_clusters_id;

		/// setup the detector
		TString setup = "../setup/setupLUXE_"+process+"_"+side+"_FullSim.txt";
		det = new KMCDetectorFwd();
		det->ReadSetup(setup, setup);
		det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VertexTelescope
		// det->SetMinITSHits( det->GetNumberOfActiveLayersITS() ); // require hit in every layer
		det->SetMinITSHits(nMinHits+1); // require hit in at least 4 layers //TODO!!!
		det->SetMinMSHits(0);		  // we don't have muon spectrometer
		det->SetMinTRHits(0);		  // we don't have muon trigger stations
		// max number of seeds on each layer to propagate (per muon track)
		det->SetMaxSeedToPropagate(1000); // relevant only if background is considered
		// det->SetMaxSeedToPropagate(4000); // relevant only if background is considered
		// set chi2 cuts
		// det->SetMaxChi2Cl(10.);  // max track to cluster chi2
		// det->SetMaxChi2Cl(10.);  // max track to cluster chi2
		det->SetMaxChi2Cl(15.); // max track to cluster chi2
		// det->SetMaxChi2NDF(3.5); // max total chi2/ndf
		// det->SetMaxChi2NDF((process=="elaser")?15.:5.); // max total chi2/ndf
		// det->SetMaxChi2NDF((process=="elaser")?15.:5.); // max total chi2/ndf
		// det->SetMaxChi2NDF((process=="elaser") ? 15. : 5.); // max total chi2/ndf
		det->SetMaxChi2NDF(15.); // max total chi2/ndf
		det->SetMaxChi2Vtx(20e9); // fiducial cut on chi2 of convergence to vtx
		// det->SetMaxChi2Vtx(50.); // fiducial cut on chi2 of convergence to vtx
		// det->SetMaxChi2Vtx(1e3);  // fiducial cut on chi2 of convergence to vtx
		// det->SetMaxChi2Vtx(500);  // fiducial cut on chi2 of convergence to vtx
		// det->SetDefStepAir(1);				 // IMPORTANT FOR NON-UNIFORM FIELDS
		det->SetDefStepAir(1);				 // IMPORTANT FOR NON-UNIFORM FIELDS
		// det->SetDefStepMat(0.1);			 // NOAM??
		// det->SetMinP2Propagate(0.3);		 // GeV (NA60+)
		det->SetMinP2Propagate(0.1); /// GeV
		det->SetIncludeVertex(kTRUE);		 // count vertex as an extra measured point
		det->ImposeVertex(0., 0., 0.);		 // the vertex position is imposed NOAM
		// det->ImposeVertex(0.0014, 0., 0.);		 // the vertex position is imposed NOAM
		det->SetApplyBransonPCorrection(-1); // Branson correction, only relevant for setup with MS
		// for reconstruction:
		// det->SetErrorScale(500.);
		// det->SetErrorScale( (process=="elaser")?500.:200. );
		// det->SetErrorScale((process=="elaser") ? 500. : 500.); // was 400 earlier, can be also anywhere up to 1000
		det->SetErrorScale((process=="elaser") ? 500. : 500.); // was 400 earlier, can be also anywhere up to 1000
		det->Print();
		// det->BookControlHistos();

		////////////////////////////////////////
		setParametersFromDet(side, process); ///
		////////////////////////////////////////
		
		
		/// monitoring histograms
		TMapTSTH1D histos;
		TMapTSTH2D histos2;
		TString hname = "";
		
		/// for the histos
		double xMinI = (side=="Eside") ? xMinEI : xMinPI;
		double xMaxI = (side=="Eside") ? xMaxEI : xMaxPI;
		double xMinO = (side=="Eside") ? xMinEO : xMinPO;
		double xMaxO = (side=="Eside") ? xMaxEO : xMaxPO;
		
		/// bookeeping of number of BXs
		hname = "h_nBX_"+side; histos.insert(make_pair(hname, new TH1D()));
		histos[hname]->SetName(hname);
		histos[hname]->GetYaxis()->SetTitle("Entries");

		/// cutflow
		hname = "h_cutflow_"+side; histos.insert(make_pair(hname, new TH1D()));
		histos[hname]->SetName(hname);
		histos[hname]->GetYaxis()->SetTitle("Entries/BX");
		
		/// E/pz tru with various binning
		hname = "h_tru_pz_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];Tracks", 68, 0, 17)));
		hname = "h_tru_pz_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];Tracks", nlogEbins0,logEbins0)));
		hname = "h_tru_pz_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];Tracks", nlogEbins1,logEbins1)));
		hname = "h_tru_pz_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];Tracks", nlogEbins2,logEbins2)));
		hname = "h_tru_pz_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];Tracks", nlogEbins3,logEbins3)));		
		
		hname = "h_tru_E_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];Tracks", 68, 0, 17)));
		hname = "h_tru_E_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];Tracks", nlogEbins0,logEbins0)));
		hname = "h_tru_E_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];Tracks", nlogEbins1,logEbins1)));
		hname = "h_tru_E_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];Tracks", nlogEbins2,logEbins2)));
		hname = "h_tru_E_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];Tracks", nlogEbins3,logEbins3)));		
		
		/// inclusive clustering info
		hname = "h_all_csize_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L1I", 40, 0, 40)));
		hname = "h_all_csize_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L2I", 40, 0, 40)));
		hname = "h_all_csize_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L3I", 40, 0, 40)));
		hname = "h_all_csize_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L4I", 40, 0, 40)));
		hname = "h_all_csize_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L1O", 40, 0, 40)));
		hname = "h_all_csize_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L2O", 40, 0, 40)));
		hname = "h_all_csize_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L3O", 40, 0, 40)));
		hname = "h_all_csize_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;All clusters in L4O", 40, 0, 40)));

		hname = "h_all_csizex_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L1I", 20, 0, 20)));
		hname = "h_all_csizex_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L2I", 20, 0, 20)));
		hname = "h_all_csizex_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L3I", 20, 0, 20)));
		hname = "h_all_csizex_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L4I", 20, 0, 20)));
		hname = "h_all_csizex_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L1O", 20, 0, 20)));
		hname = "h_all_csizex_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L2O", 20, 0, 20)));
		hname = "h_all_csizex_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L3O", 20, 0, 20)));
		hname = "h_all_csizex_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;All clusters in L4O", 20, 0, 20)));

		hname = "h_all_csizey_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L1I", 20, 0, 20)));
		hname = "h_all_csizey_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L2I", 20, 0, 20)));
		hname = "h_all_csizey_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L3I", 20, 0, 20)));
		hname = "h_all_csizey_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L4I", 20, 0, 20)));
		hname = "h_all_csizey_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L1O", 20, 0, 20)));
		hname = "h_all_csizey_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L2O", 20, 0, 20)));
		hname = "h_all_csizey_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L3O", 20, 0, 20)));
		hname = "h_all_csizey_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;All clusters in L4O", 20, 0, 20)));
		
		/// 2D occupancy
		hname = "h_all_occ_L1I_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L1I occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinI,xMaxI, 100,yDn,yUp)));
		hname = "h_all_occ_L2I_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L2I occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinI,xMaxI, 100,yDn,yUp)));
		hname = "h_all_occ_L3I_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L3I occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinI,xMaxI, 100,yDn,yUp)));
		hname = "h_all_occ_L4I_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L4I occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinI,xMaxI, 100,yDn,yUp)));
		hname = "h_all_occ_L1O_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L1O occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinO,xMaxO, 100,yDn,yUp)));
		hname = "h_all_occ_L2O_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L2O occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinO,xMaxO, 100,yDn,yUp)));
		hname = "h_all_occ_L3O_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L3O occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinO,xMaxO, 100,yDn,yUp)));
		hname = "h_all_occ_L4O_"+side; histos2.insert(make_pair(hname, new TH2D(hname, "L4O occupancy per ~pixel per BX;x [cm];y [cm];Tracks/~pixel/BX", 900,xMinO,xMaxO, 100,yDn,yUp)));
		
		vector<TString> htypes = {"rec","mat","non","sel"};
		for(size_t h=0; h<htypes.size(); ++h)
		{
			TString htype = htypes[h];
			TString ytitle = "";
			TString ytitlerat = "";
			if(htype=="rec") { ytitle = "Reconstructed Tracks"; ytitlerat = "#frac{Reconstructed}{Truth} Tracks"; }
			if(htype=="sel") { ytitle = "Selected Tracks";      ytitlerat = "#frac{Selected}{Truth} Tracks";       }
			if(htype=="mat") { ytitle = "Matched Tracks";       ytitlerat = "#frac{Matched}{Truth} Tracks";        }
			if(htype=="non") { ytitle = "Unmatched Tracks";     ytitlerat = "#frac{Unmatched}{Truth} Tracks";      }
			cout << "booking " << htype << " histos" << endl;
			
			/// reconstructed clustering info
			if(htype=="rec")
			{
				hname = "h_"+htype+"_csize_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L1I", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L2I", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L3I", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L4I", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L1O", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L2O", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L3O", 40, 0, 40)));
				hname = "h_"+htype+"_csize_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels} in cluster;"+ytitle+" clusters in L4O", 40, 0, 40)));
         	
				hname = "h_"+htype+"_csizex_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L1I", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L2I", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L3I", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L4I", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L1O", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L2O", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L3O", 20, 0, 20)));
				hname = "h_"+htype+"_csizex_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{x}} in cluster;"+ytitle+" clusters in L4O", 20, 0, 20)));
         	
				hname = "h_"+htype+"_csizey_L1I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L1I", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L2I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L2I", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L3I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L3I", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L4I_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L4I", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L1O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L1O", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L2O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L2O", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L3O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L3O", 20, 0, 20)));
				hname = "h_"+htype+"_csizey_L4O_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{pixels}^{#it{y}} in cluster;"+ytitle+" clusters in L4O", 20, 0, 20)));
			}
			
			/// KF reconstructed outputs
			hname = "h_"+htype+"_Nhits_"+side;    histos.insert(make_pair(hname, new TH1D(hname, ";N_{hits} on track;"+ytitle, 5, 3, 8)));
			hname = "h_"+htype+"_chi2dof_"+side;  histos.insert(make_pair(hname, new TH1D(hname, ";#chi^{2}/N_{DoF};"+ytitle,150,0,15)));
			hname = "h_"+htype+"_SnpSig_"+side;   histos.insert(make_pair(hname, new TH1D(hname, ";Snp/#sigma(Snp);"+ytitle,  100,-50,+50)));
			hname = "h_"+htype+"_TglSig_"+side;   histos.insert(make_pair(hname, new TH1D(hname, ";Tgl/#sigma(Tgl);"+ytitle,  100,-1000,+1000)));
			hname = "h_"+htype+"_xVtxSig_"+side;  histos.insert(make_pair(hname, new TH1D(hname, ";#it{x}_{vtx}/#sigma(#it{x}_{vtx});"+ytitle,  150,-5e-6,+5e-6)));
			hname = "h_"+htype+"_yVtxSig_"+side;  histos.insert(make_pair(hname, new TH1D(hname, ";#it{y}_{vtx}/#sigma(#it{y}_{vtx});"+ytitle,  150,-0.0015,+0.0015)));
			hname = "h_"+htype+"_px_"+side;       histos.insert(make_pair(hname, new TH1D(hname,";#it{p}_{#it{x}} [GeV];"+ytitle, 200,-0.04,+0.04)));
			hname = "h_"+htype+"_px_zoom_"+side;  histos.insert(make_pair(hname, new TH1D(hname,";#it{p}_{#it{x}} [GeV];"+ytitle, 100,-0.02,+0.02)));
			hname = "h_"+htype+"_py_"+side;       histos.insert(make_pair(hname, new TH1D(hname,";#it{p}_{#it{y}} [GeV];"+ytitle, 200,-0.02,+0.02)));
			hname = "h_"+htype+"_py_zoom_"+side;  histos.insert(make_pair(hname, new TH1D(hname,";#it{p}_{#it{y}} [GeV];"+ytitle, 100,-0.01,+0.01)));

			/// pz reco with various binning
			hname = "h_"+htype+"_pz_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitle, 68, 0, 17)));
			hname = "h_"+htype+"_pz_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitle, nlogEbins0,logEbins0)));
			hname = "h_"+htype+"_pz_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitle, nlogEbins1,logEbins1)));
			hname = "h_"+htype+"_pz_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitle, nlogEbins2,logEbins2)));
			hname = "h_"+htype+"_pz_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitle, nlogEbins3,logEbins3)));

			/// E reco with various binning
			hname = "h_"+htype+"_E_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitle, 68, 0, 17)));
			hname = "h_"+htype+"_E_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitle, nlogEbins0,logEbins0)));
			hname = "h_"+htype+"_E_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitle, nlogEbins1,logEbins1)));
			hname = "h_"+htype+"_E_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitle, nlogEbins2,logEbins2)));
			hname = "h_"+htype+"_E_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitle, nlogEbins3,logEbins3)));
			
			/// E/pz rec/sel vs truth			
			hname = "h_resol_"+htype+"_dErel_"+side; histos.insert(make_pair(hname, new TH1D(hname, ";(E_{tru}-E_{rec})/E_{tru};"+ytitle, 100, -0.05, +0.05)));
			hname = "h_resol_"+htype+"_dpzrel_"+side; histos.insert(make_pair(hname, new TH1D(hname, ";(pz_{tru}-pz_{rec})/pz_{tru};"+ytitle, 100, -0.05, +0.05)));
			
			/// "efficiencies"
			hname = "h_ratio_"+htype+"_pz_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitlerat, 68, 0, 17)));
			hname = "h_ratio_"+htype+"_pz_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitlerat, nlogEbins0,logEbins0)));
			hname = "h_ratio_"+htype+"_pz_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitlerat, nlogEbins1,logEbins1)));
			hname = "h_ratio_"+htype+"_pz_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitlerat, nlogEbins2,logEbins2)));
			hname = "h_ratio_"+htype+"_pz_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{p}_{#it{z}} [GeV];"+ytitlerat, nlogEbins3,logEbins3)));
			
			hname = "h_ratio_"+htype+"_E_"+side;         histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitlerat, 68, 0, 17)));
			hname = "h_ratio_"+htype+"_E_"+side+"_log0"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitlerat, nlogEbins0,logEbins0)));
			hname = "h_ratio_"+htype+"_E_"+side+"_log1"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitlerat, nlogEbins1,logEbins1)));
			hname = "h_ratio_"+htype+"_E_"+side+"_log2"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitlerat, nlogEbins2,logEbins2)));
			hname = "h_ratio_"+htype+"_E_"+side+"_log3"; histos.insert(make_pair(hname, new TH1D(hname, ";#it{E} [GeV];"+ytitlerat, nlogEbins3,logEbins3)));
		}		
		cout << "histograms booked for " << side << endl;
		
		
		/// get the input files
		TMapiTMapTSTS fmap;
		int nevents = getfiles((string)digdir, fmap, side);
		histos["h_nBX_"+side]->Fill("N_{BX}",nevents);

		/// start event loop here
		int nEvntsProcessed = 0;
		int ibx = 0;
		cout << "Starting loop over nbxs=" << nevents << endl;
		for(TMapiTMapTSTS::iterator it1=fmap.begin(); it1!=fmap.end(); ++it1)
		{
			ibx++;
			/// prepare the dictionaries
			prepare_cached_clusters();

			/// for timing
			Double_t av_cputime = 0;
			Double_t av_realtime = 0;
			stopwatch.Start();
			nEvntsProcessed++;

			// /// clear output vectors: truth signal physics
			for (unsigned int x = 0; x<true_rec_imatch.size(); ++x) true_rec_imatch[x].clear();
			true_rec_imatch.clear();
			for (unsigned int x = 0; x<true_clusters_id.size(); ++x) true_clusters_id[x].clear();
			true_clusters_id.clear();
			
			// // /// clear output vectors: reconstruction
			reco_q.clear();
			reco_p.clear();
			reco_x.clear();
			reco_y.clear();
			reco_z.clear();
			reco_dErel.clear();
			reco_dpzrel.clear();
			for (unsigned int x = 0; x<reco_trck_cls_r.size(); ++x) reco_trck_cls_r[x].clear();
			for (unsigned int x = 0; x<reco_trckmar.size(); ++x) delete reco_trckmar[x];
			for (unsigned int x = 0; x<reco_trcklin.size(); ++x) delete reco_trcklin[x];
			reco_trck_cls_r.clear();
			reco_trckmar.clear();
			reco_trcklin.clear();
			reco_chi2dof.clear();
			reco_ismtchd.clear();
			reco_ixmtchd.clear();
			reco_idmtchd.clear();
			for (unsigned int x = 0; x<reco_clusters_id.size(); ++x) reco_clusters_id[x].clear();
			reco_clusters_id.clear();
			
			//// clear cached clusters
			clear_cached_clusters(); /// clear for both sides

			/// rest all the layers of the detector (including inactive if any)
			reset_layers_all(); // reset both sides

			///// clear the lookup table for each event first
			clear_lookup_table();

			/// clear this side's indices
			cached_clusters_all_ids.clear();
			// cout << "done clearing up" << endl;

			/// set the charge
			float crg = (side=="Eside") ? -1 : +1;

			/// globals (per side)
			unsigned int n_truth = 0;
			unsigned int n_seeds = 0;
			unsigned int n_solve = 0;
			unsigned int n_recos = 0;
			unsigned int n_selct = 0;
			unsigned int n_match = 0;
			// unsigned int n_trumt = 0;
			vector<int> vtruid;

			// this is the loop on chips belonging to the same bx
			// cout << "Starting loop over " << it1->second.size() << " available chips from 72(=9*4*2) chips in " << side << endl;
			for(TMapTSTS::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2)
			{
				// this is the loop on the files belonging to the same stave and same event
				//cout << it1->first << " " << it2->first << " " << it2->second << endl;

				TString inFileName = it2->second;
				TFile* fIn = new TFile(inFileName, "READ");
				TTree* tIn = (TTree *)fIn->Get("clusters");
				if(!tIn) continue;

				/// for cluster tree
				int bx;
				vector<int>* pix_entry=0; //<— this is the entry of the corresponding pixels in the pixels tree
				vector<int>* tru_pdgId=0;
				vector<int>* tru_trackId=0;
				vector<int>* pixId=0;
				vector<int>* tru_type=0;
				vector<int>* tru_pixentry=0; // split the line
				vector<double>* tru_edep=0;
				vector<TLorentzVector>* tru_p=0;
				vector<TVector3>* tru_hit=0; //<— x,y,z position in the pixel 
				vector<TVector3>* tru_vertex=0;
				int encodedClsId;
				int isSignal;
				int cellx_cog;
				int celly_cog;
				int cellx_geo;
				int celly_geo;
				int size;
				int xsize;
				int ysize;
				double charge;
				TVector3* rglobal_cog=0;
				TVector3* rglobal_geo=0;
				TVector3* rlocal_cog=0;
				TVector3* rlocal_geo=0;
				vector<double>* xres=0;
				vector<double>* yres=0;
				int seed_entry;

				tIn->SetBranchAddress("bx", &bx);
				tIn->SetBranchAddress("pix_entry", &pix_entry);
				tIn->SetBranchAddress("pixId", &pixId);
				tIn->SetBranchAddress("encodedClsId", &encodedClsId);
				tIn->SetBranchAddress("isSignal", &isSignal);
				tIn->SetBranchAddress("size", &size);
				tIn->SetBranchAddress("xsize", &xsize);
				tIn->SetBranchAddress("ysize", &ysize);
				tIn->SetBranchAddress("charge", &charge);
				tIn->SetBranchAddress("rglobal_cog", &rglobal_cog);
				tIn->SetBranchAddress("rlocal_cog", &rlocal_cog);
				tIn->SetBranchAddress("rglobal_geo", &rglobal_geo);
				tIn->SetBranchAddress("rlocal_geo", &rlocal_geo);
				tIn->SetBranchAddress("cellx_cog", &cellx_cog);
				tIn->SetBranchAddress("celly_cog", &celly_cog);
				tIn->SetBranchAddress("cellx_geo", &cellx_geo);
				tIn->SetBranchAddress("celly_geo", &celly_geo);
				tIn->SetBranchAddress("seed_entry", &seed_entry);
				tIn->SetBranchAddress("tru_type", &tru_type);
				tIn->SetBranchAddress("tru_pixentry", &tru_pixentry);
				tIn->SetBranchAddress("tru_pdgId", &tru_pdgId);
				tIn->SetBranchAddress("tru_trackId", &tru_trackId);
				tIn->SetBranchAddress("tru_edep", &tru_edep);
				if(smplname=="bkg") tIn->SetBranchAddress("true_p", &tru_p);
				else                tIn->SetBranchAddress("tru_p", &tru_p);
				tIn->SetBranchAddress("tru_hit", &tru_hit);
				tIn->SetBranchAddress("tru_vertex", &tru_vertex);
				tIn->SetBranchAddress("xres", &xres);
				tIn->SetBranchAddress("yres", &yres);


				/// make a pool of all signal clusters
				/// loop over the clusters
				
				Int_t ntotclusnmbr = tIn->GetEntries();
				// cout << "Starting loop over " << ntotclusnmbr << " clusters of chip " << it2->first << endl;
				for(size_t clsnmbr=0; clsnmbr < ntotclusnmbr; clsnmbr++)
				{	
					tIn->GetEntry(clsnmbr);
					// /// fill truth signal tracks:
					for(size_t t=0; t<tru_p->size(); ++t)
					{
						if(tru_type->at(t)!=1) continue;
						if(!foundinvec(tru_trackId->at(t), vtruid))
						{
							histos["h_tru_pz_"+side]->Fill(tru_p->at(t).Pz());
							histos["h_tru_pz_"+side+"_log0"]->Fill(tru_p->at(t).Pz());
							histos["h_tru_pz_"+side+"_log1"]->Fill(tru_p->at(t).Pz());
							histos["h_tru_pz_"+side+"_log2"]->Fill(tru_p->at(t).Pz());
							histos["h_tru_pz_"+side+"_log3"]->Fill(tru_p->at(t).Pz());
							
							histos["h_tru_E_"+side]->Fill(tru_p->at(t).E());
							histos["h_tru_E_"+side+"_log0"]->Fill(tru_p->at(t).E());
							histos["h_tru_E_"+side+"_log1"]->Fill(tru_p->at(t).E());
							histos["h_tru_E_"+side+"_log2"]->Fill(tru_p->at(t).E());
							histos["h_tru_E_"+side+"_log3"]->Fill(tru_p->at(t).E());
							vtruid.push_back(tru_trackId->at(t));
							n_truth++;
						}
					}
					/// caching is happening here 
					cache_cluster(rglobal_geo, encodedClsId, isSignal, size, xsize, ysize, charge, tru_trackId, tru_type, tru_p, side, histos, histos2);
				} // end of loop on clusters
				fIn->Close(); // must close the file
			} // end loop on staves
			
			
			/// fill cutflow in the beginning of the cutflow
			histos["h_cutflow_"+side]->Fill("Truth", n_truth);

			/// run over all clusters of layer 4 in the pool --> these are the seeds for the KalmanFilter fit
			TString slyr4I = (side=="Eside") ? "EL4I" : "PL4I";
			TString slyr1I = (side=="Eside") ? "EL1I" : "PL1I";
			TString slyr4O = (side=="Eside") ? "EL4O" : "PL4O";
			TString slyr1O = (side=="Eside") ? "EL1O" : "PL1O";

			/// loop on seeds
			unsigned int n1I = cached_clusters[slyr1I].size();
			unsigned int n1O = cached_clusters[slyr1O].size();
			unsigned int n4I = cached_clusters[slyr4I].size();
			unsigned int n4O = cached_clusters[slyr4O].size();
			int n4count = 1;
			cout << "Starting loop over layer 4 clusters with " << (n4I+n4O) << " clusters" << endl;
			histos["h_cutflow_"+side]->Fill("L4 Clusters", (n4I+n4O));
			for (unsigned int i4all = 0; i4all<(n4I+n4O); ++i4all)
			{
				det->SetErrorScale((process=="elaser") ? 500 : 500);
				det->SetMaxChi2NDF(15.);
				det->SetMaxChi2Cl(15.);
				
				/// outer first
				unsigned int i4 = (i4all<n4O) ? i4all  : i4all-n4O;
				TString slyr4   = (i4all<n4O) ? slyr4O : slyr4I;
				// int ilyr4       = (i4all<n4O) ? ilyr4O : ilyr4I;
				if(slyr4==slyr4O && (side=="Pside" && cached_clusters[slyr4][i4].r.X()<xMaxPI)) continue;
				if(slyr4==slyr4O && (side=="Eside" && cached_clusters[slyr4][i4].r.X()>xMinEI)) continue;

				// /// inner first
				// unsigned int i4 = (i4all<n4I) ? i4all  : i4all-n4I;
				// TString slyr4   = (i4all<n4I) ? slyr4I : slyr4O;
				// // int ilyr4       = (i4all<n4I) ? ilyr4I : ilyr4O;
				// if(slyr4==slyr4O && (side=="Pside" && cached_clusters[slyr4][i4].r.X()<xMaxPI)) continue;
				// if(slyr4==slyr4O && (side=="Eside" && cached_clusters[slyr4][i4].r.X()>xMinEI)) continue;
				
				
				/// cluster id for the pivot 
				int clsid4 = cached_clusters[slyr4][i4].clsid;
				
				
				/// is it a cluster associated with a signal track?
				int itru     = -999;
				double Etru  = -999;
				double pztru = -999;
				bool issig4  = false;
				if(cached_clusters[slyr4][i4].issig==1)
				{
					for(size_t t=0; t<cached_clusters[slyr4][i4].trksid.size(); ++t)
					{
						if(cached_clusters[slyr4][i4].trkstype[t]!=1) continue; // not a signal track
						itru  = cached_clusters[slyr4][i4].trksid[t];
						Etru  = cached_clusters[slyr4][i4].trksp[t].E();
						pztru = cached_clusters[slyr4][i4].trksp[t].Pz();
						issig4 = true;
						break;
					}
				}

				/// rest all the layers of the detector (including inactive if any)
				reset_layers_all(); // reset both sides
				// reset all tracks from all layers
				reset_layers_tracks();
				/// clear this side's indices
				cached_clusters_all_ids.clear();

				TF1 *fEvsXL4 = 0;
				TF1 *fDx14vsX = 0;
				TF1 *fDy14vsY = 0;
				
				if(slyr4.Contains("I")) fEvsXL4 = (side=="Pside") ? fEvsX_L4I_Pside : fEvsX_L4I_Eside;
				else                    fEvsXL4 = (side=="Pside") ? fEvsX_L4O_Pside : fEvsX_L4O_Eside;


				/// add all clusters to the detector
				vector<int> embedded_clusters;
				int n1inroad,n2inroad,n3inroad;
				double rwscl1 = 1;
				double rwscl2 = 1;
				double rwscl3 = 1;
				int niterations_nMax = 4;
				int niterations_n = 1;
				
				goto embed;
				
				embed:
					add_all_clusters(side,slyr4,i4,fDx14vsXMap,fDy14vsYMap,fEvsXL4, embedded_clusters, rwscl1,rwscl2,rwscl3);
					vector<int> ncls_on_layer = nclusters_on_layer(embedded_clusters, side);
					n1inroad = ncls_on_layer[0];
					n2inroad = ncls_on_layer[1];
					n3inroad = ncls_on_layer[2];
				
				if(niterations_n<niterations_nMax && (n1inroad<1 || n2inroad<1 || n3inroad<1))
				{
					if(n1inroad<1) rwscl1 *= 1.5;
					if(n2inroad<1) rwscl2 *= 2.0;
					if(n3inroad<1) rwscl3 *= 2.5;
					niterations_n++;
					goto embed;
				}
				// print_all_clusters(side);
				if(n1inroad<1 || n2inroad<1 || n3inroad<1)
				{
					if(itru!=-999) cout << "Insufficient hits n123?=("<<(n1inroad>0)<<","<<(n2inroad>0)<<","<<(n3inroad>0)<<") after " << niterations_n << " iterations for " << "itru=" << itru << " (Etru=" << Etru << ")" << endl;
					continue;
				}
				histos["h_cutflow_"+side]->Fill("Sufficient hits",1);
				
				/// find the momentum of the seeds
				TLorentzVector pseed;
				vector<TLorentzVector> pseeds;
				float r1[3];
				float r4[3];
				r4[0] = cached_clusters[slyr4][i4].r.X();
				r4[1] = cached_clusters[slyr4][i4].r.Y();
				r4[2] = cached_clusters[slyr4][i4].r.Z();
				TString sname = (side=="Pside") ? "P" : "E";
				double xMinI = (side=="Eside") ? xMinEI : xMinPI;
				double xMaxI = (side=="Eside") ? xMaxEI : xMaxPI;
				double xMinO = (side=="Eside") ? xMinEO : xMinPO;
				double xMaxO = (side=="Eside") ? xMaxEO : xMaxPO;
				
				/// II: if x4 and x1 are in the inner layer
				if(r4[0]>xMinI && r4[0]<xMaxI)
				{
					fDx14vsX = fDx14vsXMap[sname+"L4I"];
					fDy14vsY = fDy14vsYMap[sname+"L4I"];
					double dxabs = fDx14vsX->Eval(r4[0]);
					double dy = fDy14vsY->Eval(r4[1]);
					r1[0] = (side=="Pside") ? (r4[0]-dxabs) : (r4[0]+dxabs);
					r1[1] = r4[1]-dy;
					r1[2] = (side=="Pside") ? zPL1I : zEL1I;
					if(r1[0]>xMinI && r1[0]<xMaxI)
					{
						int seedflag = makeseed_nonuniformB(process,r1,r4,side,pseed,fEvsXL4,fDx14vsX,fDy14vsY);
						if(seedflag==PASS) pseeds.push_back(pseed); 
						else
						{
							if(itru!=-999) cout << "II seed fail due to: " << seedcutnames[seedflag] << endl;
						}
					}
				}
				/// OO: if x4 and x1 are in the outer layer
				if(r4[0]>xMinO && r4[0]<xMaxO)
				{
					fDx14vsX = fDx14vsXMap[sname+"L4O"];
					fDy14vsY = fDy14vsYMap[sname+"L4O"];
					double dxabs = fDx14vsX->Eval(r4[0]);
					double dy = fDy14vsY->Eval(r4[1]);
					r1[0] = (side=="Pside") ? (r4[0]-dxabs) : (r4[0]+dxabs);
					r1[1] = r4[1]-dy;
					r1[2] = (side=="Pside") ? zPL1O : zEL1O;
					if(r1[0]>xMinO && r1[0]<xMaxO)
					{
						int seedflag = makeseed_nonuniformB(process,r1,r4,side,pseed,fEvsXL4,fDx14vsX,fDy14vsY);
						if(seedflag==PASS) pseeds.push_back(pseed);
						else
						{
							if(itru!=-999) cout << "OO seed fail due to: " << seedcutnames[seedflag] << endl;
						}
					}
				}
				/// OI: if x4 is in the outer layer and x1 is in the inner layer
				if(r4[0]>xMinO && r4[0]<xMaxO)
				{
					fDx14vsX = fDx14vsXMap[sname+"L4X"];
					fDy14vsY = fDy14vsYMap[sname+"L4X"];
					double dxabs = fDx14vsX->Eval(r4[0]);
					double dy    = fDy14vsY->Eval(r4[1]);
					r1[0] = (side=="Pside") ? (r4[0]-dxabs) : (r4[0]+dxabs);
					r1[1] = r4[1]-dy;
					r1[2] = (side=="Pside") ? zPL1I : zEL1I;
					if(r1[0]>xMinI && r1[0]<xMaxI)
					{
						int seedflag = makeseed_nonuniformB(process,r1,r4,side,pseed,fEvsXL4,fDx14vsX,fDy14vsY);
						if(seedflag==PASS) pseeds.push_back(pseed);
						else
						{
							if(itru!=-999) cout << "OI seed fail due to: " << seedcutnames[seedflag] << endl;
						}
					}
				}
				if(pseeds.size()<1)
				{
					if(itru!=-999) cout << "clusterid=" << clsid4 << " (itru=" << itru << ", Etru=" << Etru << ")" << " has no seeds!" << endl;
					continue;
				}
				else
				{
					n_seeds++;
					histos["h_cutflow_"+side]->Fill("Seeds",1);
				}
				
				/// instant summary
				int countmateff = (int)((float)n_match / (float)n4count * 100.);
				int countreceff = (int)((float)n_recos / (float)n4count * 100.);
				int countseleff = (int)((float)n_selct / (float)n4count * 100.);
				if(n4count%100==0)
				{
					cout << "ibx=" << ibx
							<< ": " << n_seeds
								<< " seeds for " << n4count 
									<< "/" << (n4I+n4O)
										<< " clusters in " << slyr4 << " inc overlap (with " 
											<< n_recos << " recos with " 
												<< n_selct << " selected and " 
													<< n_match << " matched) -> counting: " 
														<< countreceff << "%(rec), "  
															<< countseleff << "%(sel), "  
																<< countmateff << "%(mat)" << endl;
				}
				n4count++;
				
				bool doPrint = false;
				// prepare the probe from the seed and do the KF fit


				stopwatch1.Start();
				
				int nMaxIterations = 3;
				int nIterations_slv = 1;
				int nIterations_trw = 1;
				int nIterations_hit = 1;
				int nIterations_kil = 1;
				bool solved = false;
				bool wasSolved = false;

				/// reconstruction!!!
				goto reco;

				reco:
					solved = det->SolveSingleTrackViaKalmanMC_Noam_multiseed(pseeds, meGeV, crg, 99, doPrint);

				/// solve the track with the KF
				if(!solved)
				{
					if(nIterations_slv<nMaxIterations)
					{
						nIterations_slv++;
						goto reco;
					}
					continue; // reconstruction failed
				}
				if(!wasSolved)
				{
					n_solve++;
					histos["h_cutflow_"+side]->Fill("KF Solved",1);
					wasSolved = true;
				}

				
				// get the reconstructed propagated to the vertex
				KMCProbeFwd *trw = det->GetLayer(0)->GetWinnerMCTrack();
				if(!trw)
				{
					if(nIterations_trw<nMaxIterations)
					{
						nIterations_trw++;
						double err = 500;
						det->SetErrorScale((process=="elaser") ? err*2*nIterations_trw : err*2*nIterations_trw);
						det->SetMaxChi2NDF(15.+(10*nIterations_trw));
						det->SetMaxChi2Cl(15.+(10*nIterations_trw));
						goto reco;
					}
					if(itru!=-999) cout << "!trw: clusterid=" << clsid4 << " (itru=" << itru << ", E=" << Etru << ", iteration=" << nIterations_trw << " out of " << nMaxIterations << ")" << endl;
					continue; // track was not reconstructed
				}
				
				/// kill the tracks with less than minimum hits-this should be killed on the KF side
				if(trw->GetNITSHits()<nMinHits)
				{
					if(nIterations_hit<nMaxIterations)
					{
						nIterations_hit++;
						goto reco;
					}
					if(itru!=-999) cout << "trw->GetNITSHits()<nMinHits: clusterid=" << clsid4 << " (E=" << Etru << ")" << endl;
					trw->Kill(); // track has too few hits
					continue; // too few hits
				}
				
				/// if the track is killed
				if(trw->IsKilled())
				{
					if(nIterations_kil<nMaxIterations)
					{
						nIterations_kil++;
						goto reco;
					}
					if(itru!=-999) cout << "trw->IsKilled(): clusterid=" << clsid4 << " (E=" << Etru << ")" << endl;
					continue; // track was killed
				}
				n_recos++;
				histos["h_cutflow_"+side]->Fill("KF Tracks",1);
				
				stopwatch1.Stop();
				// Double_t cputime1 = stopwatch1.CpuTime();
				// Double_t realtime1 = stopwatch1.RealTime();
				// cout << "cputime1=" << cputime1 << ", realtime1=" << realtime1 << endl;
				// if(doPrint)
				// {
				// 	cout << "Track fit succeeded, associated clusters are:" << endl;
				// 	trw->Print("clid");
				// }

				/// save clusters of the track
				vector<TVector3> v3tmp;
				reco_trck_cls_r.push_back(v3tmp);

				/// the rec index-important to set it here!
				unsigned int irec = reco_trck_cls_r.size()-1;

				/// get the clusters of the winner tracK
				vector<int> win_cls_id;
				vector<int> win_cls_inx;
				TMapii win_cls_id2lr;
				const int *probeclusters = trw->GetClID();
				int nprobeclusters = sizeof(probeclusters);
				// cout << "nprobeclusters=" << nprobeclusters << endl;
				// trw->Print("clid etp");
				// trw->Print("clid");
				vector<Cluster> wincls;
				for (int l = 0; l <= nprobeclusters; ++l)
				{
					int cid = probeclusters[l];
					if(cid<0) continue;
					int cix = cached_clusters_all_ids[cid];
					int ilr = cached_clusters_id2lyr[cid]; /// essential to find the correct layer!
					TString lname = layers[ilr];
					int index_in_cached_clusters = cached_clusters_id2ix[lname][cid];
					wincls.push_back(cached_clusters[lname][index_in_cached_clusters]);
					win_cls_id.push_back(cid); // provide active layer ID, not the physical ones (most are passive)
					win_cls_inx.push_back(cix);
					win_cls_id2lr.insert(make_pair(cid, ilr));
					
					double xwin = det->GetLayer(ilr)->GetBgCluster(cix)->GetXLab();
					double ywin = det->GetLayer(ilr)->GetBgCluster(cix)->GetYLab();
					double zwin = det->GetLayer(ilr)->GetBgCluster(cix)->GetZLab();

					reco_trck_cls_r[irec].push_back(TVector3(xwin, ywin, zwin)); // fill before killing!
					
					// cout << "ilr: " << ilr << " cix: " << cix << " cid: " << cid << endl;
					// det->GetLayer(ilr)->GetBgCluster(cix)->Print();
					
					det->GetLayer(ilr)->GetBgCluster(cix)->Kill();
				}
				
				/// remove the winner cluster from the lookup table
				remove_from_lookup_table(wincls);

				/// save the clusters' id of the winner track
				reco_clusters_id.push_back(win_cls_id);
				
				int maxclssize  = 0;
				int maxclssizex = 0;
				int maxclssizey = 0;
				for(size_t c=0; c<wincls.size(); ++c)
				{
					maxclssize  = (wincls[c].npixl>maxclssize)   ? wincls[c].npixl  : maxclssize;
					maxclssizex = (wincls[c].npixlx>maxclssizex) ? wincls[c].npixlx : maxclssizex;
					maxclssizey = (wincls[c].npixly>maxclssizey) ? wincls[c].npixly : maxclssizey;
				}

				/// reco kinematics etc
				TLorentzVector prec;
				double pxyz[3];
				double xyz[3];
				trw->GetPXYZ(pxyz);
				trw->GetXYZ(xyz);
				TrackPar* trk = trw->GetTrack();
				prec.SetXYZM(pxyz[0], pxyz[1], pxyz[2], meGeV);
				// if(issig4) cout << "Etru=" << Etru << "  -->  Erec=" << pxyz[2] << " and prec.E()=" << prec.E() << endl;
				float chi2dof = trw->GetNormChi2();
				reco_chi2dof.push_back(chi2dof);
				reco_q.push_back(crg);
				reco_p.push_back(prec);
				reco_x.push_back(xyz[0]);
				reco_y.push_back(xyz[1]);
				reco_z.push_back(xyz[2]);
				reco_dErel.push_back((issig4) ? (Etru-prec.E())/Etru : -999);
				reco_dpzrel.push_back((issig4) ? (pztru-prec.Pz())/pztru : -999);
				// reco_trckmar.push_back( TrackMarker3d(trw,0,zLastLayer+1,0.1,trkcol(prec.E())) );
				reco_trckmar.push_back(TrackMarker3d(trw, 0, zLastLayer+1, 1, trkcol(prec.E())));
				reco_trcklin.push_back(TrackLine3d(trw, zLastLayer+1, 1, trkcol(prec.E())));
				
				int nHits = reco_trck_cls_r[irec].size();
				double SnpSig = (trk->GetSigmaSnp2()>0) ? trk->GetSnp()/sqrt(trk->GetSigmaSnp2()) : -1e10;
				double TglSig = (trk->GetSigmaTgl2()>0) ? trk->GetTgl()/sqrt(trk->GetSigmaTgl2()) : -1e10;
				double xVtxSig = (trk->GetSigmaY2()>0)  ? reco_x[irec]/sqrt(trk->GetSigmaY2())    : -1e10;
				double yVtxSig = (trk->GetSigmaZ2()>0)  ? reco_y[irec]/sqrt(trk->GetSigmaZ2())    : -1e10;
				double Px = reco_p[irec].Px();
				double Py = reco_p[irec].Py();
				
				
				/// fill cluster size plots
				for(size_t c=0; c<wincls.size(); ++c)
				{
					unsigned int ilayer_FS = DecodeLayer(wincls[c].clsid);
					unsigned int ilayer_KF = mapFullSim2KFLayer(side, ilayer_FS);
					TString lyrname_KF = islayers[ilayer_KF];
					if(lyrname_KF.Contains("L1I")) { histos["h_rec_csize_L1I_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L1I_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L1I_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L2I")) { histos["h_rec_csize_L2I_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L2I_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L2I_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L3I")) { histos["h_rec_csize_L3I_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L3I_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L3I_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L4I")) { histos["h_rec_csize_L4I_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L4I_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L4I_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L1O")) { histos["h_rec_csize_L1O_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L1O_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L1O_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L2O")) { histos["h_rec_csize_L2O_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L2O_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L2O_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L3O")) { histos["h_rec_csize_L3O_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L3O_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L3O_"+side]->Fill(wincls[c].npixly); }
					if(lyrname_KF.Contains("L4O")) { histos["h_rec_csize_L4O_"+side]->Fill(wincls[c].npixl); histos["h_rec_csizex_L4O_"+side]->Fill(wincls[c].npixlx); histos["h_rec_csizey_L4O_"+side]->Fill(wincls[c].npixly); }
				}
				
				// fill regardless of matching
				histos["h_rec_Nhits_"+side]->Fill(nHits);
				histos["h_rec_chi2dof_"+side]->Fill( reco_chi2dof[irec] );
				histos["h_rec_SnpSig_"+side]->Fill( SnpSig );
				histos["h_rec_TglSig_"+side]->Fill( TglSig );
				histos["h_rec_xVtxSig_"+side]->Fill( xVtxSig );
				histos["h_rec_yVtxSig_"+side]->Fill( yVtxSig );
				histos["h_rec_pz_"+side]->Fill(reco_p[irec].Pz());
				histos["h_rec_pz_"+side+"_log0"]->Fill(reco_p[irec].Pz());
				histos["h_rec_pz_"+side+"_log1"]->Fill(reco_p[irec].Pz());
				histos["h_rec_pz_"+side+"_log2"]->Fill(reco_p[irec].Pz());
				histos["h_rec_pz_"+side+"_log3"]->Fill(reco_p[irec].Pz());
				histos["h_rec_E_"+side]->Fill(reco_p[irec].E());
				histos["h_rec_E_"+side+"_log0"]->Fill(reco_p[irec].E());
				histos["h_rec_E_"+side+"_log1"]->Fill(reco_p[irec].E());
				histos["h_rec_E_"+side+"_log2"]->Fill(reco_p[irec].E());
				histos["h_rec_E_"+side+"_log3"]->Fill(reco_p[irec].E());
				histos["h_rec_px_"+side]->Fill(reco_p[irec].Px());
				histos["h_rec_px_zoom_"+side]->Fill(reco_p[irec].Px());
				histos["h_rec_py_"+side]->Fill(reco_p[irec].Py());
				histos["h_rec_py_zoom_"+side]->Fill(reco_p[irec].Py());
				if(issig4)
				{
					histos["h_resol_rec_dErel_"+side]->Fill(reco_dErel[irec]);
					histos["h_resol_rec_dpzrel_"+side]->Fill(reco_dpzrel[irec]);
				}
				
				
				/// matching
				int nmat = nmatched(itru,wincls);
				bool ismatched = (nmat==wincls.size()); // tight matching (i.e. all track's clusters have the same truid)
				if(ismatched)
				{
					n_match++;
					histos["h_mat_Nhits_"+side]->Fill(nHits);
					histos["h_mat_chi2dof_"+side]->Fill( reco_chi2dof[irec] );
					histos["h_mat_SnpSig_"+side]->Fill( SnpSig );
					histos["h_mat_TglSig_"+side]->Fill( TglSig );
					histos["h_mat_xVtxSig_"+side]->Fill( xVtxSig );
					histos["h_mat_yVtxSig_"+side]->Fill( yVtxSig );
					histos["h_mat_pz_"+side]->Fill(reco_p[irec].Pz());
					histos["h_mat_pz_"+side+"_log0"]->Fill(reco_p[irec].Pz());
					histos["h_mat_pz_"+side+"_log1"]->Fill(reco_p[irec].Pz());
					histos["h_mat_pz_"+side+"_log2"]->Fill(reco_p[irec].Pz());
					histos["h_mat_pz_"+side+"_log3"]->Fill(reco_p[irec].Pz());
					histos["h_mat_E_"+side]->Fill(reco_p[irec].E());
					histos["h_mat_E_"+side+"_log0"]->Fill(reco_p[irec].E());
					histos["h_mat_E_"+side+"_log1"]->Fill(reco_p[irec].E());
					histos["h_mat_E_"+side+"_log2"]->Fill(reco_p[irec].E());
					histos["h_mat_E_"+side+"_log3"]->Fill(reco_p[irec].E());
					histos["h_mat_px_"+side]->Fill(reco_p[irec].Px());
					histos["h_mat_px_zoom_"+side]->Fill(reco_p[irec].Px());
					histos["h_mat_py_"+side]->Fill(reco_p[irec].Py());
					histos["h_mat_py_zoom_"+side]->Fill(reco_p[irec].Py());
					if(issig4)
					{
						histos["h_resol_mat_dErel_"+side]->Fill(reco_dErel[irec]);
						histos["h_resol_mat_dpzrel_"+side]->Fill(reco_dpzrel[irec]);
					}
				}
				else
				{
					histos["h_non_Nhits_"+side]->Fill(nHits);
					histos["h_non_chi2dof_"+side]->Fill( reco_chi2dof[irec] );
					histos["h_non_SnpSig_"+side]->Fill( SnpSig );
					histos["h_non_TglSig_"+side]->Fill( TglSig );
					histos["h_non_xVtxSig_"+side]->Fill( xVtxSig );
					histos["h_non_yVtxSig_"+side]->Fill( yVtxSig );
					histos["h_non_pz_"+side]->Fill(reco_p[irec].Pz());
					histos["h_non_pz_"+side+"_log0"]->Fill(reco_p[irec].Pz());
					histos["h_non_pz_"+side+"_log1"]->Fill(reco_p[irec].Pz());
					histos["h_non_pz_"+side+"_log2"]->Fill(reco_p[irec].Pz());
					histos["h_non_pz_"+side+"_log3"]->Fill(reco_p[irec].Pz());
					histos["h_non_E_"+side]->Fill(reco_p[irec].E());
					histos["h_non_E_"+side+"_log0"]->Fill(reco_p[irec].E());
					histos["h_non_E_"+side+"_log1"]->Fill(reco_p[irec].E());
					histos["h_non_E_"+side+"_log2"]->Fill(reco_p[irec].E());
					histos["h_non_E_"+side+"_log3"]->Fill(reco_p[irec].E());
					histos["h_non_px_"+side]->Fill(reco_p[irec].Px());
					histos["h_non_px_zoom_"+side]->Fill(reco_p[irec].Px());
					histos["h_non_py_"+side]->Fill(reco_p[irec].Py());
					histos["h_non_py_zoom_"+side]->Fill(reco_p[irec].Py());
					if(issig4)
					{
						histos["h_resol_non_dErel_"+side]->Fill(reco_dErel[irec]);
						histos["h_resol_non_dpzrel_"+side]->Fill(reco_dpzrel[irec]);
					}
				}


				/////////////////////////////////////////////////////////
				//////// place holder for kinematic cuts ////////////////
				/////////////////////////////////////////////////////////
				bool pass = true;
				if(pass && maxclssize>10)                          pass = false;
				histos["h_cutflow_"+side]->Fill("Cluster size", (int)pass);
				if(pass && maxclssizex>5)                          pass = false;
				histos["h_cutflow_"+side]->Fill("Cluster size x", (int)pass);
				if(pass && maxclssizey>5)                          pass = false;
				histos["h_cutflow_"+side]->Fill("Cluster size y", (int)pass);
				if(pass && (Px<-0.0015 || Px>+0.0080))             pass = false;
				histos["h_cutflow_"+side]->Fill("p_{x}", (int)pass); /// this is maybe too tuned
				if(pass && (Py<-0.0025 || Py>+0.0025))             pass = false;
				histos["h_cutflow_"+side]->Fill("p_{y}", (int)pass);
				if(pass && reco_chi2dof[irec]>3)                   pass = false;
				histos["h_cutflow_"+side]->Fill("#chi^{2}/N_{DoF}", (int)pass);
				if(pass && SnpSig<-2)                              pass = false;
				histos["h_cutflow_"+side]->Fill("Snp/#sigma(Snp)", (int)pass);
				if(pass && abs(TglSig)>400)                        pass = false;
				histos["h_cutflow_"+side]->Fill("Tgl/#sigma(Tgl)", (int)pass);
				if(pass && (yVtxSig<-0.00025 || yVtxSig>+0.00025)) pass = false;
				histos["h_cutflow_"+side]->Fill("y_{vtx}/#sigma(y_{vtx})", (int)pass);
				if(pass)
				{
					n_selct++;
					// fill if pass regardless of matching
					histos["h_sel_Nhits_"+side]->Fill(nHits);
					histos["h_sel_chi2dof_"+side]->Fill( reco_chi2dof[irec] );
					histos["h_sel_SnpSig_"+side]->Fill( SnpSig );
					histos["h_sel_TglSig_"+side]->Fill( TglSig );
					histos["h_sel_xVtxSig_"+side]->Fill( xVtxSig );
					histos["h_sel_yVtxSig_"+side]->Fill( yVtxSig );
					histos["h_sel_pz_"+side]->Fill(reco_p[irec].Pz());
					histos["h_sel_pz_"+side+"_log0"]->Fill(reco_p[irec].Pz());
					histos["h_sel_pz_"+side+"_log1"]->Fill(reco_p[irec].Pz());
					histos["h_sel_pz_"+side+"_log2"]->Fill(reco_p[irec].Pz());
					histos["h_sel_pz_"+side+"_log3"]->Fill(reco_p[irec].Pz());
					histos["h_sel_E_"+side]->Fill(reco_p[irec].E());
					histos["h_sel_E_"+side+"_log0"]->Fill(reco_p[irec].E());
					histos["h_sel_E_"+side+"_log1"]->Fill(reco_p[irec].E());
					histos["h_sel_E_"+side+"_log2"]->Fill(reco_p[irec].E());
					histos["h_sel_E_"+side+"_log3"]->Fill(reco_p[irec].E());
					histos["h_sel_px_"+side]->Fill(reco_p[irec].Px());
					histos["h_sel_px_zoom_"+side]->Fill(reco_p[irec].Px());
					histos["h_sel_py_"+side]->Fill(reco_p[irec].Py());
					histos["h_sel_py_zoom_"+side]->Fill(reco_p[irec].Py());
					if(issig4)
					{
						histos["h_resol_sel_dErel_"+side]->Fill(reco_dErel[irec]);
						histos["h_resol_sel_dpzrel_"+side]->Fill(reco_dpzrel[irec]);
					}
				}
			} // end of loop on clusters in layer 4
			
			
			/// summarize
			// int mateff = (int)((float)n_trumt / (float)n_truth * 100.);
			int mateff = (n_truth>0) ? (int)((float)n_match / (float)n_truth * 100.) : 0. ;
			cout << "<<<<<<<<<< Performance summary for sample: " << smplname << " >>>>>>>>>>" << endl;
			cout << "Event #" << ibx << ", " << side << ": n_truth=" << n_truth
					<< ", n_cls4=" << (n4I+n4O)
					<< ", n_cls1=" << (n1I+n1O)
						<< ", n_seeds=" << n_seeds
							<< ", n_solve=" << n_solve
								<< ", n_recos=" << n_recos
									<< ", n_selct=" << n_selct
										<< ", n_match=" << n_match
											<< ", counting eff(rec,mat)=" << mateff << "%" << endl;
			cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << endl;

			fOut->cd();
			// tOut->Fill();

			stopwatch.Stop();
			Double_t cputime = stopwatch.CpuTime();
			Double_t realtime = stopwatch.RealTime();
			av_cputime += cputime;
			av_realtime += realtime;
			cout << "Event #" << ibx << ": CPU time=" << cputime << ", Real time=" << realtime << endl;
			if((ibx%outN)==0) printf("Done %d out of %d --> CPUav=%g, REAL=%g\n", ibx, nevents, av_cputime / (ibx+1), av_realtime / (ibx+1));
		} // end of loop on events

		cout << "nEvntsProcessed=" << nEvntsProcessed << endl;
		cout << "Post processing for side: " << side << endl;

		/// file for all histos and tress
		fOut->cd();
		
		// cout << "...normalising cutflow by N_BX for " << side << endl;
		// histos["h_cutflow_"+side]->Scale(1./(float)nevents);

		cout << "...plotting efficiencies for " << side << endl;
		for(size_t h=0; h<htypes.size(); ++h)
		{
			TString htype = htypes[h];
			cout << "......in " << htype << " histos" << endl;
			
			histos["h_ratio_"+htype+"_pz_"+side]->Divide(        histos["h_"+htype+"_pz_"+side],         histos["h_tru_pz_"+side]);
			histos["h_ratio_"+htype+"_pz_"+side+"_log0"]->Divide(histos["h_"+htype+"_pz_"+side+"_log0"], histos["h_tru_pz_"+side+"_log0"]);
			histos["h_ratio_"+htype+"_pz_"+side+"_log1"]->Divide(histos["h_"+htype+"_pz_"+side+"_log1"], histos["h_tru_pz_"+side+"_log1"]);
			histos["h_ratio_"+htype+"_pz_"+side+"_log2"]->Divide(histos["h_"+htype+"_pz_"+side+"_log2"], histos["h_tru_pz_"+side+"_log2"]);
			histos["h_ratio_"+htype+"_pz_"+side+"_log3"]->Divide(histos["h_"+htype+"_pz_"+side+"_log3"], histos["h_tru_pz_"+side+"_log3"]);
			
			histos["h_ratio_"+htype+"_E_"+side]->Divide(        histos["h_"+htype+"_E_"+side],         histos["h_tru_E_"+side]);
			histos["h_ratio_"+htype+"_E_"+side+"_log0"]->Divide(histos["h_"+htype+"_E_"+side+"_log0"], histos["h_tru_E_"+side+"_log0"]);
			histos["h_ratio_"+htype+"_E_"+side+"_log1"]->Divide(histos["h_"+htype+"_E_"+side+"_log1"], histos["h_tru_E_"+side+"_log1"]);
			histos["h_ratio_"+htype+"_E_"+side+"_log2"]->Divide(histos["h_"+htype+"_E_"+side+"_log2"], histos["h_tru_E_"+side+"_log2"]);
			histos["h_ratio_"+htype+"_E_"+side+"_log3"]->Divide(histos["h_"+htype+"_E_"+side+"_log3"], histos["h_tru_E_"+side+"_log3"]);
		}
		
		cout << "...normalising occupancies for " << side << endl;
		for (TMapTSTH2D::iterator it = histos2.begin(); it != histos2.end(); ++it)
		{
			if(it->first.Contains("_occ_"))
			{
				double binareacm2 = (it->second->GetXaxis()->GetBinWidth(1))*(it->second->GetYaxis()->GetBinWidth(1));
				double pixareacm2 = (27*um2cm)*(29*um2cm);
				double npixelsperbin = (binareacm2/pixareacm2);
				it->second->Scale(1./(npixelsperbin*nEvntsProcessed));
			}
		}
		
		/// 1D histos
		cout << "...writing TH1 for " << side << endl;
		for (TMapTSTH1D::iterator it = histos.begin();  it != histos.end();  ++it) it->second->Write();
		
		/// 2D histos
		cout << "...writing TH2 for " << side << endl;
		for (TMapTSTH2D::iterator it = histos2.begin(); it != histos2.end(); ++it) it->second->Write();
		
		cout << "Closing files for " << side << endl;
		/// writeout tree and file
		// tOut->Write();
		fOut->Write();
		fOut->Close();
		cout << "Done! " << side << endl;
	} // end of loop on sides

	return 0;
}
