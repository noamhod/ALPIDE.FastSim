#if !defined(__CINT__) || defined(__MAKECINT__)
#include "KMCDetectorFwd.h"
#include "KMCProbeFwd.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom.h"
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

typedef map<int, int>                TMapii;
typedef map<TString, int >           TMapTSi;
typedef map<TString, float >         TMapTSf;
typedef map<TString, double >        TMapTSd;
typedef map<TString, vector<int> >   TMapTSvi;
typedef map<TString, vector<float> > TMapTSvf;
typedef map<int,TString>             TMapiTS;
typedef map<TString, TH1D* >         TMapTSTH1D;

TString storage =  gSystem->ExpandPathName("$STORAGEDIR");

double vX=0,vY=0,vZ=0; // event vertex
KMCDetectorFwd* det = 0;
vector<double>* zlayer = new vector<double>;
double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;
double meGeV2 = meGeV*meGeV;
double cm2m = 0.01;

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 50; //27.12;   // cm
double Rbeampipe = 2.413; // cm
double RoffsetBfield = 5.7;
double x1L = -RoffsetBfield-Lstave;
double x1R = -RoffsetBfield;       
double x2L = +RoffsetBfield;       
double x2R = +RoffsetBfield+Lstave;
double yUp = +Hstave/2.;
double yDn = -Hstave/2.;
double zDipoleExit = 202.9;
double B  = 1.0; // Tesla if(proc=="trident") else 1.7 
double LB = 1.029;   // meters
double EseedMin = 1.0; // GeV
double EseedMaxBPPP = 16.0; // GeV
double EseedMaxTRIDENT = 5.0; // GeV
//// dipole geometry
double xWdipole = 33.0;
double yHdipole = 10.8;
double z1dipole = 100.0;
double z2dipole = 202.9;
//// uncertainties
float dxAlignmentXFEL = 0.005; //0.005; // cm
float dyAlignmentXFEL = 0.005; //0.005; // cm
float XvariationSign = +1.;
float YvariationSign = +1.;
float dxAlignmentInTray = 0.000; // cm
float dyAlignmentInTray = 0.000; // cm
bool doMisalignmentX = false;
bool doMisalignmentY = false;



// double yDipoleExitMin = -0.05; ## cm --> TODO: need tuning
// double yDipoleExitMax = +0.05; ## cm --> TODO: need tuning
// double xAbsMargins = 0.025; # cm --> TODO: need tuning
// double yAbsMargins = 0.025 if(proc=="bppp") else 0.1 # cm --> TODO: need tuning

vector<TString> sides{"Eside","Pside"};
vector<TString> coord{"x","y","z"};
vector<TString> attri{"id","type"};
TMapiTS         layers = { {300,"L1"}, {310,"L2"}, {320,"L3"}, {330,"L4"} };
TMapTSi         szlayers = { {"L1",300}, {"L2",310}, {"L3",320}, {"L4",330} };
TMapTSi         silayers = { {"L1",1}, {"L2",2}, {"L3",3}, {"L4",4} };
TMapiTS         islayers = { {1,"L1"}, {3,"L2"}, {5,"L3"}, {7,"L4"} };
TMapTSvf        cached_clusters_xyz; /// coordinates (x,y,z)
TMapTSvi        cached_clusters_att; /// attributes (type and id)
TMapii          cached_clusters_all_ids; /// attributes (type and id)

void resetToTridentGeometry()
{
	cout << "Resetting to Trident geometry!" << endl;
	Lstave = 50;   // cm for BPPP or 50 for Trident
	RoffsetBfield = 14; // cm for Trident in in B=1.4T
	x1L = -RoffsetBfield-Lstave;
	x1R = -RoffsetBfield;       
	x2L = +RoffsetBfield;       
	x2R = +RoffsetBfield+Lstave;
}

bool accept(double x, double y)
{
	bool failx = (x<x1L || (x>x1R && x<x2L) || x>x2R);
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
	TLegend* leg = new TLegend(0.12,0.60,0.50,0.80);
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

bool islayer(double z)
{
	for(int j=0 ; j<(int)zlayer->size() ; ++j)
	{
		double dz = abs(zlayer->at(j)-z);
		if(dz<1.e-6) return true;
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
	TPolyMarker3D *polymarker = new TPolyMarker3D(zlayer->size());
	polymarker->SetMarkerColor(col);
	int n = 0;
	for(int i=0;i<nz+1;i++) {
		if(!islayer(zp[i])) continue;
		polymarker->SetPoint(n,xp[i],yp[i],zp[i]);
		n++;
	}
	return polymarker;
}

TPolyLine3D* GetLayer(TString side, double z, Color_t col)
{
	Int_t n=5;
	Double_t xL[] = {x1L,x1L,x1R,x1R,x1L};
	Double_t xR[] = {x2R,x2R,x2L,x2L,x2R};
	Double_t y[] = {yDn,yUp,yUp,yDn,yDn};
	Double_t zC[] = {z,z,z,z,z};
	TPolyLine3D* polyline = 0;
	if(side=="L") polyline = new TPolyLine3D(n,xL,y,zC);
	if(side=="R") polyline = new TPolyLine3D(n,xR,y,zC);
	polyline->SetLineColor(col);
	return polyline;
}

TPolyLine3D* GeDipole(Color_t col)
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

bool accepttrk(vector<TPolyMarker3D*>& polm, int itrk)
{
	/// in acceptance?
	int nlayers = 4;
	int acctrk = 0;
	for (int i=0 ; i<polm[itrk]->GetN() ; i++)
	{
		Double_t xr,yr,zr;
		polm[itrk]->GetPoint(i,xr,yr,zr);
		if(zr<300) continue; //// count only the active layers
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

void prepare_cahced_clusters()
{
	vector<float> vf;
	vector<int> vi;
	for(unsigned int s=0 ; s<sides.size() ; ++s)
	{
		for(TMapiTS::iterator it=layers.begin() ; it!=layers.begin() ; ++it)
		{
			for(unsigned int c=0 ; c<coord.size() ; ++c)
			{
				cached_clusters_xyz.insert( make_pair(coord[c]+"_"+it->second+"_"+sides[s],vf) );
			}
			for(unsigned int a=0 ; a<attri.size() ; ++a)
			{
				cached_clusters_att.insert( make_pair(attri[a]+"_"+it->second+"_"+sides[s],vi) );
			}
		}
	}
}

void clear_cached_clusters()
{
	for(TMapTSvf::iterator it=cached_clusters_xyz.begin() ; it!=cached_clusters_xyz.end() ; ++it) it->second.clear();
	for(TMapTSvi::iterator it=cached_clusters_att.begin() ; it!=cached_clusters_att.end() ; ++it) it->second.clear();
}

int cache_clusters(vector<TPolyMarker3D*>* clusters_xyz, vector<vector<int> >* clusters_type, vector<vector<int> >* clusters_id, TString side, vector<int>* acc=0, int nMaxToCache=1410065407)
{
	int ncached = 0;
	int ntrks = (int)clusters_xyz->size();
	int ntrkmax = (nMaxToCache>0 && nMaxToCache<ntrks) ? nMaxToCache : ntrks;
	for(int i=0 ; i<ntrkmax ; i++)
	{
		if(acc) // if the vector is provided (for background only)
		{
			if(!acc->at(i)) continue; // check acceptnce
		}
		for(Int_t j=0 ; j<clusters_xyz->at(i)->GetN() ; ++j)
		{
			float x,y,z;
			clusters_xyz->at(i)->GetPoint(j,x,y,z); // the clusters
			if(x>0 and side=="Pside") continue;
			if(x<0 and side=="Eside") continue;
			TString sd = (x>0) ? "Eside" : "Pside";
			TString lr = layers[z];
			
			x = (doMisalignmentX) ? x+XvariationSign*dxAlignmentXFEL : x;
			y = (doMisalignmentY) ? y+YvariationSign*dyAlignmentXFEL : y;
			
			cached_clusters_xyz["x_"+lr+"_"+sd].push_back(x);
			cached_clusters_xyz["y_"+lr+"_"+sd].push_back(y);
			cached_clusters_xyz["z_"+lr+"_"+sd].push_back(z);
			cached_clusters_att["type_"+lr+"_"+sd].push_back( clusters_type->at(i)[silayers[lr]-1] );
			cached_clusters_att["id_"+lr+"_"+sd].push_back( clusters_id->at(i)[silayers[lr]-1] );
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

void add_bkg_cluster(int iLayer, float x, float y, float z, int id)
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
		TString slr = it->second;
		int     zlr = it->first;
		for(unsigned int i=0 ; i<cached_clusters_xyz["x_"+slr+"_"+side].size() ; ++i)
		{
			float x = cached_clusters_xyz["x_"+slr+"_"+side][i];
			float y = cached_clusters_xyz["y_"+slr+"_"+side][i];
			float z = cached_clusters_xyz["z_"+slr+"_"+side][i];
			int  id = cached_clusters_att["id_"+slr+"_"+side][i];
			if(zlr==300) add_bkg_cluster(1,x,y,z,id);
			if(zlr==310) add_bkg_cluster(3,x,y,z,id);
			if(zlr==320) add_bkg_cluster(5,x,y,z,id);
			if(zlr==330) add_bkg_cluster(7,x,y,z,id);
		}
	}	
	/// sort clusters
	for(int l=0 ; l<det->GetLayers()->GetEntries() ; l++)
	{
		det->GetLayer(l)->GetMCCluster()->Kill();
		det->GetLayer(l)->SortBGClusters(); /// sort!!!
		
		/// after sorting, need to map the cluster ids to their indices!!!
		bool active = (l==1 or l==3 or l==5 or l==7);
		if(!active) continue;
		for(int n=0 ; n<det->GetLayer(l)->GetNBgClusters() ; ++n)
		{
			KMCClusterFwd *cl = det->GetLayer(l)->GetBgCluster(n);
			int id = cl->GetTrID();
			cached_clusters_all_ids.insert( make_pair(id,n) );
		}
	}
}

void print_all_clusters(TString side, bool doprint = true)
{	
	if(!doprint) return;
	for(int l=1 ; l<=7 ; l+=2) // active layers
	{
		for(int c=0 ; c<det->GetLayer(l)->GetNBgClusters() ; c++)
		{
			KMCClusterFwd* cluster = det->GetLayer(l)->GetBgCluster(c);
			int  id = cluster->GetTrID();
			float x = cluster->GetXLab();
			float y = cluster->GetYLab();
			float z = cluster->GetZLab();
			cout << "side=" << side << ", layer=" << l << ", id=" << id << " --> r={" << x << ", " << y << ", " << z << "}" << endl;
		}
		cout << endl;
	}
	for(TMapii::iterator it=cached_clusters_all_ids.begin() ; it!=cached_clusters_all_ids.end() ; it++)
	{
		cout << "id=" << it->first << " --> index=" << it->second << endl;
	}
}

int fill_output_clusters(TString side, vector<TPolyMarker3D*>& cxyz, vector<int>& ctype, vector<int>& cid)
{
	int nclusters = 0;
	for(int l=1 ; l<=7 ; l+=2) // active layers
	{
		for(int c=0 ; c<det->GetLayer(l)->GetNBgClusters() ; c++)
		{
			KMCClusterFwd* cluster = det->GetLayer(l)->GetBgCluster(c);
			int  id  = cluster->GetTrID();
			float x  = cluster->GetXLab();
			float y  = cluster->GetYLab();
			float z  = cluster->GetZLab();
	
			int idx = getvecindex(id, cached_clusters_att["id_"+islayers[l]+"_"+side]);
			int typ = (idx>=0) ? cached_clusters_att["type_"+islayers[l]+"_"+side][ idx ] : -3;
			// cout << "id=" << id << ", idx=" << idx << ", typ=" << typ << endl;
			if(idx<0) cout << "WARNING: cannot find in=" << id << " in cached clusters vector" << endl;
		
			TPolyMarker3D *point = new TPolyMarker3D();
			point->SetNextPoint(x,y,z);
			
			cxyz.push_back(point);
			ctype.push_back(typ);
			cid.push_back(id);
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

bool check_clusters(unsigned int i1, unsigned int i4, TString side)
{
	float yAbsMargins = 0.02; // cm (a "road" of 200 microns around the line between r4 and r1)
	float xAbsMargins = 0.02; // cm (a "road" of 200 microns around the line between r4 and r1)
	float r1min[3] = { cached_clusters_xyz["x_L1_"+side][i1]-xAbsMargins, cached_clusters_xyz["y_L1_"+side][i1]-yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] };
	float r1max[3] = { cached_clusters_xyz["x_L1_"+side][i1]+xAbsMargins, cached_clusters_xyz["y_L1_"+side][i1]+yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] };
	float r4min[3] = { cached_clusters_xyz["x_L4_"+side][i4]-xAbsMargins, cached_clusters_xyz["y_L4_"+side][i4]-yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] };
	float r4max[3] = { cached_clusters_xyz["x_L4_"+side][i4]+xAbsMargins, cached_clusters_xyz["y_L4_"+side][i4]+yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] };

	/// check possible clusters in layer 2
	float y2min = yofz(r1min,r4min,(float)szlayers["L2"]);
	float y2max = yofz(r1max,r4max,(float)szlayers["L2"]);
	float x2min = xofz(r1min,r4min,(float)szlayers["L2"]);
	float x2max = xofz(r1max,r4max,(float)szlayers["L2"]);
	bool accept2 = false;
	for(unsigned int i2=0 ; i2<cached_clusters_xyz["x_L2_"+side].size() ; ++i2)
	{	
		bool acceptyz = ( cached_clusters_xyz["y_L2_"+side][i2]>=y2min and cached_clusters_xyz["y_L2_"+side][i2]<=y2max );
		if(!acceptyz) continue;
		bool acceptxz = ( cached_clusters_xyz["x_L2_"+side][i2]>=x2min and cached_clusters_xyz["x_L2_"+side][i2]<=x2max );
		if(!acceptxz) continue;
		accept2 = true;
		break;
	}
	if(!accept2) return false;

	float y3min = yofz(r1min,r4min,(float)szlayers["L3"]);
	float y3max = yofz(r1max,r4max,(float)szlayers["L3"]);
	float x3min = xofz(r1min,r4min,(float)szlayers["L3"]);
	float x3max = xofz(r1max,r4max,(float)szlayers["L3"]);
	bool accept3 = false;
	for(unsigned int i3=0 ; i3<cached_clusters_xyz["x_L3_"+side].size() ; ++i3)
	{
		bool acceptyz = ( cached_clusters_xyz["y_L3_"+side][i3]>=y3min and cached_clusters_xyz["y_L3_"+side][i3]<=y3max );
		if(!acceptyz) continue;
		bool acceptxz = ( cached_clusters_xyz["x_L3_"+side][i3]>=x3min and cached_clusters_xyz["x_L3_"+side][i3]<=x3max );
		if(!acceptxz) continue;
		accept3 = true;
		break;
	}
	if(!accept3) return false;

	return true;
}

bool makeseed(TString process, float* r1, float* r4, unsigned int i1, unsigned int i4, TString side, TLorentzVector& p, bool calibrate=false)
{
	if(abs(r1[0])>=abs(r4[0]))            return false; // |x1| must be smaller than |x4|
	if(r1[0]>0 and r4[0]<0)               return false;
	if(r1[0]<0 and r4[0]>0)               return false;
	if(r1[2]==r4[2])                      return false; // if z1=z4...
	float yDipoleExitAbsMax = (process=="bppp") ? 0.2 : 0.7; // cm
	float xDipoleExitAbsMin = (process=="bppp") ? 1.  : 4. ; // cm
	float xDipoleExitAbsMax = (process=="bppp") ? 25. : 30.;	// cm
	float yDipoleExit = yofz(r1,r4,zDipoleExit);
	float xDipoleExit = xofz(r1,r4,zDipoleExit);
	if(abs(yDipoleExit)>yDipoleExitAbsMax) return false; // the track should point to |y|<~0.2 at the dipole exit
	if(abs(xDipoleExit)<xDipoleExitAbsMin) return false; // the track should point to |x|<~1.0 at the dipole exit
	if(abs(xDipoleExit)>xDipoleExitAbsMax) return false;
	if(!check_clusters(i1,i4,side))        return false; // minimum one cluster at layer 2 and one at layer 3

	TRandom rnd;
	rnd.SetSeed();
	double posneg = rnd.Uniform(-1,+1);
	double pxgaus = rnd.Gaus(7.2e-4,5.0e-4);

	double x0 = 0;
	double z0 = zofx(r1,r4,x0);
	double xExit = abs(xofz(r1,r4,zDipoleExit))*cm2m;
	double H = abs((zDipoleExit-z0))*cm2m;
	double R = H*(LB)/xExit + xExit; // look this up in my slides
	double P = 0.3*B*R;
	P = (calibrate) ? P/1.001 : P;
	
	// if(i4==0 and side=="Eside") cout << "z0=" << z0 << ", xExit=" << xExit << ", H=" << H << ", R=" << R << ", P=" << P << endl;

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
	float EseedMax = (process=="bppp") ? EseedMaxBPPP : EseedMaxTRIDENT;	// GeV
	if(p.E()<EseedMin or p.E()>EseedMax) return false;

	return true;
}


int imatched(TPolyMarker3D* mrec, vector<TPolyMarker3D*>* clusters_xyz, TString side, double maxdistance=0.5)
{
	int imindistance = -1;
	double mindistance = 1e10;
	for(unsigned int i=0 ; i<clusters_xyz->size() ; i++)
	{
		double distance = 0;
		for(Int_t jTru=0 ; jTru<clusters_xyz->at(i)->GetN() ; ++jTru)
		{
			double xTru,yTru,zTru;
			clusters_xyz->at(i)->GetPoint(jTru,xTru,yTru,zTru); // the clusters
			if(xTru>0 and side=="Pside") continue;
			if(xTru<0 and side=="Eside") continue;
			for(Int_t jRec=0 ; jRec<mrec->GetN() ; ++jRec)
			{
				double xRec,yRec,zRec;
				mrec->GetPoint(jRec,xRec,yRec,zRec);
				if(abs(zRec-zTru)<1e-2)
				{
					distance += sqrt((xRec-xTru)*(xRec-xTru)+(yRec-yTru)*(yRec-yTru));
					// cout << "i=" << i << ", rTru={" << xTru << "," << yTru << "," << zTru << "}, rRec={" << xRec << "," << yRec << "," << zRec << "} --> distance=" << distance << endl;
				}
			}
		}
		// distance = distance; // sum of sqrt{dx^2+dy^2} from all 4 rec clusters and the truth position of the truth digitised clusters
		if(distance>0 and distance<mindistance and distance<maxdistance)
		{
			imindistance = i;
			mindistance  = distance; 
		}
		// if(distance>0) cout << "i=" << i << ", distance=" << distance << " --> {mindistance=" << mindistance << ", imindistance=" << imindistance << "}"<< endl;
	}
	// cout << "imindistance=" << imindistance << ", mindistance=" << mindistance << endl;
	return imindistance;
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


// void Reconstruction(TString process, int nMaxBkgTrks=-1, int Seed=12345)
// {
// 	cout << "Settings" << endl;
// 	TString setup = "../setup/setupLUXE_"+process+".txt";
// 	gROOT->LoadMacro("Loader.C+");
// 	gRandom->SetSeed(Seed);

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
	if(argc==3 and !((TString)argv[2]).Contains("-evnt=")) { printf("argc=3 but cannot parse %s\n",argv[2]); exit(-1); }
	if(argc==4 and !((TString)argv[3]).Contains("-ntrk=")) { printf("argc=4 but cannot parse %s\n",argv[3]); exit(-1); }
	if(argc==5 and !((TString)argv[5]).Contains("-seed=")) { printf("argc=5 but cannot parse %s\n",argv[4]); exit(-1); }
	//// assign inputs
	TString process = ((TString)argv[1]).ReplaceAll("-proc=",""); // mandatory
	int     evnt    = (argc>2)     ? toint(((TString)argv[2]).ReplaceAll("-evnt=","")) : -1; // job id [optional]
	int     nMaxBkgTrks = (argc>3) ? toint(((TString)argv[2]).ReplaceAll("-ntrk=","")) : -1; // job id [optional]
	int     Seed    = (argc>4)     ? toint(((TString)argv[3]).ReplaceAll("-seed=","")) : 12345; // seed [optional]
	//// print assigned inputs
	cout << "process=" << process << endl;
	cout << "evnt=" << evnt << endl;
	cout << "nMaxBkgTrks=" << nMaxBkgTrks << endl;
	cout << "Seed=" << Seed << endl;
	
	
	
	
	
	// TString proc = process;
	TString eventid = (evnt<0) ? "" : FormatEventID(evnt);
	TStopwatch stopwatch;
	TString setup = "../setup/setupLUXE_"+process+".txt";
	det = new KMCDetectorFwd();
	det->ReadSetup(setup,setup);
	det->ForceLastActiveLayer(det->GetLastActiveLayerITS()); // will not propagate beyond VertexTelescope
	det->SetMinITSHits(det->GetNumberOfActiveLayersITS()); // require hit in every layer
	det->SetMinMSHits(0); // we don't have muon spectrometer
	det->SetMinTRHits(0); // we don't have muon trigger stations
	// max number of seeds on each layer to propagate (per muon track)
	det->SetMaxSeedToPropagate(3000); // relevant only if backgrount is considered
	// set chi2 cuts
	// det->SetMaxChi2Cl(10.);  // max track to cluster chi2
	det->SetMaxChi2Cl(10.);  // max track to cluster chi2
	// det->SetMaxChi2NDF(3.5); // max total chi2/ndf
	det->SetMaxChi2NDF((process=="trident")?15.:5.); // max total chi2/ndf
	det->SetMaxChi2Vtx(20e9);  // fiducial cut on chi2 of convergence to vtx
	// det->SetMaxChi2Vtx(500);  // fiducial cut on chi2 of convergence to vtx
	// IMPORTANT FOR NON-UNIFORM FIELDS
	det->SetDefStepAir(1);
	det->SetMinP2Propagate(0.3); //NA60+
	det->SetIncludeVertex(kTRUE); // count vertex as an extra measured point
	det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
	det->SetApplyBransonPCorrection(-1); // Branson correction, only relevant for setup with MS
	// for reconstruction:
	// det->SetErrorScale(500.);
	det->SetErrorScale( (process=="trident")?500.:200. );
	det->Print();
	// det->BookControlHistos();
	
	zlayer->push_back(0);   //// NOAM --> GET FROM THE SETUP --> IP (vertex)
	zlayer->push_back(100); //// NOAM --> GET FROM THE SETUP --> start of dipol
	zlayer->push_back(200); //// NOAM --> GET FROM THE SETUP --> end of dipol
	zlayer->push_back(300); //// NOAM --> GET FROM THE SETUP --> layer 1
	zlayer->push_back(310); //// NOAM --> GET FROM THE SETUP --> layer 2
	zlayer->push_back(320); //// NOAM --> GET FROM THE SETUP --> layer 3
	zlayer->push_back(330); //// NOAM --> GET FROM THE SETUP --> layer 4

	int outN = (process=="trident") ? 10 : 10;
	B = (process=="trident") ? 1.0 : 1.7;

	if(process=="trident")
	{
		resetToTridentGeometry();
		cout << "Doing only Pside!" << endl;
		sides.clear(); sides.push_back("Pside"); /// do not reconstruct the Eside
	}

	/// get the signal clusters
	cout << "Getting signal clusters from tree" << endl;
	TFile* fSig = new TFile(storage+"/data/root/dig/dig_"+process+"_"+eventid+".root","READ");
	TTree* tSig = (TTree*)fSig->Get("dig");
	int                      sig_ngen          = 0;
	int                      sig_nslv          = 0;
	int                      sig_nacc          = 0;
	vector<double>*          sig_wgt           = 0;
	vector<int>*             sig_crg           = 0;
	vector<float>*           sig_xvtx          = 0;
	vector<float>*           sig_yvtx          = 0;
	vector<float>*           sig_zvtx          = 0;
	vector<TLorentzVector>*  sig_trkp4         = 0;
	vector<int>*             sig_acc           = 0;
	vector<vector<int> >*    sig_clusters_id   = 0;
	vector<vector<int> >*    sig_clusters_type = 0;
	vector<TPolyMarker3D*>*  sig_clusters_xyz  = 0;
	vector<TPolyMarker3D*>*  sig_trkpts        = 0;
	vector<TPolyLine3D*>*    sig_trklin        = 0;
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
	tSig->SetBranchAddress("clusters_id",  &sig_clusters_id);
	tSig->SetBranchAddress("clusters_type",&sig_clusters_type);
	tSig->SetBranchAddress("clusters_xyz", &sig_clusters_xyz);
	tSig->SetBranchAddress("trkpts",       &sig_trkpts);
	tSig->SetBranchAddress("trklin",       &sig_trklin);
	
	/// get the background clusters
	cout << "Getting background clusters from tree" << endl;
	// TChain* tBkg = new TChain("dig");
	// tBkg->Add(storage+"/data/root/dig_"+process+"_bkg_0*.root");
	// cout << "---- TChain content ----" << endl;
	// tBkg->ls();
	// cout << "------------------------" << endl;
	TFile* fBkg = new TFile(storage+"/data/root/dig/dig_"+process+"_bkg_"+eventid+".root","READ");
	TTree* tBkg = (TTree*)fBkg->Get("dig");
	int                      bkg_ngen          = 0;
	int                      bkg_nslv          = 0;
	int                      bkg_nacc          = 0;
	vector<double>*          bkg_wgt           = 0;
	vector<int>*             bkg_crg           = 0;
	vector<float>*           bkg_xvtx          = 0;
	vector<float>*           bkg_yvtx          = 0;
	vector<float>*           bkg_zvtx          = 0;
	vector<TLorentzVector>*  bkg_trkp4         = 0;
	vector<int>*             bkg_acc           = 0;
	vector<vector<int> >*    bkg_clusters_id   = 0;
	vector<vector<int> >*    bkg_clusters_type = 0;
	vector<TPolyMarker3D*>*  bkg_clusters_xyz  = 0;
	vector<TPolyMarker3D*>*  bkg_trkpts        = 0;
	vector<TPolyLine3D*>*    bkg_trklin        = 0;
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
	tBkg->SetBranchAddress("clusters_id",  &bkg_clusters_id);
	tBkg->SetBranchAddress("clusters_type",&bkg_clusters_type);
	tBkg->SetBranchAddress("clusters_xyz", &bkg_clusters_xyz);
	tBkg->SetBranchAddress("trkpts",       &bkg_trkpts);
	tBkg->SetBranchAddress("trklin",       &bkg_trklin);
	

	// output tree
	cout << "Setting the output tree" << endl;
	gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",   "vector");
	gSystem->Exec("mkdir -p "+storage+"/data/root/rec");
	TFile* fOut = new TFile(storage+"/data/root/rec/rec_"+process+"_"+eventid+".root","RECREATE");
	TTree* tOut = new TTree("reco","reco");
	/// all clusters output branches
	vector<TPolyMarker3D*>  all_clusters_xyz;
	vector<int>             all_clusters_type;
	vector<int>             all_clusters_id;
	tOut->Branch("all_clusters_xyz",  &all_clusters_xyz);
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
	vector<float>            reco_q;
	vector<TLorentzVector>   reco_p;
	vector<float>            reco_x;
	vector<float>            reco_y;
	vector<float>            reco_z;
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
 
	/// prepare the dictionaries
	prepare_cahced_clusters();
 
	/// for timing
	Double_t av_cputime  = 0;
	Double_t av_realtime = 0;
 
	/// loop on events
	Int_t nsigevents = tSig->GetEntries();
	Int_t nbkgevents = tBkg->GetEntries();
	if(nbkgevents<nsigevents)
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
		tBkg->GetEntry(iev); /// background
		
		////////////////////////////////////////////
		/// clear output vectors: digitized clusters
		for(unsigned int x=0 ; x<all_clusters_xyz.size() ; ++x) delete all_clusters_xyz[x];
		all_clusters_xyz.clear();
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
		for(unsigned int x=0 ; x<reco_trckmar.size() ; ++x) delete reco_trckmar[x];
		for(unsigned int x=0 ; x<reco_trcklin.size() ; ++x) delete reco_trcklin[x];
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
			vector<int> vtruid{ sig_clusters_id->at(t)[0],sig_clusters_id->at(t)[1],sig_clusters_id->at(t)[2],sig_clusters_id->at(t)[3] };
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
		int nbtrks = (int)bkg_crg->size();
		int nbmax = (nMaxBkgTrks>0 && nMaxBkgTrks<nbtrks) ? nMaxBkgTrks : nbtrks;
		for(int b=0 ; b<nbmax ; ++b)
		{
			if(!bkg_acc->at(b)) continue; // ignore tracks out of acceptance!
			
			vector<int> vbkgid{ bkg_clusters_id->at(b)[0],bkg_clusters_id->at(b)[1],bkg_clusters_id->at(b)[2],bkg_clusters_id->at(b)[3] };
			bkgr_clusters_id.push_back( vbkgid );
			bkgr_acc.push_back( bkg_acc->at(b) );
			bkgr_wgt.push_back( bkg_wgt->at(b) );
			bkgr_x.push_back( bkg_xvtx->at(b) );
			bkgr_y.push_back( bkg_yvtx->at(b) );
			bkgr_z.push_back( bkg_zvtx->at(b) );
			bkgr_q.push_back( bkg_crg->at(b) );
			bkgr_p.push_back( bkg_trkp4->at(b) );
		   bkgr_trckmar.push_back( bkg_trkpts->at(b) );			
			// for(unsigned int xx=0 ; xx<bkg_trkpts->at(b)->GetN() ; xx++)
			// {
			// 	Double_t x,y,z;
			// 	bkg_trkpts->at(b)->GetPoint(xx,x,y,z);
			// 	cout << "background track: #" << b << " point #" << xx << ", xyz={"<<x<<","<<y<<","<<z<<"}" << endl;
			// }
		   bkgr_trcklin.push_back( bkg_trklin->at(b) );
		}
		
		/////////////////////
		/// loop on the sides
		for(unsigned int s=0 ; s<sides.size() ; ++s)
		{
			TString side = sides[s];
			
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
			int ncached_signal_clusters = cache_clusters(sig_clusters_xyz,sig_clusters_type,sig_clusters_id,side);

			/// make a pool of all background and noise clusters
			int ncached_background_clusters = cache_clusters(bkg_clusters_xyz,bkg_clusters_type,bkg_clusters_id,side,bkg_acc,nMaxBkgTrks);
			
			/// rest all the layers of the detector (including inactive if any)
			reset_layers_all(); // reset both sides 
			
			/// add all clusters to the detector
			add_all_clusters(side);
			print_all_clusters(side,false);
			
			/// write out all clusters when these are sorted
			int all_clusters = fill_output_clusters(side,all_clusters_xyz,all_clusters_type,all_clusters_id);
			// cout << "ncached_signal_clusters=" << ncached_signal_clusters << ", ncached_background_clusters=" << ncached_background_clusters << ", all_clusters=" << all_clusters << endl;
			
			/// offset for signal id's !!!
			int sigoffset = 100000; // should be multiplied by the layer number
			
			/// run over all clusters of layer 4 in the pool --> these are the seeds for the KalmanFilter fit
			for(unsigned int i4=0 ; i4<cached_clusters_xyz["x_L4_"+side].size() ; ++i4)
			// for(unsigned int i4=0 ; i4<1 ; ++i4)
			{
				// reset all tracks from all layers
				reset_layers_tracks();
				
				vector<TLorentzVector> pseeds;
				for(unsigned int i1=0 ; i1<cached_clusters_xyz["x_L1_"+side].size() ; ++i1)
				{	
					// reset all tracks from all layers but layer 0
					reset_layers_tracks(0);
					
					/// find the momentum of the seed
					TLorentzVector pseed;
					bool docalibrate = true;
					float r1[3] = {cached_clusters_xyz["x_L1_"+side][i1], cached_clusters_xyz["y_L1_"+side][i1], cached_clusters_xyz["z_L1_"+side][i1]};
					float r4[3] = {cached_clusters_xyz["x_L4_"+side][i4], cached_clusters_xyz["y_L4_"+side][i4], cached_clusters_xyz["z_L4_"+side][i4]};
					bool seed = makeseed(process,r1,r4,i1,i4,side,pseed,docalibrate);
					if(!seed) continue; // cannot make a meaningful seed
					pseeds.push_back(pseed);
					bool issig  = ( cached_clusters_att["type_L1_"+side][i1]==1 and cached_clusters_att["type_L4_"+side][i4]==1 );
					bool sameid = ( (cached_clusters_att["id_L1_"+side][i1]-1*sigoffset)==(cached_clusters_att["id_L4_"+side][i4]-4*sigoffset) );
					seed_type.push_back( issig and sameid );
					vector<int> vidseed{cached_clusters_att["id_L1_"+side][i1],-1,-1,cached_clusters_att["id_L4_"+side][i4]};
					seed_clusters_id.push_back(vidseed);
					n_seeds++;
		   		
					seed_q.push_back(crg);
					seed_p.push_back(pseed);
				} // end of loop on clusters in layer 1
				if(n_seeds<1) continue;
				
				// bool doPrint0 = (cached_clusters_att["id_L4_"+side][i4]==0+4*sigoffset and cached_clusters_att["type_L4_"+side][i4]==1 and side=="Eside");
				// bool doPrint2 = (cached_clusters_att["id_L4_"+side][i4]==2+4*sigoffset and cached_clusters_att["type_L4_"+side][i4]==1 and side=="Eside");
				// bool doPrint = (doPrint0 or doPrint2);
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
				
				/// get the clusters of the winner tracK
				int win_cls_id1 = trw->GetClID(1); // provide active layer ID, not the physical ones (most are passive)
				int win_cls_id2 = trw->GetClID(2); // provide active layer ID, not the physical ones (most are passive)
				int win_cls_id3 = trw->GetClID(3); // provide active layer ID, not the physical ones (most are passive)
				int win_cls_id4 = trw->GetClID(4); // provide active layer ID, not the physical ones (most are passive)
				int win_cls_inx1 = cached_clusters_all_ids[win_cls_id1];
				int win_cls_inx2 = cached_clusters_all_ids[win_cls_id2];
				int win_cls_inx3 = cached_clusters_all_ids[win_cls_id3];
				int win_cls_inx4 = cached_clusters_all_ids[win_cls_id4];
				if(doPrint) cout << "going to kill: id1="<<win_cls_id1<<", id2="<<win_cls_id2<<", id3="<<win_cls_id3<<", id4="<<win_cls_id4<<endl;
				if(doPrint) cout << "               ix1="<<win_cls_inx1<<", ix2="<<win_cls_inx2<<", ix3="<<win_cls_inx3<<", ix4="<<win_cls_inx4<<endl;
				if(win_cls_id1>0) { det->GetLayer(1)->GetBgCluster( win_cls_inx1 )->Kill(); }
				if(win_cls_id2>0) { det->GetLayer(3)->GetBgCluster( win_cls_inx2 )->Kill(); }
				if(win_cls_id3>0) { det->GetLayer(5)->GetBgCluster( win_cls_inx3 )->Kill(); }
				if(win_cls_id4>0) { det->GetLayer(7)->GetBgCluster( win_cls_inx4 )->Kill(); }
				if(doPrint) det->CheckClusters(win_cls_inx1,win_cls_inx2,win_cls_inx3,win_cls_inx4);
				
				vector<int> vrecid{ win_cls_id1,win_cls_id2,win_cls_id3,win_cls_id4 };
				reco_clusters_id.push_back( vrecid );
				
				TLorentzVector prec;
				double pxyz[3];
				double xyz[3];
				trw->GetPXYZ(pxyz);
				trw->GetXYZ(xyz);
				prec.SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
				int ismatched =  0; // TODO: GET THE WINNER CLUSTERS AND CHECK ALL LAYERS
				int ixmatched = -1; // TODO: GET THE WINNER CLUSTERS AND CHECK ALL LAYERS
				int idmatched = -1; // TODO: GET THE WINNER CLUSTERS AND CHECK ALL LAYERS
				float chi2dof = trw->GetNormChi2();
				reco_chi2dof.push_back( chi2dof );
				reco_q.push_back( crg );
				reco_p.push_back( prec );
				reco_x.push_back( xyz[0] );
				reco_y.push_back( xyz[1] );
				reco_z.push_back( xyz[2] );
				reco_trckmar.push_back( TrackMarker3d(trw,0,361,0.1,trkcol(prec.E())) );
				reco_trcklin.push_back( TrackLine3d(trw,361,1,trkcol(prec.E())) );

				/// TODO: this is a test
				unsigned int irec = reco_trckmar.size()-1;
				Int_t imatch = imatched(reco_trckmar[irec],sig_clusters_xyz,side);
				if(imatch>=0)
				{
					ismatched = 1;
					ixmatched = imatch;
					idmatched = true_clusters_id[imatch][0];
					true_rec_imatch[imatch].push_back( irec );
					n_match++;
					// cout << "Ntru=" << n_truth << ", Nclsperlyr=" << ncached_signal_clusters/4 << ", Etru=" << sig_trkp4->at(imatch).E() << " GeV, Erec=" << prec.E() << "GeV --> imatch=" << imatch << ": win_cls_id1=" << win_cls_id1 << ", win_cls_id2=" << win_cls_id2 << ", win_cls_id3=" << win_cls_id3 << ", win_cls_id4=" << win_cls_id4 << endl;
				}
				else
				{
					ismatched =  0;
					ixmatched = -1;
					idmatched = -1;
					// cout << "Ntru=" << n_truth << ", Nclsperlyr=" << ncached_signal_clusters/4 << ", Etru=!!!NOT MATCHED!!!, Erec=" << prec.E() << "GeV --> win_cls_id1=" << win_cls_id1 << ", win_cls_id2=" << win_cls_id2 << ", win_cls_id3=" << win_cls_id3 << ", win_cls_id4=" << win_cls_id4 << endl;
				}
				reco_ismtchd.push_back( ismatched );
				reco_ixmtchd.push_back( ixmatched );
				reco_idmtchd.push_back( idmatched );
				
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
			
			////////////////////////////////////////////////////////////
			/// post-processing per side histos to fill
			for(unsigned int t=0 ; t<true_q.size() ; ++t)
			{
				if(side=="Eside" and true_q[t]>0) continue;
				if(side=="Pside" and true_q[t]<0) continue;
				if(true_rec_imatch[t].size()>0) n_trumt++;
				// else cout << "This track is not matched: Etru[" << t << "]=" << true_p[t].E() << " GeV" << endl;
				histos["h_E_tru_all_"+side]->Fill( true_p[t].E() );
				int truid1 = true_clusters_id[t][0];
				int truid4 = true_clusters_id[t][3];
				for(unsigned int s=0 ; s<seed_p.size() ; ++s)
				{
					if(seed_type[s]!=1)                continue; // has to be signal track
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
		} // end of loop on sides
		
		fOut->cd();
		tOut->Fill();
		
		stopwatch.Stop();
		Double_t cputime  = stopwatch.CpuTime();
		Double_t realtime = stopwatch.RealTime();
		av_cputime  += cputime;
		av_realtime += realtime;
		cout << "Event #" << iev << ": CPU time=" << cputime << ", Real time=" << realtime << endl;
		if((iev%outN)==0) printf("Done %d out of %d --> CPUav=%g, REAL=%g\n",iev,nsigevents,av_cputime/(iev+1),av_realtime/(iev+1));
	}
	
	histos["h_E_eff_sed_Eside"]->Divide(histos["h_E_tru_sed_mat_Eside"],histos["h_E_tru_all_Eside"]);
	histos["h_E_eff_sed_Pside"]->Divide(histos["h_E_tru_sed_mat_Pside"],histos["h_E_tru_all_Pside"]);
	histos["h_E_eff_rec_Eside"]->Divide(histos["h_E_tru_rec_mat_Eside"],histos["h_E_tru_all_Eside"]);
	histos["h_E_eff_rec_Pside"]->Divide(histos["h_E_tru_rec_mat_Pside"],histos["h_E_tru_all_Pside"]);
	fOut->cd();
	tOut->Write();
	for(TMapTSTH1D::iterator it=histos.begin() ; it!=histos.end() ; ++it) it->second->Write();
	fOut->Write();
	fOut->Close();
	
	return 0;
}
