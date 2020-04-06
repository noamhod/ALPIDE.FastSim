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
#include "TTreeStream.h"
#include "TList.h"
#include "TVector2.h"
#include "TRandom.h"
#include "TObject.h"
#include <iterator> 
#include <map>
#endif

typedef map<TString, int >           TMapTSi;
typedef map<TString, float >         TMapTSf;
typedef map<TString, double >        TMapTSd;
typedef map<TString, vector<int> >   TMapTSvi;
typedef map<TString, vector<float> > TMapTSvf;
typedef map<int,TString>             TMapiTS;

double vX=0,vY=0,vZ=0; // event vertex
KMCDetectorFwd* det = 0;
vector<double>* zlayer = new vector<double>;
double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;
double meGeV2 = meGeV*meGeV;
double cm2m = 0.01;

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 27;   // cm for BPPP or 50 for Trident
double Rbeampipe = 4; // cm for BPPP or 14 for Trident
double RoffsetBfield22BPPP = 7.0; // cm for BPPP in B=2.2T
double RoffsetBfield20BPPP = 5.7; // cm for BPPP in B=2.0T
double RoffsetBfield14BPPP = 4.0; // cm for BPPP in B=1.4T
double RoffsetBfield = RoffsetBfield20BPPP;
double x1L = -RoffsetBfield-Lstave;
double x1R = -RoffsetBfield;       
double x2L = +RoffsetBfield;       
double x2R = +RoffsetBfield+Lstave;
double yUp = +Hstave/2.;
double yDn = -Hstave/2.;
double zDipoleExit = 202.9;
double B  = 1.4; // Tesla if(proc=="trident") else 2.0 
double LB = 1;   // meters
double EseedMin = 1.0; // GeV
double EseedMax = 18.0; // GeV

// double yDipoleExitMin = -0.05; ## cm --> TODO: need tuning
// double yDipoleExitMax = +0.05; ## cm --> TODO: need tuning
// double xAbsMargins = 0.025; # cm --> TODO: need tuning
// double yAbsMargins = 0.025 if(proc=="bppp") else 0.1 # cm --> TODO: need tuning

//// dipole geometry
double xW = 120;
double yH = 67.2;
double z1 = 102.9;
double z2 = 202.9;

vector<TString> sides{"Eside","Pside"};
vector<TString> coord{"x","y","z"};
vector<TString> attri{"id","type"};
TMapiTS         layers = {{300,"L1"}, {310,"L2"}, {320,"L3"}, {330,"L4"} };
TMapTSvf        cached_clusters_xyz; /// coordinates (x,y,z)
TMapTSvi        cached_clusters_att; /// attributes (type and id)

void resetToTridentGeometry()
{
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
	polyline->SetPoint(0,-xW/2,-yH/2,z1);
	polyline->SetPoint(1,-xW/2,+yH/2,z1);
	polyline->SetPoint(2,+xW/2,+yH/2,z1);
	polyline->SetPoint(3,+xW/2,-yH/2,z1);
	polyline->SetPoint(4,-xW/2,-yH/2,z1);

	polyline->SetPoint(5,-xW/2,-yH/2,z2); // go up

	polyline->SetPoint(6,-xW/2,+yH/2,z2); // move
	polyline->SetPoint(7,-xW/2,+yH/2,z1); // go down
	polyline->SetPoint(8,-xW/2,+yH/2,z2); // up again

	polyline->SetPoint(9,+xW/2,+yH/2,z2); // move
	polyline->SetPoint(10,+xW/2,+yH/2,z1); // go down
	polyline->SetPoint(11,+xW/2,+yH/2,z2); // up again

	polyline->SetPoint(12,+xW/2,-yH/2,z2); // move
	polyline->SetPoint(13,+xW/2,-yH/2,z1); // go down
	polyline->SetPoint(14,+xW/2,-yH/2,z2); // up again

	polyline->SetPoint(15,-xW/2,-yH/2,z2); // move
	polyline->SetPoint(16,-xW/2,-yH/2,z1); // go down
	polyline->SetPoint(17,-xW/2,-yH/2,z2); // up again

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

void cache_signal_clusters(vector<TPolyMarker3D*>* polm_clusters, TString side)
{
	for(unsigned int i=0 ; i<polm_clusters->size() ; i++)
	{
		for(Int_t j=0 ; j<polm_clusters->at(i)->GetN() ; ++j)
		{
			float x,y,z;
			polm_clusters->at(i)->GetPoint(j,x,y,z); // the clusters
			if(x>0 and side=="Pside") continue;
			if(x<0 and side=="Eside") continue;
			TString sd = (x>0) ? "Eside" : "Pside";
			TString lr = layers[z];
			cached_clusters_xyz["x_"+lr+"_"+sd].push_back(x);
			cached_clusters_xyz["y_"+lr+"_"+sd].push_back(y);
			cached_clusters_xyz["z_"+lr+"_"+sd].push_back(z);
			cached_clusters_att["type_"+lr+"_"+sd].push_back(1);
			cached_clusters_att["id_"+lr+"_"+sd].push_back(i);
		}
	}
}

void cache_background_clusters(vector<TPolyMarker3D*>* clusters_xyz, vector<int>* clusters_type, vector<int>* clusters_id, TString side)
{	
	for(unsigned int i=0 ; i<clusters_xyz->size() ; i++)
	{
		for(Int_t j=0 ; j<clusters_xyz->at(i)->GetN() ; ++j)
		{
			float x,y,z;
			clusters_xyz->at(i)->GetPoint(j,x,y,z); // the clusters
			if(x>0 and side=="Pside") continue;
			if(x<0 and side=="Eside") continue;
			TString sd = (x>0) ? "Eside" : "Pside";
			TString lr = layers[z];
			cached_clusters_xyz["x_"+lr+"_"+sd].push_back(x);
			cached_clusters_xyz["y_"+lr+"_"+sd].push_back(y);
			cached_clusters_xyz["z_"+lr+"_"+sd].push_back(z);
			cached_clusters_att["type_"+lr+"_"+sd].push_back( clusters_type->at(i) );
			cached_clusters_att["id_"+lr+"_"+sd].push_back( clusters_id->at(i) );
		}
	}
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
			if(zlr==300) add_bkg_cluster(1,x,y,z,i);
			if(zlr==310) add_bkg_cluster(3,x,y,z,i);
			if(zlr==320) add_bkg_cluster(5,x,y,z,i);
			if(zlr==330) add_bkg_cluster(7,x,y,z,i);
		}
	}
	/// sort clusters
	for(int l=0 ; l<det->GetLayers()->GetEntries() ; l++)
	{
		det->GetLayer(l)->GetMCCluster()->Kill();
		det->GetLayer(l)->SortBGClusters();
	}
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

vector<float> getxminmax(float* r1min, float* r1max, float* r4min, float* r4max, float zTest)
{
	vector<float> xminmax{ xofz(r1min,r4min,zTest), xofz(r1max,r4max,zTest) };
	return xminmax;
}

vector<float> getyminmax(float* r1min, float* r1max, float* r4min, float* r4max, float zTest)
{
	vector<float> yminmax{ yofz(r1min,r4min,zTest), yofz(r1max,r4max,zTest) };
	return yminmax;
}

bool check_clusters(unsigned int i1, unsigned int i4, TString side)
{
	float yAbsMargins = 0.005; /// cm
	float r1min[3] = { cached_clusters_xyz["x_L1_"+side][i1], cached_clusters_xyz["y_L1_"+side][i1]-yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] };
	float r1max[3] = { cached_clusters_xyz["x_L1_"+side][i1], cached_clusters_xyz["y_L1_"+side][i1]+yAbsMargins, cached_clusters_xyz["z_L1_"+side][i1] };
	float r4min[3] = { cached_clusters_xyz["x_L4_"+side][i4], cached_clusters_xyz["y_L4_"+side][i4]-yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] };
	float r4max[3] = { cached_clusters_xyz["x_L4_"+side][i4], cached_clusters_xyz["y_L4_"+side][i4]+yAbsMargins, cached_clusters_xyz["z_L4_"+side][i4] };

	/// check possible clusters in layer 2
	bool accept2 = false;
	for(unsigned int i2=0 ; i2<cached_clusters_xyz["x_L2_"+side].size() ; ++i2)
	{	
		float y2min = yofz(r1min,r4min, cached_clusters_xyz["z_L2_"+side][i2]);
		float y2max = yofz(r1max,r4max, cached_clusters_xyz["z_L2_"+side][i2]);
		bool acceptyz = ( cached_clusters_xyz["y_L2_"+side][i2]>=y2min and cached_clusters_xyz["y_L2_"+side][i2]<=y2max );
		if(!acceptyz) continue;

		accept2 = true;
		break;
	}
	if(!accept2) return false;

	bool accept3 = false;
	for(unsigned int i3=0 ; i3<cached_clusters_xyz["x_L3_"+side].size() ; ++i3)
	{
		float y3min = yofz(r1min,r4min, cached_clusters_xyz["z_L3_"+side][i3]);
		float y3max = yofz(r1max,r4max, cached_clusters_xyz["z_L3_"+side][i3]);
		bool acceptyz = ( cached_clusters_xyz["y_L3_"+side][i3]>=y3min and cached_clusters_xyz["y_L3_"+side][i3]<=y3max );
		if(!acceptyz) continue;

		accept3 = true;
		break;
	}
	if(!accept3) return false;

	return true;
}

bool makeseed(TString process, float* r1, float* r4, unsigned int i1, unsigned int i4, TString side, TLorentzVector& p)
{
	if(abs(r1[0])>=abs(r4[0]))            return false; // |x1| must be smaller than |x4|
	if(r1[0]>0 and r4[0]<0)               return false;
	if(r1[0]<0 and r4[0]>0)               return false;
	if(r1[2]==r4[2])                      return false; // if z1=z4...
	float yDipoleExitAbsMax = 0.1; // cm
	float xDipoleExitAbsMin = 0.1; // cm
	float xDipoleExitAbsMax = (process=="bppp") ? 30. : 35.;	// cm
	float yDipoleExit = yofz(r1,r4,zDipoleExit);
	float xDipoleExit = xofz(r1,r4,zDipoleExit);
	if(abs(yDipoleExit)>yDipoleExitAbsMax) return false; // the track should point to y~0 at the dipole exit
	if(abs(xDipoleExit)<xDipoleExitAbsMin) return false;
	if(abs(xDipoleExit)>xDipoleExitAbsMax) return false;
	// if(!check_clusters(i1,i4,side))        return false; // minimum one cluster at layer 2 and one at layer 3

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
	// cout << "z0=" << z0 << ", xExit=" << xExit << ", H=" << H << ", R=" << R << ", P=" << P << endl;

	TVector2 v1(r1[2],r1[1]);
	TVector2 v4(r4[2],r4[1]);
	TVector2 u = rUnit2(v1,v4);
	double uz = u.X();
	double uy = u.Y();
	double px = (posneg>=0) ? pxgaus : -pxgaus;
	double py = P*uy;
	double pz = P*uz;
	p.SetPxPyPzE(px,py,pz,TMath::Sqrt(px*px + py*py + pz*pz + meGeV2));
	// cout << "px=" << px << ", py=" << py << ", pz=" << pz << endl;
	
	if(p.E()<EseedMin or p.E()>EseedMax) return false;

	return true;
}



void runLUXEeeRecoFromClusters(TString process, int Seed=12345) //, const char* setup="setup/setupLUXE.txt")
{
	cout << "Settings" << endl;
	TString setup = "setup/setupLUXE_"+process+".txt";
	gROOT->LoadMacro("Loader.C+");
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
	// det->SetMaxChi2Vtx(20e9);  // fiducial cut on chi2 of convergence to vtx
	det->SetMaxChi2Vtx(200);  // fiducial cut on chi2 of convergence to vtx
	// IMPORTANT FOR NON-UNIFORM FIELDS
	det->SetDefStepAir(1);
	det->SetMinP2Propagate(0.3); //NA60+
	det->SetIncludeVertex(kTRUE); // count vertex as an extra measured point
	// det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
	det->SetApplyBransonPCorrection(-1); // Branson correction, only relevant for setup with MS
	det->Print();
	// det->BookControlHistos();
   
	zlayer->push_back(0);   //// NOAM --> GET FROM THE SETUP   --> IP (vertex)
	zlayer->push_back(100); //// NOAM --> GET FROM THE SETUP --> start of dipol
	zlayer->push_back(200); //// NOAM --> GET FROM THE SETUP --> end of dipol
	zlayer->push_back(300); //// NOAM --> GET FROM THE SETUP
	zlayer->push_back(310); //// NOAM --> GET FROM THE SETUP
	zlayer->push_back(320); //// NOAM --> GET FROM THE SETUP
	zlayer->push_back(330); //// NOAM --> GET FROM THE SETUP

	int outN = (process=="trident") ? 10 : 10;
	B = (process=="trident") ? B : 2.0;

	// TString process = "bppp";  /// trident or bppp or bppp_bkg or trident_bkg
	if(process=="trident") resetToTridentGeometry();

	/// get the signal clusters
	cout << "Getting signal clusters from tree" << endl;
	TFile* fSig = new TFile("data/root/rec_"+process+".root","READ");
	TTree* tSig = (TTree*)fSig->Get("res");
	vector<double>*          wgtgen = 0;
	vector<TLorentzVector>*  pgen = 0;
	vector<int>*             qgen = 0;
	vector<TPolyMarker3D*>*  polm_gen = 0;
	vector<TPolyLine3D*>*    poll_gen = 0;
	vector<int>*             acctrkgen = 0;
	vector<int>*             igentrk = 0;
	vector<TPolyMarker3D*>*  polm_clusters = 0;
	tSig->SetBranchAddress("wgtgen",                 &wgtgen);
	tSig->SetBranchAddress("qgen",                   &qgen);
	tSig->SetBranchAddress("pgen",                   &pgen);
	tSig->SetBranchAddress("igentrk",                &igentrk);
	tSig->SetBranchAddress("acctrkgen",              &acctrkgen);
	tSig->SetBranchAddress("polm_clusters",          &polm_clusters);
	tSig->SetBranchAddress("polm_gen",               &polm_gen);
	tSig->SetBranchAddress("poll_gen",               &poll_gen);
	
	/// get the background and noise clusters
	cout << "Getting background clusters from tree" << endl;
	TFile* fBkg = new TFile("data/root/background_clusters_"+process+".root","READ");
	TTree* tBkg = (TTree*)fBkg->Get("clusters");
	vector<TPolyMarker3D*> *clusters_xyz  = 0;
	vector<int>            *clusters_type = 0;
	vector<int>            *clusters_id   = 0;
	tBkg->SetBranchAddress("clusters_xyz",  &clusters_xyz);
	tBkg->SetBranchAddress("clusters_type", &clusters_type);
	tBkg->SetBranchAddress("clusters_id",   &clusters_id);

	// output tree
	cout << "Setting the output tree" << endl;
	gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
	TFile* fOut = new TFile("data/root/rec_from_clusters_"+process+".root","RECREATE");
	TTree* tOut = new TTree("reco","reco");
	/// all clusters output branches
	// vector<TPolyMarker3D*>  all_clusters_xyz;
	// vector<int>             all_clusters_type;
	// vector<int>             all_clusters_id;
	// tOut->Branch("all_clusters_xyz",  &all_clusters_xyz);
	// tOut->Branch("all_clusters_type", &all_clusters_type);
	// tOut->Branch("all_clusters_id",   &all_clusters_id);
	/// truth output branches
	vector<float>            true_q;
	vector<TLorentzVector>   true_p;
	vector<TPolyMarker3D*>   true_trckmar;
	vector<TPolyLine3D*>     true_trcklin;
	// vector<int>              true_cluster4_id;
	// vector<int>              true_cluster3_id;
	// vector<int>              true_cluster2_id;
	// vector<int>              true_cluster1_id;
	tOut->Branch("true_q",       &true_q);
	tOut->Branch("true_p",       &true_p);
	tOut->Branch("true_trckmar", &true_trckmar);
	tOut->Branch("true_trcklin", &true_trcklin);
	// tOut->Branch("true_cluster4_id", &true_cluster4_id);
	// tOut->Branch("true_cluster3_id", &true_cluster3_id);
	// tOut->Branch("true_cluster2_id", &true_cluster2_id);
	// tOut->Branch("true_cluster1_id", &true_cluster1_id);
	/// background tracks output branches
	vector<TPolyMarker3D*>   bkgr_trckmar;
	vector<TPolyLine3D*>     bkgr_trcklin;
	// vector<int>              bkgr_cluster4_id;
	// vector<int>              bkgr_cluster3_id;
	// vector<int>              bkgr_cluster2_id;
	// vector<int>              bkgr_cluster1_id;
	tOut->Branch("bkgr_trckmar", &bkgr_trckmar);
	tOut->Branch("bkgr_trcklin", &bkgr_trcklin);
	// tOut->Branch("bkgr_cluster4_id", &bkgr_cluster4_id);
	// tOut->Branch("bkgr_cluster3_id", &bkgr_cluster3_id);
	// tOut->Branch("bkgr_cluster2_id", &bkgr_cluster2_id);
	// tOut->Branch("bkgr_cluster1_id", &bkgr_cluster1_id);
	/// seeds output branches
	vector<float>          seed_q;
	vector<TLorentzVector> seed_p;
	// vector<int>            seed_cluster4_id;
	// vector<int>            seed_cluster1_id;
	tOut->Branch("seed_q", &seed_q);
	tOut->Branch("seed_p", &seed_p);
	// tOut->Branch("seed_cluster4_id", &seed_cluster4_id);
	// tOut->Branch("seed_cluster1_id", &seed_cluster1_id);
	/// reconstructed clusters output branches
	vector<float>            reco_q;
	vector<TLorentzVector>   reco_p;
	vector<TPolyMarker3D*>   reco_trckmar;
	vector<TPolyLine3D*>     reco_trcklin;
	vector<float>            reco_chi2dof;
	vector<int>              reco_ismtchd;
	vector<int>              reco_idmtchd;
	// vector<int>              reco_cluster4_id;
	// vector<int>              reco_cluster1_id;
	// vector<int>              reco_cluster2_id;
	// vector<int>              reco_cluster3_id;
	tOut->Branch("reco_q",       &reco_q);
	tOut->Branch("reco_p",       &reco_p);
	tOut->Branch("reco_trckmar", &reco_trckmar);
	tOut->Branch("reco_trcklin", &reco_trcklin);
	tOut->Branch("reco_chi2dof", &reco_chi2dof);
	tOut->Branch("reco_ismtchd", &reco_ismtchd);
	tOut->Branch("reco_idmtchd", &reco_idmtchd);
	// tOut->Branch("reco_cluster4_id", &reco_cluster4_id);
	// tOut->Branch("reco_cluster1_id", &reco_cluster1_id);
	// tOut->Branch("reco_cluster2_id", &reco_cluster2_id);
	// tOut->Branch("reco_cluster3_id", &reco_cluster3_id);
	
	/// monitoring histograms
	TH1D* h_dErel_rec_gen    = new TH1D("h_dErel_rec_gen","Rec vs Gen;(E_{rec}-E_{gen})/E_{gen};Tracks",500,-1.,+1.);
	TH1D* h_chi2             = new TH1D("h_chi2",";#chi^2;Tracks",50,0,15);
	TH1D* h_chi2_matched     = new TH1D("h_chi2_matched",";#chi^2;Tracks",50,0,15);
	TH1D* h_chi2_nonmatched  = new TH1D("h_chi2_nonmatched",";#chi^2;Tracks",50,0,15);
	TH1D* h_E_tru_all        = new TH1D("h_E_tru_all",";#it{E}_{tru}^{all} [GeV];Tracks",34,0,17);
	TH1D* h_E_tru_mat        = new TH1D("h_E_tru_mat",";#it{E}_{tru}^{mat} [GeV];Tracks",34,0,17);
	TH1D* h_E_eff            = new TH1D("h_E_eff",";#it{E}_{tru} [GeV];Tracks",34,0,17);
 
	/// prepare the dictionaries
	prepare_cahced_clusters();
 
	/// loop on events
	Int_t nsigevents = tSig->GetEntries();
	cout << "Starting loop over signal events with nsigevents=" << nsigevents << endl;
	for(int iev=0 ; iev<tSig->GetEntries() and iev<1 ; iev++)
		// for(int iev=0 ; iev<tSig->GetEntries() ; iev++)
	{
		/// clear output vectors
		// for(unsigned int x=0 ; x<all_clusters_xyz.size() ; ++x) delete all_clusters_xyz[x];
		// all_clusters_xyz.clear();
		// all_clusters_type.clear();
		// all_clusters_id.clear();
		
		true_q.clear();
		true_p.clear();
		// for(unsigned int x=0 ; x<true_trckmar.size() ; ++x) delete true_trckmar[x];
		// for(unsigned int x=0 ; x<true_trcklin.size() ; ++x) delete true_trcklin[x];
		true_trckmar.clear();
		true_trcklin.clear();
		// true_cluster4_id.clear();
		// true_cluster3_id.clear();
		// true_cluster2_id.clear();
		// true_cluster1_id.clear();
	
		for(unsigned int x=0 ; x<bkgr_trckmar.size() ; ++x) delete bkgr_trckmar[x];
		for(unsigned int x=0 ; x<bkgr_trcklin.size() ; ++x) delete bkgr_trcklin[x];
		bkgr_trckmar.clear();
		bkgr_trcklin.clear();
		// bkgr_cluster4_id.clear();
		// bkgr_cluster3_id.clear();
		// bkgr_cluster2_id.clear();
		// bkgr_cluster1_id.clear();
	
		seed_q.clear();
		seed_p.clear();
		// seed_cluster4_id.clear();
		// seed_cluster1_id.clear();
	
		reco_q.clear();
		reco_p.clear();
		for(unsigned int x=0 ; x<reco_trckmar.size() ; ++x) delete reco_trckmar[x];
		for(unsigned int x=0 ; x<reco_trcklin.size() ; ++x) delete reco_trcklin[x];
		reco_trckmar.clear();
		reco_trcklin.clear();
		reco_chi2dof.clear();
		reco_ismtchd.clear();
		reco_idmtchd.clear();
		// reco_cluster4_id.clear();
		// reco_cluster3_id.clear();
		// reco_cluster2_id.clear();
		// reco_cluster1_id.clear();
		
		//// clear cached clusters
		clear_cached_clusters(); /// clear for both sides
		
		/// rest all the layers of the detector (including inactive if any)
		reset_layers_all(); // reset both sides 
		
		//// get the next entry
		tSig->GetEntry(iev);
		tBkg->GetEntry(iev);
		
		if((iev%outN)==0) printf("Done %d out of %d\n",iev,nsigevents);
		
		for(unsigned int s=0 ; s<sides.size() ; ++s)
		{
			TString side = sides[s];
			
			/// set the charge
			float crg = (side=="Eside") ? -1 : +1;
		   
			/// globals (per side)
			unsigned int n_truth = 0;
			unsigned int n_seeds = 0;
			unsigned int n_solve = 0;
			unsigned int n_recos = 0;
			unsigned int n_match = 0;
		   
			/// make a pool of all signal clusters
			cache_signal_clusters(polm_clusters,side);
			
			// cout << "After signal:" << endl;
			// cout << "cached_clusters_xyz[x_L1_'+side+'].size()=" << cached_clusters_xyz["x_L1_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L2_'+side+'].size()=" << cached_clusters_xyz["x_L2_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L3_'+side+'].size()=" << cached_clusters_xyz["x_L3_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L4_'+side+'].size()=" << cached_clusters_xyz["x_L4_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L1_'+side+'].size()=" << cached_clusters_att["id_L1_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L2_'+side+'].size()=" << cached_clusters_att["id_L2_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L3_'+side+'].size()=" << cached_clusters_att["id_L3_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L4_'+side+'].size()=" << cached_clusters_att["id_L4_"+side].size() << endl;		
			
			/// make a pool of all background and noise clusters
			cache_background_clusters(clusters_xyz,clusters_type,clusters_id,side);
			
			// cout << "After background:" << endl;
			// cout << "cached_clusters_xyz[x_L1_'+side+'].size()=" << cached_clusters_xyz["x_L1_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L2_'+side+'].size()=" << cached_clusters_xyz["x_L2_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L3_'+side+'].size()=" << cached_clusters_xyz["x_L3_"+side].size() << endl;
			// cout << "cached_clusters_xyz[x_L4_'+side+'].size()=" << cached_clusters_xyz["x_L4_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L1_'+side+'].size()=" << cached_clusters_att["id_L1_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L2_'+side+'].size()=" << cached_clusters_att["id_L2_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L3_'+side+'].size()=" << cached_clusters_att["id_L3_"+side].size() << endl;
			// cout << "cached_clusters_att[id_L4_'+side+'].size()=" << cached_clusters_att["id_L4_"+side].size() << endl;
		   
			/// rest all the layers of the detector (including inactive if any)
			reset_layers_all(); // reset both sides 
			
			/// add all clusters to the detector (Eside/Pside/both)
			add_all_clusters(side);
			
			// /// write out all clusters
			// write_out_clusters("Eside");
		   
			/// truth tracks:
			for(unsigned int t=0 ; t<qgen->size() ; ++t)
			{
				if(side=="Eside" and qgen->at(t)>0) continue;
				if(side=="Pside" and qgen->at(t)<0) continue;
				h_E_tru_all->Fill( pgen->at(t).E() );
				n_truth++;
				true_q.push_back( qgen->at(t) );
				true_p.push_back( pgen->at(t) );
				true_trckmar.push_back( polm_gen->at(t) );
				true_trcklin.push_back( poll_gen->at(t) );
				// true_cluster1_id.push_back( t );
				// true_cluster2_id.push_back( t );
				// true_cluster3_id.push_back( t );
				// true_cluster4_id.push_back( t );
			}
			
			// /// background tracks (allways appear )
			// for(unsigned int b=0 ; b<trkid1Cluster->size() ; ++b)
			// {
			//    bkgr_trckmar.push_back(  );
			//    bkgr_trcklin.push_back(  );
			//    bkgr_cluster4_id.push_back( trkid4Cluster->at(b) );
			//    bkgr_cluster3_id.push_back( trkid3Cluster->at(b) );
			//    bkgr_cluster2_id.push_back( trkid2Cluster->at(b) );
			//    bkgr_cluster1_id.push_back( trkid1Cluster->at(b) );
			// 	   }
	      
		   
			/// run over all clusters of layer 4 in the pool --> these are the seed for the KalmanFilter fit
			for(unsigned int i4=0 ; i4<cached_clusters_xyz["x_L4_"+side].size() ; ++i4)
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
					float r1[3] = {cached_clusters_xyz["x_L1_"+side][i1], cached_clusters_xyz["y_L1_"+side][i1], cached_clusters_xyz["z_L1_"+side][i1]};
					float r4[3] = {cached_clusters_xyz["x_L4_"+side][i4], cached_clusters_xyz["y_L4_"+side][i4], cached_clusters_xyz["z_L4_"+side][i4]};
					bool seed = makeseed(process,r1,r4,i1,i4,side,pseed);
					if(!seed) continue; // cannot make a meaningful seed
					pseeds.push_back(pseed);
					n_seeds++;
		   		
					seed_q.push_back(crg);
					seed_p.push_back(pseed);
				} // end of loop on clusters in layer 1
				if(n_seeds<1) continue;
		   	
				// prepare the probe from the seed and do the KF fit
				bool solved = det->SolveSingleTrackViaKalmanMC_Noam_multiseed(pseeds,meGeV,crg,99);
				if(!solved) continue; // reconstruction failed
				n_solve++;
				
				// cout << "For " << side << ": layer 4 is filled with " << det->GetLayer(7)->GetNBgClusters() << " clusters" << endl;
				
				// get the reconstructed propagated to the vertex 
				KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack(); 
				if(!trw) continue; // track was not reconstructed
				n_recos++;
				
				/// get the clusters of the winner track
				int win_cls_id1 = trw->GetClID(1);
				int win_cls_id2 = trw->GetClID(3);
				int win_cls_id3 = trw->GetClID(5);
				int win_cls_id4 = trw->GetClID(7);
				float r1_ntrksys[3] = { cached_clusters_xyz["x_L1_"+side][win_cls_id1],cached_clusters_xyz["y_L1_"+side][win_cls_id1], cached_clusters_xyz["z_L1_"+side][win_cls_id1] };
				float r2_ntrksys[3] = { cached_clusters_xyz["x_L2_"+side][win_cls_id2],cached_clusters_xyz["y_L2_"+side][win_cls_id2], cached_clusters_xyz["z_L2_"+side][win_cls_id2] };
				float r3_ntrksys[3] = { cached_clusters_xyz["x_L3_"+side][win_cls_id3],cached_clusters_xyz["y_L3_"+side][win_cls_id3], cached_clusters_xyz["z_L3_"+side][win_cls_id3] };
				float r4_ntrksys[3] = { cached_clusters_xyz["x_L4_"+side][win_cls_id4],cached_clusters_xyz["y_L4_"+side][win_cls_id4], cached_clusters_xyz["z_L4_"+side][win_cls_id4] };
								
				TLorentzVector prec;
				double pxyz[3];
				trw->GetPXYZ(pxyz);
				prec.SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
				
				cout << "Etru=" << pgen->at(cached_clusters_att["id_L4_"+side][i4]).E() << " GeV, Erec=" << prec.E() << "GeV --> i4=" << i4 << ": win_cls_id1=" << win_cls_id1 << ", win_cls_id2=" << win_cls_id2 << ", win_cls_id3=" << win_cls_id3 << ", win_cls_id4=" << win_cls_id4 << endl;
				
				int ismatched = cached_clusters_att["type_L4_"+side][i4]; // TODO: NEED TO GET THE WINNER CLUSTERS AND CHECK ALL LAYERS
				int idmatched = cached_clusters_att["id_L4_"+side][i4];
				float chi2dof = trw->GetNormChi2();
				
				/// fill the reco branches
				reco_q.push_back( crg );
				reco_p.push_back( prec );
				reco_trckmar.push_back( TrackMarker3d(trw,0,361,1,trkcol(prec.E())) );
				reco_trcklin.push_back( TrackLine3d(trw,361,1,trkcol(prec.E())) );
				reco_chi2dof.push_back( chi2dof );
				reco_ismtchd.push_back( ismatched );
				reco_idmtchd.push_back( idmatched );
				
				h_chi2->Fill( chi2dof );
				if(ismatched==1 and idmatched>=0) // i4 is a signal cluster
				{
					bool accept = (chi2dof<3. and (prec.E()>1. and prec.E()<18.));
					if(accept)
					{
						h_E_tru_mat->Fill( pgen->at(idmatched).E() );
						h_dErel_rec_gen->Fill( (prec.E()-pgen->at(idmatched).E())/pgen->at(idmatched).E() );
					}
					h_chi2_matched->Fill( chi2dof );
					n_match++;
				}
				else h_chi2_nonmatched->Fill( chi2dof );
				
				/// this is maybe redundant
				pseeds.clear();
				
			} // end of loop on clusters in layer 4

			cout << "Event #" << iev << ", "<< side << ": n_truth=" << n_truth
				<< ", n_seeds=" << n_seeds
					<< ", n_solve=" << n_solve
						<< ", n_recos=" << n_recos
							<< ", n_match=" << n_match << endl;
		} // end of loop on sides
		fOut->cd();
		tOut->Fill();
	}
	
	h_E_eff->Divide(h_E_tru_mat,h_E_tru_all);
	fOut->cd();
	tOut->Write();
	h_chi2->Write();
	h_chi2_matched->Write();
	h_chi2_nonmatched->Write();
	h_dErel_rec_gen->Write();
	fOut->Write();
	fOut->Close();
}
