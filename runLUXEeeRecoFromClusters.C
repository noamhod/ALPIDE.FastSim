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
// #include "TIter.h"
#endif

typedef map<TString, int >             TMapTSi;
typedef map<TString, float >           TMapTSf;
typedef map<TString, double >          TMapTSd;
typedef map<TString, vector<int>* >    TMapTSP2vi;
typedef map<TString, vector<float>* >  TMapTSP2vf;
typedef map<TString, vector<double>* > TMapTSP2vd;
typedef map<TString, TH1*>             TMapTSP2TH1;
typedef map<TString, TH2*>             TMapTSP2TH2;


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

//// dipole geometry
double xW = 120;
double yH = 67.2;
double z1 = 102.9;
double z2 = 202.9;

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

vector<float> cls_x_L1_Eside, cls_y_L1_Eside, cls_z_L1_Eside, cls_id_L1_Eside, cls_type_L1_Eside;
vector<float> cls_x_L2_Eside, cls_y_L2_Eside, cls_z_L2_Eside, cls_id_L2_Eside, cls_type_L2_Eside;
vector<float> cls_x_L3_Eside, cls_y_L3_Eside, cls_z_L3_Eside, cls_id_L3_Eside, cls_type_L3_Eside;
vector<float> cls_x_L4_Eside, cls_y_L4_Eside, cls_z_L4_Eside, cls_id_L4_Eside, cls_type_L4_Eside;

vector<float> cls_x_L1_Pside, cls_y_L1_Pside, cls_z_L1_Pside, cls_id_L1_Pside, cls_type_L1_Pside;
vector<float> cls_x_L2_Pside, cls_y_L2_Pside, cls_z_L2_Pside, cls_id_L2_Pside, cls_type_L2_Pside;
vector<float> cls_x_L3_Pside, cls_y_L3_Pside, cls_z_L3_Pside, cls_id_L3_Pside, cls_type_L3_Pside;
vector<float> cls_x_L4_Pside, cls_y_L4_Pside, cls_z_L4_Pside, cls_id_L4_Pside, cls_type_L4_Pside;

void clear_clusters()
{
	cls_x_L1_Eside.clear();
	cls_x_L2_Eside.clear();
	cls_x_L3_Eside.clear();
	cls_x_L4_Eside.clear();
	
	cls_y_L1_Eside.clear();
	cls_y_L2_Eside.clear();
	cls_y_L3_Eside.clear();
	cls_y_L4_Eside.clear();
	
	cls_z_L1_Eside.clear();
	cls_z_L2_Eside.clear();
	cls_z_L3_Eside.clear();
	cls_z_L4_Eside.clear();
	
	cls_id_L1_Eside.clear();
	cls_id_L2_Eside.clear();
	cls_id_L3_Eside.clear();
	cls_id_L4_Eside.clear();
	
	cls_type_L1_Eside.clear();
	cls_type_L2_Eside.clear();
	cls_type_L3_Eside.clear();
	cls_type_L4_Eside.clear();
	
	
	cls_x_L1_Pside.clear();
	cls_x_L2_Pside.clear();
	cls_x_L3_Pside.clear();
	cls_x_L4_Pside.clear();
	
	cls_y_L1_Pside.clear();
	cls_y_L2_Pside.clear();
	cls_y_L3_Pside.clear();
	cls_y_L4_Pside.clear();
	
	cls_z_L1_Pside.clear();
	cls_z_L2_Pside.clear();
	cls_z_L3_Pside.clear();
	cls_z_L4_Pside.clear();
	
	cls_id_L1_Pside.clear();
	cls_id_L2_Pside.clear();
	cls_id_L3_Pside.clear();
	cls_id_L4_Pside.clear();
	
	cls_type_L1_Pside.clear();
	cls_type_L2_Pside.clear();
	cls_type_L3_Pside.clear();
	cls_type_L4_Pside.clear();
}

void cache_signal_clusters(vector<TPolyMarker3D*>* polm_clusters, TString side)
{
	/// fill signam clusters
	for(unsigned int i=0 ; i<polm_clusters->size() ; i++)
	{
		for(Int_t j=0 ; j<polm_clusters->at(i)->GetN() ; ++j)
		{
			float x,y,z;
			polm_clusters->at(i)->GetPoint(j,x,y,z); // the clusters
			cout << "signal xyz[" << i << "] " << x << "," << y << ", " << z << endl;
			if(x>0 and (side=="Eside" or side=="both")) /// Eside
			{	
				if(z==300) { cls_x_L1_Eside.push_back(x); cls_y_L1_Eside.push_back(y); cls_z_L1_Eside.push_back(z); cls_id_L1_Eside.push_back(i); cls_type_L1_Eside.push_back(1); }
				if(z==310) { cls_x_L2_Eside.push_back(x); cls_y_L2_Eside.push_back(y); cls_z_L2_Eside.push_back(z); cls_id_L2_Eside.push_back(i); cls_type_L2_Eside.push_back(1); }
				if(z==320) { cls_x_L3_Eside.push_back(x); cls_y_L3_Eside.push_back(y); cls_z_L3_Eside.push_back(z); cls_id_L3_Eside.push_back(i); cls_type_L3_Eside.push_back(1); }
				if(z==330) { cls_x_L4_Eside.push_back(x); cls_y_L4_Eside.push_back(y); cls_z_L4_Eside.push_back(z); cls_id_L4_Eside.push_back(i); cls_type_L4_Eside.push_back(1); }
			}
			if(x<0 and (side=="Pside" or side=="both")) /// Pside
			{
				if(z==300) { cls_x_L1_Pside.push_back(x); cls_y_L1_Pside.push_back(y); cls_z_L1_Pside.push_back(z); cls_id_L1_Pside.push_back(i); cls_type_L1_Pside.push_back(1); }
				if(z==310) { cls_x_L2_Pside.push_back(x); cls_y_L2_Pside.push_back(y); cls_z_L2_Pside.push_back(z); cls_id_L2_Pside.push_back(i); cls_type_L2_Pside.push_back(1); }
				if(z==320) { cls_x_L3_Pside.push_back(x); cls_y_L3_Pside.push_back(y); cls_z_L3_Pside.push_back(z); cls_id_L3_Pside.push_back(i); cls_type_L3_Pside.push_back(1); }
				if(z==330) { cls_x_L4_Pside.push_back(x); cls_y_L4_Pside.push_back(y); cls_z_L4_Pside.push_back(z); cls_id_L4_Pside.push_back(i); cls_type_L4_Pside.push_back(1); }
			}
		}
	}
}

void cache_background_clusters(
	vector<float> *x1Cluster, vector<float> *y1Cluster, vector<float> *z1Cluster, vector<int> *type1Cluster, vector<int> *trkid1Cluster,
	vector<float> *x2Cluster, vector<float> *y2Cluster, vector<float> *z2Cluster, vector<int> *type2Cluster, vector<int> *trkid2Cluster,
	vector<float> *x3Cluster, vector<float> *y3Cluster, vector<float> *z3Cluster, vector<int> *type3Cluster, vector<int> *trkid3Cluster,
	vector<float> *x4Cluster, vector<float> *y4Cluster, vector<float> *z4Cluster, vector<int> *type4Cluster, vector<int> *trkid4Cluster,
	TString side)
{
	for(unsigned int i=0 ; i<x1Cluster->size() ; i++)
	{
		if(x1Cluster->at(i)>0 and (side=="Eside" or side=="both")) /// Eside
		{
			cls_x_L1_Eside.push_back(x1Cluster->at(i)); cls_y_L1_Eside.push_back(y1Cluster->at(i)); cls_z_L1_Eside.push_back(z1Cluster->at(i)); cls_id_L1_Eside.push_back(trkid1Cluster->at(i)); cls_type_L1_Eside.push_back(type1Cluster->at(i));
			cls_x_L2_Eside.push_back(x2Cluster->at(i)); cls_y_L2_Eside.push_back(y2Cluster->at(i)); cls_z_L2_Eside.push_back(z2Cluster->at(i)); cls_id_L2_Eside.push_back(trkid2Cluster->at(i)); cls_type_L2_Eside.push_back(type2Cluster->at(i));
			cls_x_L3_Eside.push_back(x3Cluster->at(i)); cls_y_L3_Eside.push_back(y3Cluster->at(i)); cls_z_L3_Eside.push_back(z3Cluster->at(i)); cls_id_L3_Eside.push_back(trkid3Cluster->at(i)); cls_type_L3_Eside.push_back(type3Cluster->at(i));
			cls_x_L4_Eside.push_back(x4Cluster->at(i)); cls_y_L4_Eside.push_back(y4Cluster->at(i)); cls_z_L4_Eside.push_back(z4Cluster->at(i)); cls_id_L4_Eside.push_back(trkid4Cluster->at(i)); cls_type_L4_Eside.push_back(type4Cluster->at(i));
		}
		if(x1Cluster->at(i)<0 and (side=="Pside" or side=="both")) /// Pside
		{
			cls_x_L1_Pside.push_back(x1Cluster->at(i)); cls_y_L1_Pside.push_back(y1Cluster->at(i)); cls_z_L1_Pside.push_back(z1Cluster->at(i)); cls_id_L1_Pside.push_back(trkid1Cluster->at(i)); cls_type_L1_Pside.push_back(type1Cluster->at(i));
			cls_x_L2_Pside.push_back(x2Cluster->at(i)); cls_y_L2_Pside.push_back(y2Cluster->at(i)); cls_z_L2_Pside.push_back(z2Cluster->at(i)); cls_id_L2_Pside.push_back(trkid2Cluster->at(i)); cls_type_L2_Pside.push_back(type2Cluster->at(i));
			cls_x_L3_Pside.push_back(x3Cluster->at(i)); cls_y_L3_Pside.push_back(y3Cluster->at(i)); cls_z_L3_Pside.push_back(z3Cluster->at(i)); cls_id_L3_Pside.push_back(trkid3Cluster->at(i)); cls_type_L3_Pside.push_back(type3Cluster->at(i));
			cls_x_L4_Pside.push_back(x4Cluster->at(i)); cls_y_L4_Pside.push_back(y4Cluster->at(i)); cls_z_L4_Pside.push_back(z4Cluster->at(i)); cls_id_L4_Pside.push_back(trkid4Cluster->at(i)); cls_type_L4_Pside.push_back(type4Cluster->at(i));
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
	if(side=="Eside" or side=="both")
	{
   	for(unsigned int i1=0 ; i1<cls_x_L1_Eside.size() ; ++i1) add_bkg_cluster(1,cls_x_L1_Eside[i1],cls_y_L1_Eside[i1],cls_z_L1_Eside[i1],i1);
	   for(unsigned int i2=0 ; i2<cls_x_L2_Eside.size() ; ++i2) add_bkg_cluster(3,cls_x_L2_Eside[i2],cls_y_L2_Eside[i2],cls_z_L2_Eside[i2],i2);
	   for(unsigned int i3=0 ; i3<cls_x_L3_Eside.size() ; ++i3) add_bkg_cluster(5,cls_x_L3_Eside[i3],cls_y_L3_Eside[i3],cls_z_L3_Eside[i3],i3);
	   for(unsigned int i4=0 ; i4<cls_x_L4_Eside.size() ; ++i4) add_bkg_cluster(7,cls_x_L4_Eside[i4],cls_y_L4_Eside[i4],cls_z_L4_Eside[i4],i4);
	}
	if(side=="Pside" or side=="both")
	{
	   for(unsigned int i1=0 ; i1<cls_x_L1_Pside.size() ; ++i1) add_bkg_cluster(1,cls_x_L1_Pside[i1],cls_y_L1_Pside[i1],cls_z_L1_Pside[i1],i1);
	   for(unsigned int i2=0 ; i2<cls_x_L2_Pside.size() ; ++i2) add_bkg_cluster(3,cls_x_L2_Pside[i2],cls_y_L2_Pside[i2],cls_z_L2_Pside[i2],i2);
	   for(unsigned int i3=0 ; i3<cls_x_L3_Pside.size() ; ++i3) add_bkg_cluster(5,cls_x_L3_Pside[i3],cls_y_L3_Pside[i3],cls_z_L3_Pside[i3],i3);
	   for(unsigned int i4=0 ; i4<cls_x_L4_Pside.size() ; ++i4) add_bkg_cluster(7,cls_x_L4_Pside[i4],cls_y_L4_Pside[i4],cls_z_L4_Pside[i4],i4);
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
		cout << "ERROR in xofz: dx=0" << endl;
		exit(-1);
	}
	float a = dz/dx;
	float b = r1[2]-a*r1[0];
	float z = a*x+b;
	return z;
}

bool makeseed(float* r1, float* r4, TString side, TLorentzVector& p)
{
	if(side=="Eside" and r1[0]>=r4[0]) return false;
	if(side=="Pside" and r1[0]<=r4[0]) return false;
	if(r1[2]==r4[2])                   return false;

	TRandom rnd;
	rnd.SetSeed();
	double posneg = rnd.Uniform(-1,+1);
	double pxgaus = rnd.Gaus(7.2e-4,5.e-4);

	double x0 = 0;
	double z0 = zofx(r1,r4,x0);
	double xExit = xofz(r1,r4,zDipoleExit)*cm2m;
	double H = (zDipoleExit-z0)*cm2m;
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
	det->SetMaxChi2Vtx(20e9);  // fiducial cut on chi2 of convergence to vtx
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
	vector<TPolyMarker3D*>*  polm_clusters_intrksys = 0;
	tSig->SetBranchAddress("wgtgen",                 &wgtgen);
	tSig->SetBranchAddress("qgen",                   &qgen);
	tSig->SetBranchAddress("pgen",                   &pgen);
	tSig->SetBranchAddress("igentrk",                &igentrk);
	tSig->SetBranchAddress("acctrkgen",              &acctrkgen);
	tSig->SetBranchAddress("polm_clusters",          &polm_clusters);
	tSig->SetBranchAddress("polm_clusters_intrksys", &polm_clusters_intrksys);
	tSig->SetBranchAddress("polm_gen",               &polm_gen);
	tSig->SetBranchAddress("poll_gen",               &poll_gen);
	
	/// get the background and noise clusters
	cout << "Getting background clusters from tree" << endl;
	TFile* fBkg = new TFile("data/root/background_clusters_"+process+".root","READ");
	TTree* tBkg = (TTree*)fBkg->Get("clusters");
	vector<int>     *type1Cluster = 0;
	vector<int>     *type2Cluster = 0;
	vector<int>     *type3Cluster = 0;
	vector<int>     *type4Cluster = 0;
	vector<int>     *trkid1Cluster = 0;
	vector<int>     *trkid2Cluster = 0;
	vector<int>     *trkid3Cluster = 0;
	vector<int>     *trkid4Cluster = 0;
	vector<float>   *x1Cluster = 0;
	vector<float>   *y1Cluster = 0;
	vector<float>   *z1Cluster = 0;
	vector<float>   *x2Cluster = 0;
	vector<float>   *y2Cluster = 0;
	vector<float>   *z2Cluster = 0;
	vector<float>   *x3Cluster = 0;
	vector<float>   *y3Cluster = 0;
	vector<float>   *z3Cluster = 0;
	vector<float>   *x4Cluster = 0;
	vector<float>   *y4Cluster = 0;
	vector<float>   *z4Cluster = 0;
	tBkg->SetBranchAddress("type1Cluster",  &type1Cluster);
	tBkg->SetBranchAddress("type2Cluster",  &type2Cluster);
	tBkg->SetBranchAddress("type3Cluster",  &type3Cluster);
	tBkg->SetBranchAddress("type4Cluster",  &type4Cluster);
	tBkg->SetBranchAddress("trkid1Cluster", &trkid1Cluster);
	tBkg->SetBranchAddress("trkid2Cluster", &trkid2Cluster);
	tBkg->SetBranchAddress("trkid3Cluster", &trkid3Cluster);
	tBkg->SetBranchAddress("trkid4Cluster", &trkid4Cluster);
	tBkg->SetBranchAddress("x1Cluster",     &x1Cluster);
	tBkg->SetBranchAddress("y1Cluster",     &y1Cluster);
	tBkg->SetBranchAddress("z1Cluster",     &z1Cluster);
	tBkg->SetBranchAddress("x2Cluster",     &x2Cluster);
	tBkg->SetBranchAddress("y2Cluster",     &y2Cluster);
	tBkg->SetBranchAddress("z2Cluster",     &z2Cluster);
	tBkg->SetBranchAddress("x3Cluster",     &x3Cluster);
	tBkg->SetBranchAddress("y3Cluster",     &y3Cluster);
	tBkg->SetBranchAddress("z3Cluster",     &z3Cluster);
	tBkg->SetBranchAddress("x4Cluster",     &x4Cluster);
	tBkg->SetBranchAddress("y4Cluster",     &y4Cluster);
	tBkg->SetBranchAddress("z4Cluster",     &z4Cluster);


	// output tree
	cout << "Setting the output tree" << endl;
	gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
	gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
	TFile* fOut = new TFile("data/root/rec_from_clusters_"+process+".root","RECREATE");
	TTree* tOut = new TTree("res","res");
	int                      n_gen = 0;
	vector<int>              q_gen;
	vector<TLorentzVector>   p_gen;
	int                      n_rec = 0;
	vector<int>              q_rec;
	vector<TLorentzVector>   p_rec;
	vector<int>              acc_rec;
	vector<int>              jgen_rec;
	tOut->Branch("n_gen",    &n_gen);
	tOut->Branch("q_gen",    &q_gen);
	tOut->Branch("p_gen",    &p_gen);
	tOut->Branch("n_rec",    &n_rec);
	tOut->Branch("jgen_rec", &jgen_rec);
	tOut->Branch("q_rec",    &q_rec);
	tOut->Branch("p_rec",    &p_rec);
	tOut->Branch("acc_rec",  &acc_rec);
	
   TH1D* h_dErel_rec_gen = new TH1D("h_dErel_rec_gen","Rec vs Gen;(E_{rec}-E_{gen})/E_{gen};Tracks",500,-1.,+1.);
   TH1D* h_chi2          = new TH1D("h_chi2",";#chi^2;Tracks",100,0,25);
   TH1D* h_chi2_matched  = new TH1D("h_chi2_matched",";#chi^2;Tracks",100,0,25);
 
	/// loop on events
	Int_t nsigevents = tSig->GetEntries();
	cout << "Starting loop over signal events with nsigevents=" << nsigevents << endl;
	for(int iev=0 ; iev<tSig->GetEntries() and iev<1 ; iev++)
	{
		/// global
		TLorentzVector ptmp;
		
		//// clear
		clear_clusters();
      
		//// get the next entry
		tSig->GetEntry(iev);
		tBkg->GetEntry(iev);
		if((iev%outN)==0) printf("Done %d out of %d\n",iev,nsigevents);
		
		/// make a pool of all signal clusters
		cache_signal_clusters(polm_clusters, "Eside");
		
		/// make a pool of all background and noise clusters
		cache_background_clusters(
			x1Cluster,y1Cluster,z1Cluster,type1Cluster,trkid1Cluster,
         x2Cluster,y2Cluster,z2Cluster,type2Cluster,trkid2Cluster,
         x3Cluster,y3Cluster,z3Cluster,type3Cluster,trkid3Cluster,
         x4Cluster,y4Cluster,z4Cluster,type4Cluster,trkid4Cluster,
			"Eside"
		);
		
		/// rest all the layers of the detector (including inactive if any):
		reset_layers_all();
		
		/// add all clusters to the detector (Eside/Pside/both)
		add_all_clusters("Eside");
		
		/// run over all clusters of layer 4 in the pool --> these are the seed for the KalmanFilter fit
		float crg = -1; // Eside...
		// for(unsigned int i4=0 ; i4<cls_x_L4_Eside.size() ; ++i4)
		for(unsigned int i4=0 ; i4<1 ; ++i4)
		{
			int n_seeds = 0;
			// cout << "starting test of i4=" << i4 << " with N1=" << cls_x_L1_Eside.size() << endl;
			reset_layers_tracks(); // reset all tracks from all layers
			cout << "All seeds for i4=" << i4 << " (type="<< (cls_type_L4_Eside[i4]==1)<< ", itru=" << cls_id_L4_Eside[i4] << ")" << ":" << endl;
			vector<TLorentzVector> pseeds;
			// for(unsigned int i1=0 ; i1<cls_x_L1_Eside.size() ; ++i1)
			for(unsigned int i1=0 ; i1<1 ; ++i1)
			{
				reset_layers_tracks(0); // reset all tracks from all layers but layer 0
				/// find the momentum of the seed
				TLorentzVector pseed;
				float r1[3] = {cls_x_L1_Eside[i1], cls_y_L1_Eside[i1], cls_z_L1_Eside[i1]};
				float r4[3] = {cls_x_L4_Eside[i4], cls_y_L4_Eside[i4], cls_z_L4_Eside[i4]};
				bool seed = makeseed(r1,r4,"Eside",pseed);
				if(!seed) continue; // cannot make a meaningful seed
				pseeds.push_back(pseed);
				n_seeds++;
			}
			
			// prepare the probe from the seed and do the KF fit
		   bool solved = det->SolveSingleTrackViaKalmanMC_Noam_multiseed(pseeds,meGeV,crg,99);
			if(!solved) continue; // reconstruction failed
			cout << "Solved! (i4=" << i4 << ")" << endl;
			
			// get the reconstructed propagated to the vertex 
			KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack(); 
			if(!trw) continue; // track was not reconstructed
			cout << "Reconstructed! (i4=" << i4 << ")" << endl;
			
			TLorentzVector prec;
			double pxyz[3]; trw->GetPXYZ(pxyz);
			prec.SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
			Double_t normChi2    = trw->GetNormChi2();
			h_chi2->Fill(normChi2);
			cout << " - i4=" << i4 << " --> Passed with normChi2=" << normChi2 << endl;
			if(cls_type_L4_Eside[i4]==1) // it is a signal cluster
			{
				int itru = cls_id_L4_Eside[i4];
				if(itru>=0) cout << "         itru=" << itru << " --> Etru=" << pgen->at(itru).E() << ", Erec=" << prec.E() << endl;
				h_dErel_rec_gen->Fill( (prec.E()-pgen->at(itru).E())/pgen->at(itru).E() );
				h_chi2_matched->Fill(normChi2);
			}
			pseeds.clear();
			cout << "End of all seeds for i4=" << i4 << "\n\n" << endl;
		}
		fOut->cd();
		tOut->Fill();
	}
	fOut->cd();
	tOut->Write();
	h_chi2->Write();
	h_chi2_matched->Write();
	h_dErel_rec_gen->Write();
	fOut->Write();
	fOut->Close();
}
