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

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 27;   // cm for BPPP or 50 for Trident
double Rbeampipe = 2.413; // cm
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


void ReconstructionFromSeeds(TString process, int Seed=12345) //, const char* setup="setup/setupLUXE.txt")
{
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
	// for reconstruction:
	det->SetErrorScale(200.);
   det->Print();
   // det->BookControlHistos();
   
   zlayer->push_back(0); //// NOAM --> GET FROM THE SETUP   --> IP (vertex)
   zlayer->push_back(100); //// NOAM --> GET FROM THE SETUP --> start of dipol
   zlayer->push_back(200); //// NOAM --> GET FROM THE SETUP --> end of dipol
   zlayer->push_back(300); //// NOAM --> GET FROM THE SETUP
   zlayer->push_back(310); //// NOAM --> GET FROM THE SETUP
   zlayer->push_back(320); //// NOAM --> GET FROM THE SETUP
   zlayer->push_back(330); //// NOAM --> GET FROM THE SETUP

   int outN = (process=="trident") ? 10 : 100;

   // TString process = "bppp";  /// trident or bppp or bppp_bkg or trident_bkg
	if(process=="trident") resetToTridentGeometry();

   /// get the particles from a ttree
   TFile* fIn = new TFile("data/root/seeds_"+process+".root","READ");
   TTree* tIn = (TTree*)fIn->Get("seeds");
   int nev = tIn->GetEntries();
	
   vector<float>   *svd0Seed = 0;
   vector<float>   *svd1Seed = 0;
   vector<float>   *svd2Seed = 0;
   vector<float>   *chi2xzSeed = 0;
   vector<float>   *chi2yzSeed = 0;
   vector<float>   *residxzSeed = 0;
   vector<float>   *residyzSeed = 0;
   vector<int>     *issigSeed = 0;
   vector<int>     *iGenMatch = 0;
   vector<float>   *x1Seed = 0;
   vector<float>   *y1Seed = 0;
   vector<float>   *z1Seed = 0;
   vector<float>   *x2Seed = 0;
   vector<float>   *y2Seed = 0;
   vector<float>   *z2Seed = 0;
   vector<float>   *x3Seed = 0;
   vector<float>   *y3Seed = 0;
   vector<float>   *z3Seed = 0;
   vector<float>   *x4Seed = 0;
   vector<float>   *y4Seed = 0;
   vector<float>   *z4Seed = 0;
   vector<float>   *pxSeed = 0;
   vector<float>   *pySeed = 0;
   vector<float>   *pzSeed = 0;
   vector<float>   *eSeed = 0;
   vector<float>   *pxGen = 0;
   vector<float>   *pyGen = 0;
   vector<float>   *pzGen = 0;
   vector<float>   *eGen = 0;
   vector<float>   *qGen = 0;
   vector<int>     *iGen = 0;
	
   tIn->SetBranchAddress("svd0Seed",    &svd0Seed);
   tIn->SetBranchAddress("svd1Seed",    &svd1Seed);
   tIn->SetBranchAddress("svd2Seed",    &svd2Seed);
   tIn->SetBranchAddress("chi2xzSeed",  &chi2xzSeed);
   tIn->SetBranchAddress("chi2yzSeed",  &chi2yzSeed);
   tIn->SetBranchAddress("residxzSeed", &residxzSeed);
   tIn->SetBranchAddress("residyzSeed", &residyzSeed);
   tIn->SetBranchAddress("issigSeed",   &issigSeed);
   tIn->SetBranchAddress("iGenMatch",   &iGenMatch);
   tIn->SetBranchAddress("x1Seed",      &x1Seed);
   tIn->SetBranchAddress("y1Seed",      &y1Seed);
   tIn->SetBranchAddress("z1Seed",      &z1Seed);
   tIn->SetBranchAddress("x2Seed",      &x2Seed);
   tIn->SetBranchAddress("y2Seed",      &y2Seed);
   tIn->SetBranchAddress("z2Seed",      &z2Seed);
   tIn->SetBranchAddress("x3Seed",      &x3Seed);
   tIn->SetBranchAddress("y3Seed",      &y3Seed);
   tIn->SetBranchAddress("z3Seed",      &z3Seed);
   tIn->SetBranchAddress("x4Seed",      &x4Seed);
   tIn->SetBranchAddress("y4Seed",      &y4Seed);
   tIn->SetBranchAddress("z4Seed",      &z4Seed);
   tIn->SetBranchAddress("pxSeed",      &pxSeed);
   tIn->SetBranchAddress("pySeed",      &pySeed);
   tIn->SetBranchAddress("pzSeed",      &pzSeed);
   tIn->SetBranchAddress("eSeed",       &eSeed);
   tIn->SetBranchAddress("pxGen",       &pxGen);
   tIn->SetBranchAddress("pyGen",       &pyGen);
   tIn->SetBranchAddress("pzGen",       &pzGen);
   tIn->SetBranchAddress("eGen",        &eGen);
   tIn->SetBranchAddress("qGen",        &qGen);
   tIn->SetBranchAddress("iGen",        &iGen);
	
   // output tree
   gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");

   int                     n_gen = 0;
	vector<int>             j_gen;
	vector<int>             q_gen;
   vector<TLorentzVector>  p_gen;
	
   int                     n_seed = 0;
   vector<int>             jgen_seed;
   vector<int>             q_seed;
   vector<TLorentzVector>  p_seed;
   vector<int>             acc_seed;
   vector<TPolyMarker3D*>  polm_seed;
   vector<TPolyLine3D*>    poll_seed;
   vector<float>           svd0_seed;
   vector<float>           svd1_seed;
   vector<float>           svd2_seed;
   vector<float>           chi2xz_seed;
   vector<float>           chi2yz_seed;
   vector<float>           residxz_seed;
   vector<float>           residyz_seed;
   vector<int>             issig_seed;
   vector<float>           x1_seed;
   vector<float>           y1_seed;
   vector<float>           z1_seed;
   vector<float>           x2_seed;
   vector<float>           y2_seed;
   vector<float>           z2_seed;
   vector<float>           x3_seed;
   vector<float>           y3_seed;
   vector<float>           z3_seed;
   vector<float>           x4_seed;
   vector<float>           y4_seed;
   vector<float>           z4_seed;
	

   int                     n_rec = 0;
   vector<int>             jgen_rec;
   vector<int>             jseed_rec;
   vector<int>             q_rec;
   vector<TLorentzVector>  p_rec;
   vector<int>             acc_rec;
   vector<TPolyMarker3D*>  polm_rec;
   vector<TPolyLine3D*>    poll_rec;
   vector<int>             issig_rec;

   TFile* fOut = new TFile("data/root/rec_from_seeds_"+process+".root","RECREATE");
   TTree* tOut = new TTree("res","res");
   tOut->Branch("n_gen",    &n_gen);
   tOut->Branch("j_gen",    &j_gen);
	tOut->Branch("q_gen",    &q_gen);
   tOut->Branch("p_gen",    &p_gen);
	
   tOut->Branch("n_seed",      &n_seed);
	tOut->Branch("jgen_seed",   &jgen_seed);
	tOut->Branch("q_seed",      &q_seed);
   tOut->Branch("p_seed",      &p_seed);
   tOut->Branch("acc_seed",    &acc_seed);
   tOut->Branch("polm_seed",   &polm_seed);
   tOut->Branch("poll_seed",   &poll_seed);
	tOut->Branch("svd0_seed",    &svd0_seed);
	tOut->Branch("svd1_seed",    &svd1_seed);
	tOut->Branch("svd2_seed",    &svd2_seed);
	tOut->Branch("chi2xz_seed",  &chi2xz_seed);
	tOut->Branch("chi2yz_seed",  &chi2yz_seed);
	tOut->Branch("residxz_seed", &residxz_seed);
	tOut->Branch("residyz_seed", &residyz_seed);
	tOut->Branch("issig_seed",   &issig_seed);
	tOut->Branch("x1_seed",      &x1_seed);
	tOut->Branch("y1_seed",      &y1_seed);
	tOut->Branch("z1_seed",      &z1_seed);
	tOut->Branch("x2_seed",      &x2_seed);
	tOut->Branch("y2_seed",      &y2_seed);
	tOut->Branch("z2_seed",      &z2_seed);
	tOut->Branch("x3_seed",      &x3_seed);
	tOut->Branch("y3_seed",      &y3_seed);
	tOut->Branch("z3_seed",      &z3_seed);
	tOut->Branch("x4_seed",      &x4_seed);
	tOut->Branch("y4_seed",      &y4_seed);
	tOut->Branch("z4_seed",      &z4_seed);
   tOut->Branch("n_rec",     &n_rec);
	tOut->Branch("jgen_rec",  &jgen_rec);
	tOut->Branch("jseed_rec", &jseed_rec);
	tOut->Branch("q_rec",     &q_rec);
   tOut->Branch("p_rec",     &p_rec);
   tOut->Branch("acc_rec",   &acc_rec);
   tOut->Branch("polm_rec",  &polm_rec);
   tOut->Branch("poll_rec",  &poll_rec);
	tOut->Branch("issig_rec", &issig_rec);
	

   TH1D* h_dPrel_seed_gen = new TH1D("h_dPrel_seed_gen","Seed vs Gen;(P_{seed}-P_{gen})/P_{gen};Tracks",500,-0.05,+0.05);
   TH1D* h_dErel_seed_gen = new TH1D("h_dErel_seed_gen","Seed vs Gen;(E_{seed}-E_{gen})/E_{gen};Tracks",500,-0.05,+0.05);

   TH1D* h_dPrel_rec_gen = new TH1D("h_dPrel_rec_gen","Rec vs Gen;(P_{rec}-P_{gen})/P_{gen};Tracks",500,-0.05,+0.05);
   TH1D* h_dErel_rec_gen = new TH1D("h_dErel_rec_gen","Rec vs Gen;(E_{rec}-E_{gen})/E_{gen};Tracks",500,-0.05,+0.05);

   TH1D* h_dPrel_rec_seed = new TH1D("h_dPrel_rec_seed","Rec vs Seed;(P_{rec}-P_{seed})/P_{seed};Tracks",500,-0.05,+0.05);
   TH1D* h_dErel_rec_seed = new TH1D("h_dErel_rec_seed","Rec vs Seed;(E_{rec}-E_{seed})/E_{seed};Tracks",500,-0.05,+0.05);
 
   /// loop on events
   for(int iev=0;iev<nev;iev++)
   // for(int iev=0;iev<nev and iev<10;iev++)
   {
		/// global
		int n_res = 0;
		TLorentzVector ptmp;
		
      //// clear
	   n_gen = 0;
	   j_gen.clear();
		q_gen.clear();
	   p_gen.clear();
		
	   n_seed = 0;
		jgen_seed.clear();
		q_seed.clear();
	   p_seed.clear();
	   acc_seed.clear();
 	   for(int i=0;i<(int)polm_seed.size();++i) delete polm_seed[i];
 	   for(int i=0;i<(int)poll_seed.size();++i) delete poll_seed[i];
	   polm_seed.clear();
	   poll_seed.clear();
		svd0_seed.clear();
		svd1_seed.clear();
		svd2_seed.clear();
		chi2xz_seed.clear();
		chi2yz_seed.clear();
		residxz_seed.clear();
		residyz_seed.clear();
		issig_seed.clear();
		x1_seed.clear();
		y1_seed.clear();
		z1_seed.clear();
		x2_seed.clear();
		y2_seed.clear();
		z2_seed.clear();
		x3_seed.clear();
		y3_seed.clear();
		z3_seed.clear();
		x4_seed.clear();
		y4_seed.clear();
		z4_seed.clear();
		
	   n_rec = 0;
		jgen_rec.clear();
	   jseed_rec.clear();
		q_rec.clear();
	   p_rec.clear();
	   acc_rec.clear();
 	   for(int i=0;i<(int)polm_rec.size();++i) delete polm_rec[i];
 	   for(int i=0;i<(int)poll_rec.size();++i) delete poll_rec[i];
	   polm_rec.clear();
	   poll_rec.clear();
	   issig_rec.clear();
      
 	   //// get the next entry
 	   tIn->GetEntry(iev);
      if((iev%outN)==0) printf("Done %d out of %d\n",iev,nev);

      /// write all generated tracks
      for(unsigned int k=0 ; k<pxGen->size() ; k++)
      {
			n_gen++;
			ptmp.SetXYZM(pxGen->at(k), pyGen->at(k), pzGen->at(k), meGeV);
			p_gen.push_back(ptmp);
			q_gen.push_back(qGen->at(k));
			j_gen.push_back(iGen->at(k));
		}

      /// loop on all seeds
      for(unsigned int iseed=0 ; iseed<pxSeed->size() ; iseed++)
      {
			/// cut on some quality
         bool goodseed = (svd1Seed->at(iseed)<0.005 and svd2Seed->at(iseed)<0.0025);
			
         n_seed++;
			int pdgId = (x4Seed->at(iseed)>0) ? 11 : -11;
         int crg   = (x4Seed->at(iseed)>0) ? -1 : +1;
			ptmp.SetXYZM(pxSeed->at(iseed), pySeed->at(iseed), pzSeed->at(iseed), meGeV);
         p_seed.push_back(ptmp);
         q_seed.push_back(crg);
         acc_seed.push_back(0);
			int igenmatch_seed = iGenMatch->at(iseed);
         jgen_seed.push_back( igenmatch_seed );
			svd0_seed.push_back(svd0Seed->at(iseed));
			svd1_seed.push_back(svd1Seed->at(iseed));
			svd2_seed.push_back(svd2Seed->at(iseed));
			chi2xz_seed.push_back(chi2xzSeed->at(iseed));
			chi2yz_seed.push_back(chi2yzSeed->at(iseed));
			residxz_seed.push_back(residxzSeed->at(iseed));
			residyz_seed.push_back(residyzSeed->at(iseed));
			issig_seed.push_back(issigSeed->at(iseed));
			x1_seed.push_back(x1Seed->at(iseed));
			y1_seed.push_back(y1Seed->at(iseed));
			z1_seed.push_back(z1Seed->at(iseed));
			x2_seed.push_back(x2Seed->at(iseed));
			y2_seed.push_back(y2Seed->at(iseed));
			z2_seed.push_back(z2Seed->at(iseed));
			x3_seed.push_back(x3Seed->at(iseed));
			y3_seed.push_back(y3Seed->at(iseed));
			z3_seed.push_back(z3Seed->at(iseed));
			x4_seed.push_back(x4Seed->at(iseed));
			y4_seed.push_back(y4Seed->at(iseed));
			z4_seed.push_back(z4Seed->at(iseed));
			
         if(goodseed && igenmatch_seed>-1)
			{
            h_dErel_seed_gen->Fill((p_seed[iseed].E()-p_gen[igenmatch_seed].E())/p_gen[igenmatch_seed].E());
            h_dPrel_seed_gen->Fill((p_seed[iseed].P()-p_gen[igenmatch_seed].P())/p_gen[igenmatch_seed].P());
			}
			

         /// rest all the layers of the detector (including inactive if any):
			for(Int_t i=0 ; i<det->GetLayers()->GetEntries() ; i++) det->GetLayer(i)->Reset();
			
			/// set the clusters of the seed
         double clxyzTrk1[3],clxyzLab1[3], clxyzTrk2[3],clxyzLab2[3], clxyzTrk3[3],clxyzLab3[3], clxyzTrk4[3],clxyzLab4[3];
         clxyzLab1[0]=x1Seed->at(iseed); clxyzLab1[1]=y1Seed->at(iseed); clxyzLab1[2]=z1Seed->at(iseed);
         clxyzLab2[0]=x2Seed->at(iseed); clxyzLab2[1]=y2Seed->at(iseed); clxyzLab2[2]=z2Seed->at(iseed);
         clxyzLab3[0]=x3Seed->at(iseed); clxyzLab3[1]=y3Seed->at(iseed); clxyzLab3[2]=z3Seed->at(iseed);
         clxyzLab4[0]=x4Seed->at(iseed); clxyzLab4[1]=y4Seed->at(iseed); clxyzLab4[2]=z4Seed->at(iseed);
         KMCProbeFwd::Lab2Trk(clxyzLab1, clxyzTrk1);
         KMCProbeFwd::Lab2Trk(clxyzLab2, clxyzTrk2);
         KMCProbeFwd::Lab2Trk(clxyzLab3, clxyzTrk3);
         KMCProbeFwd::Lab2Trk(clxyzLab4, clxyzTrk4);
         det->GetLayer(1)->GetMCCluster()->Set( clxyzTrk1[0], clxyzTrk1[1], clxyzTrk1[2]);
         det->GetLayer(3)->GetMCCluster()->Set( clxyzTrk2[0], clxyzTrk2[1], clxyzTrk2[2]);
         det->GetLayer(5)->GetMCCluster()->Set( clxyzTrk3[0], clxyzTrk3[1], clxyzTrk3[2]);
         det->GetLayer(7)->GetMCCluster()->Set( clxyzTrk4[0], clxyzTrk4[1], clxyzTrk4[2]);
         
			// prepare the probe from the seed and do the KF fit
			bool res = det->SolveSingleTrackViaKalmanMC_Noam(p_seed[iseed].Pt(),p_seed[iseed].Rapidity(),p_seed[iseed].Phi(), meGeV, crg, vX,vY,vZ,99);
         if(!res) { cout << "SolveSingleTrackViaKalmanMC_Noam failed" << endl; continue; } // reconstruction failed
         n_res++;

         // get the truth trajectory in the B-field
         KMCProbeFwd* seedtrk = det->GetProbe();
         poll_seed.push_back( TrackLine3d(seedtrk,361,1,trkcol(p_seed[iseed].E())) );
         polm_seed.push_back( TrackMarker3d(seedtrk,0,361,1,trkcol(p_seed[iseed].E())) );


         // check truth acceptance
         acc_seed[iseed] = accepttrk(polm_seed,iseed);

         // get the reconstructed propagated to the vertex 
         KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack(); 
         if(!trw) continue; // track was not reconstructed
         n_rec++;
			
			int irec = n_rec-1;
         double pxyz[3];
         double xyz[3];
         trw->GetPXYZ(pxyz);
         trw->GetXYZ(xyz);
			
			jgen_rec.push_back( igenmatch_seed );
			jseed_rec.push_back( iseed );
         acc_rec.push_back(0);
         ptmp.SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
         p_rec.push_back(ptmp);
         q_rec.push_back(crg);
			issig_rec.push_back(issigSeed->at(iseed));

         if(goodseed)
			{
            h_dErel_rec_seed->Fill((p_rec[irec].E()-p_seed[iseed].E())/p_seed[iseed].E());
            h_dPrel_rec_seed->Fill((p_rec[irec].P()-p_seed[iseed].P())/p_seed[iseed].P());
            if(igenmatch_seed>-1)
			   {
               h_dErel_rec_gen->Fill((p_rec[irec].E()-p_gen[igenmatch_seed].E())/p_gen[igenmatch_seed].E());
               h_dPrel_rec_gen->Fill((p_rec[irec].P()-p_gen[igenmatch_seed].P())/p_gen[igenmatch_seed].P());
            }
			}
         			
         /// get the tracks
         poll_rec.push_back( TrackLine3d(trw,361,1,trkcol(p_rec[irec].E())) );
         polm_rec.push_back( TrackMarker3d(trw,0,361,1,trkcol(p_rec[irec].E())) );

         /// in acceptance?
         acc_rec[irec] = accepttrk(polm_rec,irec);
      }
      if(n_res!=n_seed) cout << "Warning: n_res=" << n_res << ", n_seed=" << n_seed << " --> problem" << endl;
      fOut->cd();
      tOut->Fill();
		printf("Event %d: Ngen=%d, Nseed=%d and Nrec=%d\n",iev,n_gen,n_seed,n_rec);
   }
   fOut->cd();
   tOut->Write();
   h_dErel_seed_gen->Write();
   h_dPrel_seed_gen->Write();
   h_dErel_rec_gen->Write();
   h_dPrel_rec_gen->Write();
   h_dErel_rec_seed->Write();
   h_dPrel_rec_seed->Write();
   fOut->Write();
   fOut->Close();
}
