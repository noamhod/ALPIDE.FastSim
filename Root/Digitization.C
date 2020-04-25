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
#endif

KMCDetectorFwd* det = 0;
vector<double>* zlayer = new vector<double>;
double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 27;   // cm for BPPP or 50 for Trident
double Rbeampipe = 4; // cm for BPPP or 14 for Trident
double RoffsetBfield22BPPP = 7.0; // cm for BPPP in B=2.2T
double RoffsetBfield20BPPP = 5.7; // cm for BPPP in B=2.0T
double RoffsetBfield14BPPP = 4.0; // cm for BPPP in B=1.4T
double RoffsetBfield = RoffsetBfield20BPPP;
double xPsideL = -RoffsetBfield-Lstave;
double xPsideR = -RoffsetBfield;       
double xEsideL = +RoffsetBfield;       
double xEsideR = +RoffsetBfield+Lstave;
double yUp = +Hstave/2.;
double yDn = -Hstave/2.;

//// dipole geometry
double xW = 120;
double yH = 67.2;
double z1 = 100;
double z2 = 200;

void resetToTridentGeometry()
{
	Lstave = 50;   // cm for BPPP or 50 for Trident
	RoffsetBfield = 14; // cm for Trident in in B=1.4T
	xPsideL = -RoffsetBfield-Lstave;
	xPsideR = -RoffsetBfield;       
	xEsideL = +RoffsetBfield;       
	xEsideR = +RoffsetBfield+Lstave;
}

int acceptcls(double x, double y, double z)
{
   // bool failx = (x<xPsideL || (x>xPsideR && x<xEsideL) || x>xEsideR);
	bool failx = (abs(x)>xEsideR or abs(x)<xEsideL);
	if(failx) return 0;
   // bool faily = (y>yUp || y<yDn);
   bool faily = (abs(y)>yUp);
	if(faily) return 0;
	bool failz = (z!=300 && z!=310 && z!=320 && z!=330);
   if(failz) return 0;
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
   TLegend* leg = new TLegend(0.12,0.60,0.50,0.80);
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
    zp[0] = round(xyz[2]); // z-vertex
	 int nz = 0;
    for(int iz=1 ; iz<nZ ; iz++)
	 {
       if(!det->PropagateToZBxByBz(&tmp, TMath::Min(round(tmp.GetZ())+step, zMax), step)) break;
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
		double dz = abs(zlayer->at(layerindex)-z);
		if(dz<stepsize) return true;
	}
	else
	{
		for(int j=0 ; j<(int)zlayer->size() ; ++j)
		{
			double dz = abs(zlayer->at(j)-z);
			if(dz<stepsize/2.)
			{
				return true;
			}
		}
	}
	return false;
}

TPolyMarker3D* TrackMarker3d(const KMCProbeFwd* source, double zmin, double zmax, double zstep, Color_t col, bool doPrint=false)
{
    KMCProbeFwd tmp(*source);
	 int nZ = (int)(zmax-zmin)/zstep;
    double xp[nZ],yp[nZ],zp[nZ];
    double xyz[3];
    tmp.GetXYZ(xyz);
    xp[0] = xyz[0];
    yp[0] = xyz[1];
    zp[0] = round(xyz[2]);
	 int nz = 0;
    for(int iz=1;iz<nZ;iz++)
	 {
		 if(!det->PropagateToZBxByBz(&tmp, round(tmp.GetZ())+zstep, zstep)) break;
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
       if(!islayer(zp[i],-1,zstep)) continue;
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

TPolyLine3D* GetLayer(TString side, double z, Color_t col)
{
   Int_t n=5;
   Double_t xL[] = {xPsideL,xPsideL,xPsideR,xPsideR,xPsideL};
   Double_t xR[] = {xEsideR,xEsideR,xEsideL,xEsideL,xEsideR};
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

void WriteGeometry(vector<TPolyMarker3D*>& polm, vector<TPolyLine3D*>& poll, TString process, vector<int>& inacc, vector<TPolyMarker3D*>& clusters, TString suff="")
{
   TCanvas* cnv_pl3d = new TCanvas("cnv_pl3d"+suff,"",500,500);
   TView* view_pl3d = TView::CreateView(1);
   view_pl3d->SetRange(-80,-50,0, +80,+50,350);
   view_pl3d->ShowAxis();
   
   TCanvas* cnv_pm3d = new TCanvas("cnv_pm3d"+suff,"",500,500);
   TView* view_pm3d = TView::CreateView(1);
   view_pm3d->SetRange(-80,-50,0, +80,+50,350);
   view_pm3d->ShowAxis();
   
   TPolyLine3D* stave1L = GetLayer("L",300,kGreen+3);
   TPolyLine3D* stave1R = GetLayer("R",300,kGreen+3);
   TPolyLine3D* stave2L = GetLayer("L",310,kGreen+3);
   TPolyLine3D* stave2R = GetLayer("R",310,kGreen+3);
   TPolyLine3D* stave3L = GetLayer("L",320,kGreen+3);
   TPolyLine3D* stave3R = GetLayer("R",320,kGreen+3);
   TPolyLine3D* stave4L = GetLayer("L",330,kGreen+3);
   TPolyLine3D* stave4R = GetLayer("R",330,kGreen+3);
   TPolyLine3D* dipole  = GeDipole(kGray);
   
   cnv_pl3d->cd();
   dipole->Draw();
   stave1L->Draw();
   stave1R->Draw();
   stave2L->Draw();
   stave2R->Draw();
   stave3L->Draw();
   stave3R->Draw();
   stave4L->Draw();
   stave4R->Draw();
   
   cnv_pm3d->cd();
   dipole->Draw();
   stave1L->Draw();
   stave1R->Draw();
   stave2L->Draw();
   stave2R->Draw();
   stave3L->Draw();
   stave3R->Draw();
   stave4L->Draw();
   stave4R->Draw();
   
	cnv_pl3d->cd();
	vector<int> problems;
   for(int i=0 ; i<(int)poll.size() ; ++i)
	{
		if(!inacc[i]) continue;
		poll[i]->Draw();
		clusters[i]->Draw();
	}
	
	cnv_pm3d->cd();
   for(int i=0 ; i<(int)polm.size() ; ++i)
	{
		if(!inacc[i]) continue;
		polm[i]->Draw();
	}
   
   TLegend* leg = trkcolleg();
   cnv_pl3d->cd();
   leg->Draw("same");
   cnv_pm3d->cd();
   leg->Draw("same");
   
   cnv_pl3d->SaveAs("../output/root/"+process+"_tracks_pl3d"+suff+".root");
   cnv_pl3d->SaveAs("../output/pdf/"+process+"_tracks_pl3d"+suff+".pdf");
   cnv_pm3d->SaveAs("../output/root/"+process+"_tracks_pm3d"+suff+".root");
   cnv_pm3d->SaveAs("../output/pdf/"+process+"_tracks_pm3d"+suff+".pdf");
   
   TFile* flines = new TFile("../data/root/"+process+"_geometry"+suff+".root","RECREATE");
   flines->cd();
   dipole->Write();
   stave1L->Write();
   stave1R->Write();
   stave2L->Write();
   stave2R->Write();
   stave3L->Write();
   stave3R->Write();
   stave4L->Write();
   stave4R->Write();
   leg->Write();
   // flines->Write();
   flines->Close();
}

bool accepttrk(TPolyMarker3D* clusters, bool fullacc)
{
	/// in acceptance?
   int nlayers = 4;
   int acc = 0;
   Double_t xr,yr,zr;
   clusters->GetPoint(0,xr,yr,zr);
   acc += acceptcls(xr,yr,zr);
   clusters->GetPoint(1,xr,yr,zr);
   acc += acceptcls(xr,yr,zr);
   clusters->GetPoint(2,xr,yr,zr);
   acc += acceptcls(xr,yr,zr);
   clusters->GetPoint(3,xr,yr,zr);
   acc += acceptcls(xr,yr,zr);
   return (fullacc) ? (acc==nlayers) : (acc>0);
}

// bool accepttrk(TPolyMarker3D* marker, bool fullacc, double inflate=1, double zstep=1)
// {
// 	/// in acceptance?
//    int nlayers = 4;
//    int acc = 0;
//    for(int i=0 ; i<marker->GetN() ; i++)
//    {
//       Double_t xr,yr,zr;
//       marker->GetPoint(i,xr,yr,zr);
//       // if(zr<300) continue; //// count only the active layers
// 		if(!islayer(zr,-1,zstep)) continue;
//       int inacclayer = acceptcls(xr,yr,inflate,inflate);
//       acc += inacclayer;
//    }
//    return (fullacc) ? (acc==nlayers) : (acc>0);
// }

// bool accepttrk4(TPolyMarker3D* marker, double zstep=1)
// {
// 	/// in acceptance of layer 4?
//    for(int i=marker->GetN()-1 ; i>=0 ; --i)
//    {
//       Double_t xr,yr,zr;
//       marker->GetPoint(i,xr,yr,zr);
//       if(islayer(zr,zlayer->size()-1,zstep))
// 		{
// 			bool inacc = acceptcls(xr,yr,1.,1.);
// 			return inacc;
// 		}
//    }
//    return false;
// }

void Digitization(TString process, int Seed=12345) //, const char* setup="setup/setupLUXE.txt")
{
	TString proc = process;
	proc.ReplaceAll("_bkg","");
	TString setup = "../setup/setupLUXE_"+proc+".txt";
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
   det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
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

   int outN = 100;

   // TString process = "bppp";  /// trident or bppp or bppp_bkg or trident_bkg
	if(process=="trident") resetToTridentGeometry();
	
	int index_offset = (process.Contains("bkg")) ? 10000 : 100000;

   /// get the particles from a ttree
   TFile* fIn = new TFile("../data/root/raw_"+process+".root","READ");
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
   gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",   "vector");
   int ngen = 0;
   int nslv = 0;
	int nacc = 0;
   vector<double>          wgt;
   vector<float>           xvtx;
   vector<float>           yvtx;
   vector<float>           zvtx;
   vector<TLorentzVector>  trkp4;
   vector<int>             crg;
	vector<vector<int> >    clusters_id;
	vector<vector<int> >    clusters_type;
	vector<TPolyMarker3D*>  clusters_xyz;
   vector<TPolyMarker3D*>  trkpts;
   vector<TPolyLine3D*>    trklin;
   vector<int>             acc;

   TFile* fOut = new TFile("../data/root/dig_"+process+".root","RECREATE");
   TTree* tOut = new TTree("dig","dig");
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
   tOut->Branch("clusters_type",&clusters_type);
   tOut->Branch("clusters_xyz", &clusters_xyz);
   tOut->Branch("trkpts",       &trkpts);
   tOut->Branch("trklin",       &trklin);
	 
   /// loop on events
   for(int iev=0;iev<nev;iev++)
   // for(int iev=0;iev<nev and iev<50;iev++)
   {	
      //// clear
      ngen = 0;    
      nslv = 0;
      nacc = 0;
		
 	   for(int i=0;i<(int)clusters_id.size();++i) clusters_id[i].clear();
 	   for(int i=0;i<(int)clusters_type.size();++i) clusters_type[i].clear();
 	   for(int i=0;i<(int)clusters_xyz.size();++i) delete clusters_xyz[i];
 	   for(int i=0;i<(int)trkpts.size();++i) delete trkpts[i];
 	   for(int i=0;i<(int)trklin.size();++i) delete trklin[i];
      wgt.clear();
      xvtx.clear();
      yvtx.clear();
      zvtx.clear();
      trkp4.clear();
      crg.clear();
		clusters_id.clear();
		clusters_type.clear();
      clusters_xyz.clear();
      trkpts.clear();
      trklin.clear();
      acc.clear();
      
 	   //// get the next entry
 	   tIn->GetEntry(iev);
      if((iev%outN)==0) printf("Done %d out of %d\n",iev,nev);
      int nfakeHits = 0;
 	   TLorentzVector ptmp;
 	   int ngenall = pdgId->size();
      int pair = 1;
      /// loop on particles
      for(int igen=0 ; igen<ngenall ; igen++)
      {
         if(abs(pdgId->at(igen))!=11)
         {
			   cout << "illegal pdgId: " << pdgId->at(igen) << endl;
            break;
         }
			
         //// all the rest
         ngen++;
			vector<int> vtmp;
         int q = (pdgId->at(igen)==11) ? -1 : +1;
			ptmp.SetXYZM(px->at(igen), py->at(igen), pz->at(igen), meGeV);
			
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
			trklin.push_back( TrackLine3d(trutrk,361,1,trkcol(ptmp.E())) );
			trkpts.push_back( TrackMarker3d(trutrk,0,361,1,trkcol(ptmp.E())) );
         
			clusters_id.push_back( vtmp );
			clusters_type.push_back( vtmp );
			clusters_xyz.push_back( new TPolyMarker3D() );
			acc.push_back( 0 );
			
			// get the reconstructed propagated to the vertex
         KMCClusterFwd* cluster1 = det->GetLayer(1)->GetMCCluster();
         KMCClusterFwd* cluster2 = det->GetLayer(3)->GetMCCluster();
         KMCClusterFwd* cluster3 = det->GetLayer(5)->GetMCCluster();
         KMCClusterFwd* cluster4 = det->GetLayer(7)->GetMCCluster();
	  	   clusters_xyz[slvidx]->SetNextPoint(cluster1->GetXLab(),cluster1->GetYLab(),cluster1->GetZLab());
	  	   clusters_xyz[slvidx]->SetNextPoint(cluster2->GetXLab(),cluster2->GetYLab(),cluster2->GetZLab());
	  	   clusters_xyz[slvidx]->SetNextPoint(cluster3->GetXLab(),cluster3->GetYLab(),cluster3->GetZLab());
	  	   clusters_xyz[slvidx]->SetNextPoint(cluster4->GetXLab(),cluster4->GetYLab(),cluster4->GetZLab());
			clusters_id[slvidx].push_back( 1*index_offset+slvidx ); // assuming no chance to have >index_offset tracks
			clusters_id[slvidx].push_back( 2*index_offset+slvidx ); // assuming no chance to have >index_offset tracks
			clusters_id[slvidx].push_back( 3*index_offset+slvidx ); // assuming no chance to have >index_offset tracks
			clusters_id[slvidx].push_back( 4*index_offset+slvidx ); // assuming no chance to have >index_offset tracks
			clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
			clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
			clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
			clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
			
			/// check acceptance
			acc[slvidx] = accepttrk(clusters_xyz[slvidx],false);
			if(acc[slvidx]) nacc++;
			
      }
		if(iev==0) WriteGeometry(trkpts,trklin,process,acc,clusters_xyz,"_truth");
		if(iev%1==0) cout << "iev=" << iev << " --> ngen=" << ngen << ", nslv=" << nslv << ", nacc=" << nacc << endl;
      if(nslv!=ngen and !process.Contains("bkg")) cout << "Warning: nslv=" << nslv << ", ngen=" << ngen << " --> problem" << endl;
		
      fOut->cd();
      tOut->Fill();
   }
   fOut->cd();
   tOut->Write();
   fOut->Write();
   fOut->Close();
}