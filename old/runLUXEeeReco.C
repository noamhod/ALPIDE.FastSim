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

double vX=0,vY=0,vZ=0; // event vertex
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
    double zCurr = xyz[2]; //source->GetZ();
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
	 // vector<double> xp,yp,zp;
	 // xp.push_back(xyz[0]);
	 // yp.push_back(xyz[1]);
	 // zp.push_back(xyz[2]);
	 // int iz = 0;
	 // while(zp[zp.size()-1]<zMax)
    for(int iz=1 ; iz<nZ ; iz++)
	 {
       if(!det->PropagateToZBxByBz(&tmp, TMath::Min(tmp.GetZ()+step, zMax), step)) break;
       tmp.GetXYZ(xyz);
       xp[iz] = xyz[0];
       yp[iz] = xyz[1];
       zp[iz] = xyz[2];
       nz++;
		 // xp.push_back(xyz[0]);
		 // yp.push_back(xyz[1]);
		 // zp.push_back(xyz[2]);
    }
    TPolyLine3D *polyline = new TPolyLine3D(nz+1);
    // TPolyLine3D *polyline = new TPolyLine3D();
	 polyline->SetLineColor(col);
    for(int i=0 ; i<nz+1 ; i++)
    // for(unsigned int i=0 ; i<xp.size() ; i++)
	 {
       polyline->SetPoint(i,xp[i],yp[i],zp[i]);
       // polyline->SetNextPoint(xp[i],yp[i],zp[i]);
    }
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
			if(dz<stepsize) return true;
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

TPolyLine3D* GetTruthTrack(double* xyz0, TLorentzVector p, double R, Color_t col)
{
   double pHat[] = { p.Px()/p.Vect().Mag(), p.Py()/p.Vect().Mag(), p.Pz()/p.Vect().Mag() };
   double xyz1[] = { pHat[0]*R, pHat[1]*R, pHat[2]*R };
   Int_t n=2;
   Double_t x[] = {xyz0[0],xyz1[0]};
   Double_t y[] = {xyz0[1],xyz1[1]};
   Double_t z[] = {xyz0[2],xyz1[2]};
   TPolyLine3D* polyline = new TPolyLine3D(n,x,y,z);
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
   
	cout << "poll.size()=" << poll.size() << ", clusters.size()=" << clusters.size() <<endl;
   for(int i=0 ; i<(int)poll.size() ; ++i)
	{
		cnv_pl3d->cd();
		if(inacc[i])
		{
			poll[i]->Draw();
			clusters[i]->Draw();
		}
	}
   for(int i=0 ; i<(int)polm.size() ; ++i) { cnv_pm3d->cd(); if(inacc[i]) polm[i]->Draw(); }
   
   TLegend* leg = trkcolleg();
   cnv_pl3d->cd();
   leg->Draw("same");
   cnv_pm3d->cd();
   leg->Draw("same");
   
   cnv_pl3d->SaveAs("output/root/"+process+"_tracks_pl3d"+suff+".root");
   cnv_pl3d->SaveAs("output/pdf/"+process+"_tracks_pl3d"+suff+".pdf");
   cnv_pm3d->SaveAs("output/root/"+process+"_tracks_pm3d"+suff+".root");
   cnv_pm3d->SaveAs("output/pdf/"+process+"_tracks_pm3d"+suff+".pdf");
   
   TFile* flines = new TFile("data/root/"+process+"_geometry"+suff+".root","RECREATE");
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

bool accepttrk4(TPolyMarker3D* marker, double zstep=1)
{
	/// in acceptance of layer 4?
   for(int i=marker->GetN()-1 ; i>=0 ; --i)
   {
      Double_t xr,yr,zr;
      marker->GetPoint(i,xr,yr,zr);
      if(islayer(zr,zlayer->size()-1,zstep))
		{
			bool inacc = accept(xr,yr);
			// cout << "r["<<i<<"]={"<<xr<<","<<yr<<","<<zr<<"} --> inacc? " << inacc << endl;
			return inacc;
		}
   }
   return false;
}

/// direct track parameters
vector<double> trk_rec_Tgl;
vector<double> trk_rec_Snp; // the slope in X direction: probe->GetTrack()->GetSnp()
vector<double> trk_rec_alpha;
vector<double> trk_rec_signedinvpT; // new: the curvature (q/Pyz): probe->GetTrack()->GetSigned1Pt()
vector<double> trk_rec_sigmaY2;
vector<double> trk_rec_sigmaZY;
vector<double> trk_rec_sigmaZ2;
vector<double> trk_rec_sigmaSnpY;
vector<double> trk_rec_sigmaSnpZ;
vector<double> trk_rec_sigmaSnp2; // probe->GetTrack()->GetSigmaSnp2()
vector<double> trk_rec_sigmaTglY;
vector<double> trk_rec_sigmaTglZ;
vector<double> trk_rec_sigmaTglSnp;
vector<double> trk_rec_sigmaTgl2;
vector<double> trk_rec_sigma1PtY;
vector<double> trk_rec_sigma1PtZ;
vector<double> trk_rec_sigma1PtSnp;
vector<double> trk_rec_sigma1PtTgl;
vector<double> trk_rec_sigma1Pt2;
vector<double> trk_rec_invpT;
vector<double> trk_rec_signedpT;

void settrkvecbranches(TTree* t)
{
   t->Branch("trk_rec_Tgl",         &trk_rec_Tgl        );
   t->Branch("trk_rec_Snp",         &trk_rec_Snp        );
   t->Branch("trk_rec_alpha",       &trk_rec_alpha      );
   t->Branch("trk_rec_signedinvpT", &trk_rec_signedinvpT);
   t->Branch("trk_rec_sigmaY2",     &trk_rec_sigmaY2    );
   t->Branch("trk_rec_sigmaZY",     &trk_rec_sigmaZY    );
   t->Branch("trk_rec_sigmaZ2",     &trk_rec_sigmaZ2    );
   t->Branch("trk_rec_sigmaSnpY",   &trk_rec_sigmaSnpY  );
   t->Branch("trk_rec_sigmaSnpZ",   &trk_rec_sigmaSnpZ  );
   t->Branch("trk_rec_sigmaSnp2",   &trk_rec_sigmaSnp2  );
   t->Branch("trk_rec_sigmaTglY",   &trk_rec_sigmaTglY  );
   t->Branch("trk_rec_sigmaTglZ",   &trk_rec_sigmaTglZ  );
   t->Branch("trk_rec_sigmaTglSnp", &trk_rec_sigmaTglSnp);
   t->Branch("trk_rec_sigmaTgl2",   &trk_rec_sigmaTgl2  );
   t->Branch("trk_rec_sigma1PtY",   &trk_rec_sigma1PtY  );
   t->Branch("trk_rec_sigma1PtZ",   &trk_rec_sigma1PtZ  );
   t->Branch("trk_rec_sigma1PtSnp", &trk_rec_sigma1PtSnp);
   t->Branch("trk_rec_sigma1PtTgl", &trk_rec_sigma1PtTgl);
   t->Branch("trk_rec_sigma1Pt2",   &trk_rec_sigma1Pt2  );
   t->Branch("trk_rec_invpT",       &trk_rec_invpT      );
   t->Branch("trk_rec_signedpT",    &trk_rec_signedpT   );
}

void cleartrkvec()
{
   trk_rec_Tgl.clear();          
   trk_rec_Snp.clear();
   trk_rec_alpha.clear();
   trk_rec_signedinvpT.clear();   
   trk_rec_sigmaY2.clear();
   trk_rec_sigmaZY.clear();
   trk_rec_sigmaZ2.clear();
   trk_rec_sigmaSnpY.clear();
   trk_rec_sigmaSnpZ.clear();
   trk_rec_sigmaSnp2.clear();
   trk_rec_sigmaTglY.clear();
   trk_rec_sigmaTglZ.clear();
   trk_rec_sigmaTglSnp.clear();
   trk_rec_sigmaTgl2.clear();
   trk_rec_sigma1PtY.clear();
   trk_rec_sigma1PtZ.clear();
   trk_rec_sigma1PtSnp.clear();
   trk_rec_sigma1PtTgl.clear();
   trk_rec_sigma1Pt2.clear();
   trk_rec_invpT.clear();
   trk_rec_signedpT.clear();
}

void filltrkvec(KMCProbeFwd* probe)
{
	/// TODO: this is not good --> the PxPyPz at least are given in the track frame and not in the lab.
	/// TODO: instead, need to do something like: double pxyztmp[3]; probe->GetPXYZ(pxyztmp); and take the values from there
	TrackPar* trk = probe->GetTrack();
   trk_rec_Tgl.push_back( trk->GetTgl() );         
   trk_rec_Snp.push_back( trk->GetSnp() );         
   trk_rec_alpha.push_back( trk->GetAlpha() );       
   trk_rec_signedinvpT.push_back( trk->GetSigned1Pt() ); 
   trk_rec_sigmaY2.push_back( trk->GetSigmaY2() );     
   trk_rec_sigmaZY.push_back( trk->GetSigmaZY() );     
   trk_rec_sigmaZ2.push_back( trk->GetSigmaZ2() );     
   trk_rec_sigmaSnpY.push_back( trk->GetSigmaSnpY() );   
   trk_rec_sigmaSnpZ.push_back( trk->GetSigmaSnpZ() );   
   trk_rec_sigmaSnp2.push_back( trk->GetSigmaSnp2() );   
   trk_rec_sigmaTglY.push_back( trk->GetSigmaTglY() );   
   trk_rec_sigmaTglZ.push_back( trk->GetSigmaTglZ() );   
   trk_rec_sigmaTglSnp.push_back( trk->GetSigmaTglSnp() ); 
   trk_rec_sigmaTgl2.push_back( trk->GetSigmaTgl2() );   
   trk_rec_sigma1PtY.push_back( trk->GetSigma1PtY() );   
   trk_rec_sigma1PtZ.push_back( trk->GetSigma1PtZ() );   
   trk_rec_sigma1PtSnp.push_back( trk->GetSigma1PtSnp() ); 
   trk_rec_sigma1PtTgl.push_back( trk->GetSigma1PtTgl() ); 
   trk_rec_sigma1Pt2.push_back( trk->GetSigma1Pt2() );   
   trk_rec_invpT.push_back( trk->OneOverPt() );                
   trk_rec_signedpT.push_back( trk->GetSignedPt() );    
}


void runLUXEeeReco(TString process, int Seed=12345) //, const char* setup="setup/setupLUXE.txt")
{
	TString proc = process;
	proc.ReplaceAll("_bkg","");
	TString setup = "setup/setupLUXE_"+proc+".txt";
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
   TFile* fIn = new TFile("data/root/raw_"+process+".root","READ");
   TTree* tIn = (TTree*)fIn->Get("tt");
   int nev = tIn->GetEntries();
   vector<double>* vx    = 0;
   vector<double>* vy    = 0;
   vector<double>* vz    = 0;
   vector<double>* px    = 0;
   vector<double>* py    = 0;
   vector<double>* pz    = 0;
   vector<double>* E     = 0;
   vector<double>* wgt   = 0;
   vector<int>*    pdgId = 0;
   tIn->SetBranchAddress("vx",&vx);
   tIn->SetBranchAddress("vy",&vy);
   tIn->SetBranchAddress("vz",&vz);
   tIn->SetBranchAddress("px",&px);
   tIn->SetBranchAddress("py",&py);
   tIn->SetBranchAddress("pz",&pz);
   tIn->SetBranchAddress("E",&E);
   tIn->SetBranchAddress("wgt",&wgt);
   tIn->SetBranchAddress("pdgId",&pdgId);
 
   // output tree
   gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",   "vector");
   int ngen = 0;
   int nres = 0;
   // int nrec = 0;
	int nacc_gen = 0;
	// int nacc_rec = 0;
   vector<double>          wgtgen;
   vector<float>           xgen;
   vector<float>           ygen;
   vector<float>           zgen;
   vector<TLorentzVector>  pgen;
   // vector<TLorentzVector>  prec;
   vector<int>             jpair;
   // vector<int>             jrec;
   vector<int>             qgen;
   // vector<int>             qrec;
   vector<int>             igentrk;
	vector<vector<int> >    polm_clusters_id;
	vector<TPolyMarker3D*>  polm_clusters;
   // vector<TPolyMarker3D*>  polm;
   // vector<TPolyLine3D*>    poll;
   vector<TPolyMarker3D*>  polm_gen;
   vector<TPolyLine3D*>    poll_gen;
   vector<int>             acctrkgen;
   // vector<int>             acctrkrec;
   // vector<float>           vtxtrkrec;
   /// direct detector parameters
   // vector<int>    nITSHits;

   TFile* fOut = new TFile("data/root/rec_"+process+".root","RECREATE");
   TTree* tOut = new TTree("res","res");
   tOut->Branch("wgtgen",&wgtgen);
   tOut->Branch("xgen",&xgen);
   tOut->Branch("ygen",&ygen);
   tOut->Branch("zgen",&zgen);
   tOut->Branch("jpair",&jpair);
   // tOut->Branch("jrec",&jrec);
   tOut->Branch("ngen",&ngen);
   // tOut->Branch("nrec",&nrec);
   tOut->Branch("qgen",&qgen);
   tOut->Branch("pgen",&pgen);
   // tOut->Branch("prec",&prec);
   tOut->Branch("igentrk",&igentrk);
   tOut->Branch("acctrkgen",&acctrkgen);
   // tOut->Branch("acctrkrec",&acctrkrec);
   // tOut->Branch("vtxtrkrec",&vtxtrkrec);
   tOut->Branch("polm_clusters_id",&polm_clusters_id);
   tOut->Branch("polm_clusters",&polm_clusters);
   // tOut->Branch("polm",&polm);
   tOut->Branch("polm_gen",&polm_gen);
   tOut->Branch("poll_gen",&poll_gen);
   // tOut->Branch("poll",&poll);
   // tOut->Branch("nITSHits",&nITSHits);
   // settrkvecbranches(tOut);

   // TH1D* h_dPrel_Ruben = new TH1D("h_dPrel_Ruben",";(P_{rec}-P_{gen})/P_{gen};Tracks",500,-0.01,+0.01);
   // TH1D* h_dPrel = new TH1D("h_dPrel",";(P_{rec}-P_{gen})/P_{gen};Tracks",500,-0.01,+0.01);
   // TH1D* h_dErel = new TH1D("h_dErel",";(E_{rec}-E_{gen})/E_{gen};Tracks",500,-0.01,+0.01);
   // TH1D* h_dxrel = new TH1D("h_dxrel",";(x_{rec}-x_{gen})/x_{gen};Tracks",200,-0.01,+0.01);
 
   /// loop on events
   for(int iev=0;iev<nev;iev++)
   // for(int iev=0;iev<nev and iev<50;iev++)
   {	
      //// clear
      ngen = 0;    
      nres = 0;
      // nrec = 0;
      nacc_gen = 0;
      // nacc_rec = 0;
		
 	   for(int i=0;i<(int)polm_clusters_id.size();++i) polm_clusters_id[i].clear();
 	   for(int i=0;i<(int)polm_clusters.size();++i) delete polm_clusters[i];
 	   for(int i=0;i<(int)polm_gen.size();++i) delete polm_gen[i];
 	   for(int i=0;i<(int)poll_gen.size();++i) delete poll_gen[i];
 	   // for(int i=0;i<(int)polm.size();++i)     delete polm[i];
      // for(int i=0;i<(int)poll.size();++i)     delete poll[i];
      wgtgen.clear();
      xgen.clear();
      ygen.clear();
      zgen.clear();
      pgen.clear();
      // prec.clear();
      qgen.clear();
      // qrec.clear();
      igentrk.clear();
		polm_clusters_id.clear();
      polm_clusters.clear();
      polm_gen.clear();
      poll_gen.clear();
      // polm.clear();
      // poll.clear();
      acctrkgen.clear();
      // acctrkrec.clear();
      // vtxtrkrec.clear();
      jpair.clear();
      // jrec.clear();
      // nITSHits.clear();
      // cleartrkvec();

      int npairs = 0;
      int npositrons = 0;
      
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
         int crg = (pdgId->at(igen)==11) ? -1 : +1;
         wgtgen.push_back(wgt->at(igen));
         xgen.push_back(vx->at(igen));
         ygen.push_back(vy->at(igen));
         zgen.push_back(vz->at(igen));
         pgen.push_back(ptmp);
         qgen.push_back(crg);
         acctrkgen.push_back(0);
			polm_clusters_id.push_back( vtmp );
			polm_clusters.push_back( new TPolyMarker3D() );
         // jrec.push_back(-999);
         pgen[igen].SetXYZM(px->at(igen), py->at(igen), pz->at(igen), meGeV);
			// cout << "pgen["<<igen<<"].Pz()="<<pgen[igen].Pz()<<", m="<<pgen[igen].M()<<endl;
         if(pdgId->at(igen)==-11) npositrons ++; // count positrons
         /// truth e+e- indices pairing for later analysis
         if(igen<(ngenall-1))
         {
            bool isel = (pdgId->at(igen)==11); // avoid double counting
            bool is0q = (pdgId->at(igen)+pdgId->at(igen+1)==0); // oposite sign
            bool is1w = (((int)wgt->at(igen)==1) and ((int)wgt->at(igen+1)==1));
            if(isel and is0q and is1w)
            {
               npairs++; // count pairs
               jpair.push_back(igen+1);
               jpair.push_back(igen);
            }
         }
         if((int)jpair.size()<(igen+1)) jpair.push_back(-999);
			poll_gen.push_back( 0 );
			polm_gen.push_back( 0 );
			
         
         // prepare the probe
			// cout << "pgen["<<igen<<"].Pt()="<<pgen[igen].Pt()<<", pgen["<<igen<<"].Rapidity()="<<pgen[igen].Rapidity()<<", pgen["<<igen<<"].Phi()="<<pgen[igen].Phi()<<", crg="<<crg<<", vtx={"<<vx->at(igen)<<","<<vy->at(igen)<<","<<vz->at(igen)<<"}"<<endl;
         bool res = det->SolveSingleTrack(pgen[igen].Pt(),pgen[igen].Rapidity(),pgen[igen].Phi(), meGeV, crg, vx->at(igen),vy->at(igen),vz->at(igen), 0,1,99);
         // bool res = det->SolveSingleTrack(pgen[igen].Pt(),pgen[igen].Rapidity(),pgen[igen].Phi(), meGeV, crg, vX,vY,vZ, 0,1,99);
         if(!res) continue; // reconstruction failed
         nres++;
			
			
         // get the truth trajectory in the B-field
         KMCProbeFwd* trutrk = det->GetProbe();
			double pxyztmp[3];
			trutrk->GetPXYZ(pxyztmp);
			// cout << "igen="<<igen<<" --> pxyz[2]=" << pxyztmp[2] << endl;
			if(pxyztmp[2]<0) continue;
			
			
			poll_gen[igen] = TrackLine3d(trutrk,361,1,trkcol(pgen[igen].E()));
			polm_gen[igen] = TrackMarker3d(trutrk,0,361,1,trkcol(pgen[igen].E()));

			/// check acceptance
         acctrkgen[igen] = accepttrk4(polm_gen[igen]);
			if(acctrkgen[igen]) nacc_gen++;
			
			// get the reconstructed propagated to the vertex
         KMCClusterFwd* cluster1 = det->GetLayer(1)->GetMCCluster();
         KMCClusterFwd* cluster2 = det->GetLayer(3)->GetMCCluster();
         KMCClusterFwd* cluster3 = det->GetLayer(5)->GetMCCluster();
         KMCClusterFwd* cluster4 = det->GetLayer(7)->GetMCCluster();
			
	  	   polm_clusters[igen]->SetNextPoint(cluster1->GetXLab(),cluster1->GetYLab(),cluster1->GetZLab());
	  	   polm_clusters[igen]->SetNextPoint(cluster2->GetXLab(),cluster2->GetYLab(),cluster2->GetZLab());
	  	   polm_clusters[igen]->SetNextPoint(cluster3->GetXLab(),cluster3->GetYLab(),cluster3->GetZLab());
	  	   polm_clusters[igen]->SetNextPoint(cluster4->GetXLab(),cluster4->GetYLab(),cluster4->GetZLab());

			polm_clusters_id[igen].push_back( 1*index_offset+igen ); // assuming no chance to have >index_offset tracks
			polm_clusters_id[igen].push_back( 2*index_offset+igen ); // assuming no chance to have >index_offset tracks
			polm_clusters_id[igen].push_back( 3*index_offset+igen ); // assuming no chance to have >index_offset tracks
			polm_clusters_id[igen].push_back( 4*index_offset+igen ); // assuming no chance to have >index_offset tracks
			
			//          // get the reconstructed propagated to the vertex
			//          KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack();
			//          if(!trw) continue; // track was not reconstructed
			//
			// /// TODO: need to understand why this is happenning
			// /// TODO: i.e. for which truth particles the vertex cannot be at {x,y}={0,0}
			//          double pxyz[3];
			//          double xyz[3];
			//          trw->GetPXYZ(pxyz);
			//          trw->GetXYZ(xyz);
			// // cout << "xyz={"<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<"}" << endl;
			// // if(abs(xyz[0])>10 or abs(xyz[1])>10) continue; // force vertex.
			//
			//          nrec++;
			//
			//
			//          int irec = nrec-1;
			//          jrec[igen] = irec;
			//          acctrkrec.push_back(0);
			//          vtxtrkrec.push_back(-9999);
			//          prec.push_back(ptmp);
			//          qrec.push_back(crg);
			//          igentrk.push_back(igen);
			//
			//
			//
			//          prec[irec].SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);
			//          nfakeHits += trw->GetNFakeITSHits(); // in absence of bkg, there should be no fake hits, but count anyway
			//          nITSHits.push_back(trw->GetNITSHits());
			//          filltrkvec(trw); // fill rec trk parameters
			//
			//          h_dErel->Fill((prec[irec].E()-pgen[igen].E())/pgen[igen].E());
			//          h_dPrel->Fill((prec[irec].P()-pgen[igen].P())/pgen[igen].P());
			//          h_dPrel_Ruben->Fill(  (trw->GetP()-trutrk->GetP())/trutrk->GetP() );
			//
			//          /// get the tracks
			//          poll.push_back( TrackLine3d(trw,361,1,trkcol(prec[irec].E())) );
			//          polm.push_back( TrackMarker3d(trw,0,361,1,trkcol(prec[irec].E())) );
			//
			// /// check acceptance
			// vtxtrkrec[irec] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
			// acctrkrec[irec] = accepttrk4(polm[irec]);
			// if(acctrkrec[irec]) nacc_rec++;
      }
		
		// if(iev%1==0) cout << "iev=" << iev << " --> nacc_gen=" << nacc_gen << ", nacc_rec=" << nacc_rec << endl;
		if(iev%1==0) cout << "iev=" << iev << " --> nacc_gen=" << nacc_gen << endl;
		
      if(nres!=ngen and !process.Contains("bkg")) cout << "Warning: nres=" << nres << ", ngen=" << ngen << " --> problem" << endl;
      // if(iev==0) WriteGeometry(polm,poll,process,acctrkrec,polm_clusters);
      if(iev==0) WriteGeometry(polm_gen,poll_gen,process,acctrkgen,polm_clusters,"_truth");
		
      fOut->cd();
      tOut->Fill();
   }
   fOut->cd();
   tOut->Write();
   // h_dxrel->Write();
   // h_dErel->Write();
   // h_dPrel->Write();
   // h_dPrel_Ruben->Write();
   fOut->Write();
   fOut->Close();
}
