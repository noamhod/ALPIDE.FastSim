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

//// stave geometry
double Hstave = 1.5;  // cm
double Lstave = 27;   // cm
double Rbeampipe = 4; // cm
double x1L = -Rbeampipe-Lstave; // = -31
double x1R = -Rbeampipe;        // = -4
double x2L = +Rbeampipe;        // = +4
double x2R = +Rbeampipe+Lstave; // = +31
double yUp = +Hstave/2.;        // = +0.75
double yDn = -Hstave/2.;        // = -0.75

//// dipole geometry
double xW = 120;
double yH = 67.2;
double z1 = 102.9;
double z2 = 202.9;


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
       if (!tmp.PropagateToZBxByBz( TMath::Min(tmp.GetZ()+step, zMax), step)) break; // for different reasons the propagation may fail
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
       if (!tmp.PropagateToZBxByBz( tmp.GetZ()+zstep, zstep )) break; // for different reasons the propagation may fail
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

double m2(int i, vector<double>* px, vector<double>* py, vector<double>* pz, vector<double>* E)
{
   return (E->at(i)*E->at(i)-(px->at(i)*px->at(i)+py->at(i)*py->at(i)+pz->at(i)*pz->at(i)));
}

void WriteGeometry(vector<TPolyMarker3D*>& polm, vector<TPolyLine3D*>& poll, vector<int>& accepttrk)
{
   TCanvas* cnv_pl3d = new TCanvas("cnv_pl3d","",500,500);
   TView* view_pl3d = TView::CreateView(1);
   view_pl3d->SetRange(-80,-50,0, +80,+50,350);
   view_pl3d->ShowAxis();
   
   TCanvas* cnv_pm3d = new TCanvas("cnv_pm3d","",500,500);
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
   
   for(int i=0 ; i<(int)poll.size()    ; ++i) { /*if(accepttrk[irec])*/ cnv_pl3d->cd(); poll[i]->Draw(); }
   for(int i=0 ; i<(int)polm.size()    ; ++i) { /*if(accepttrk[irec])*/ cnv_pm3d->cd(); polm[i]->Draw(); }
   
   TLegend* leg = trkcolleg();
   cnv_pl3d->cd();
   leg->Draw("same");
   cnv_pm3d->cd();
   leg->Draw("same");
   
   cnv_pl3d->SaveAs("tracks_pl3d.root");
   cnv_pl3d->SaveAs("tracks_pl3d.pdf");
   cnv_pm3d->SaveAs("tracks_pm3d.root");
   cnv_pm3d->SaveAs("tracks_pm3d.pdf");
   
   TFile* flines = new TFile("geometry.root","RECREATE");
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


void runLUXEee(int Seed=12345, const char* setup="setup/setupLUXE.txt")
{
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
   //  det->SetMinP2Propagate(1); //NA60+
   det->SetIncludeVertex(kTRUE); // count vertex as an extra measured point
   // det->ImposeVertex(0.,0.,0.); // the vertex position is imposed NOAM
   // det->SetApplyBransonPCorrection(); // Branson correction, only relevant for setup with MS
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
 
   /// get the particles from a ttree
   TFile* fIn = new TFile("test_tony.root","READ");
   TTree* tIn = (TTree*)fIn->Get("tt");
   int nev = tIn->GetEntries();
   vector<double>* vx    = 0;
   vector<double>* vy    = 0;
   vector<double>* vz    = 0;
   vector<double>* px    = 0;
   vector<double>* py    = 0;
   vector<double>* pz    = 0;
   vector<double>* E     = 0;
   vector<int>*    pdgId = 0;
   tIn->SetBranchAddress("vx",&vx);
   tIn->SetBranchAddress("vy",&vy);
   tIn->SetBranchAddress("vz",&vz);
   tIn->SetBranchAddress("px",&px);
   tIn->SetBranchAddress("py",&py);
   tIn->SetBranchAddress("pz",&pz);
   tIn->SetBranchAddress("E",&E);
   tIn->SetBranchAddress("pdgId",&pdgId);
 
   // output tree
   gInterpreter->GenerateDictionary("vector<TLorentzVector>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>", "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",   "vector");
   int ngen = 0;
   int nres = 0;
   int nrec = 0;
   vector<TLorentzVector>  pgen;
   vector<TLorentzVector>  prec;
   vector<int>             ipair;
   vector<int>             qgen;
   vector<int>             qrec;
   vector<int>             igentrk;
   vector<TPolyMarker3D*>  polm;
   vector<TPolyLine3D*>    poll;
   vector<TPolyMarker3D*>  polm_gen;
   vector<TPolyLine3D*>    poll_gen;
   vector<int>             accepttrk;
   vector<int>     trkNITSHits;
   vector<double>  trkpT;
   vector<double>  trkSnp;
   vector<double>  trkTgl;


   TFile* fOut = new TFile("tony_rec.root","RECREATE");
   TTree* tOut = new TTree("res","res");
   tOut->Branch("ipair",&ipair);
   tOut->Branch("ngen",&ngen);
   tOut->Branch("nrec",&nrec);
   tOut->Branch("qgen",&qgen);
   tOut->Branch("pgen",&pgen);
   tOut->Branch("prec",&prec);
   tOut->Branch("igentrk",&igentrk);
   tOut->Branch("accepttrk",&accepttrk);
   tOut->Branch("polm",&polm);
   tOut->Branch("polm_gen",&polm_gen);
   tOut->Branch("poll_gen",&poll_gen);
   tOut->Branch("poll",&poll);
   tOut->Branch("trkNITSHits",&trkNITSHits);
   tOut->Branch("trkpT",&trkpT);
   tOut->Branch("trkSnp",&trkSnp);
   tOut->Branch("trkTgl",&trkTgl);

   TH1D* h_dPrel_Ruben = new TH1D("h_dPrel_Ruben",";(P_{rec}-P_{gen})/P_{gen};Tracks",500,-0.01,+0.01);
   TH1D* h_dPrel = new TH1D("h_dPrel",";(P_{rec}-P_{gen})/P_{gen};Tracks",500,-0.01,+0.01);
   TH1D* h_dErel = new TH1D("h_dErel",";(E_{rec}-E_{gen})/E_{gen};Tracks",500,-0.01,+0.01);
   TH1D* h_dxrel = new TH1D("h_dxrel",";(x_{rec}-x_{gen})/x_{gen};Tracks",200,-0.01,+0.01);
 
   /// loop on events
   for (int iev=0;iev<nev;iev++)
   {
      //// clear
      ngen = 0;    
      nres = 0;
      nrec = 0;
 	  for(int i=0;i<(int)polm_gen.size();++i) delete polm_gen[i];
 	  for(int i=0;i<(int)poll_gen.size();++i) delete poll_gen[i];
 	  for(int i=0;i<(int)polm.size();++i)     delete polm[i];
      for(int i=0;i<(int)poll.size();++i)     delete poll[i];
      pgen.clear();
      prec.clear();
      qgen.clear();
      qrec.clear();
      igentrk.clear();
      polm_gen.clear();
      poll_gen.clear();
      polm.clear();
      poll.clear();
      accepttrk.clear();
      ipair.clear();
      trkNITSHits.clear();
      trkpT.clear();
      trkSnp.clear();
      trkTgl.clear();
      
 	  //// get the next entry
 	  tIn->GetEntry(iev);
      if ((iev%outN)==0) printf("Done %d out of %d\n",iev,nev);
      int nfakeHits = 0;
 	  TLorentzVector ptmp;
 	  int ngenall = pdgId->size();
 	  int nacc = 0;
      int pair = 1;
      /// loop on particles
      for(int igen=0;igen<ngenall;igen++)
      {
         if(pdgId->at(igen)==22) continue; //// skip photons!         

         //// all the rest
         ngen++;
         int outigen = ngen-1;
         int crg = (pdgId->at(igen)==11) ? -1 : +1;
         pgen.push_back(ptmp);
         qgen.push_back(crg);
		 ipair.push_back(pair);
		 if(igen%2==0) pair++;
         if(m2(igen,px,py,pz,E)>0) pgen[outigen].SetPxPyPzE(px->at(igen), py->at(igen), pz->at(igen), E->at(igen));
         else                      pgen[outigen].SetXYZM(px->at(igen),    py->at(igen), pz->at(igen), meGeV);
         



         // prepare the probe
         bool res = det->SolveSingleTrack(pgen[outigen].Pt(),pgen[outigen].Rapidity(),pgen[outigen].Phi(), meGeV, crg, vX,vY,vZ, 0,1,99);
         if(!res) continue; // reconstruction failed
         nres++;

         // get the truth trajectory in the B-field
         KMCProbeFwd* trutrk = det->GetProbe();
         // reset the gen particle from the probe
         double pxyz_gen[3];
         double e_gen = TMath::Sqrt(trutrk->GetP()*trutrk->GetP()+meGeV*meGeV);
         trutrk->GetPXYZ(pxyz_gen);
         pgen[outigen].SetPxPyPzE(pxyz_gen[0],pxyz_gen[1],pxyz_gen[2],e_gen);
         poll_gen.push_back( TrackLine3d(trutrk,330,1,trkcol(pgen[outigen].E())) );
         polm_gen.push_back( TrackMarker3d(trutrk,0,331,1,trkcol(pgen[outigen].E())) );
         
         // get the reconstructed propagated to the vertex 
         KMCProbeFwd* trw = det->GetLayer(0)->GetWinnerMCTrack(); 
         if(!trw) break; // track was not reconstructed
         
         nrec++;         
         int irec = nrec-1;
         accepttrk.push_back(0);
         prec.push_back(ptmp);
         qrec.push_back(crg);
         igentrk.push_back(outigen);
         
         nfakeHits += trw->GetNFakeITSHits(); // in absence of background, there should be no fake hits, but we count just in case...
         double pxyz[3];
         double xyz[3];
         trw->GetPXYZ(pxyz);
         trw->GetXYZ(xyz);
         prec[irec].SetXYZM(pxyz[0],pxyz[1],pxyz[2],meGeV);

         trkNITSHits.push_back(trw->GetNITSHits());
         trkpT.push_back(trw->GetTrack()->Pt());
         trkSnp.push_back(trw->GetTrack()->GetSnp());
         trkTgl.push_back(trw->GetTrack()->GetTgl());
         h_dErel->Fill((prec[irec].E()-pgen[outigen].E())/pgen[outigen].E());
         h_dPrel->Fill((prec[irec].P()-pgen[outigen].P())/pgen[outigen].P());
         h_dPrel_Ruben->Fill(  (trw->GetP()-trutrk->GetP())/trutrk->GetP() );
         
         /// get the tracks
         poll.push_back( TrackLine3d(trw,330,1,trkcol(prec[irec].E())) );
         polm.push_back( TrackMarker3d(trw,0,331,1,trkcol(prec[irec].E())) );

         /// in acceptance?
         int nlayers = 4;
         int acctrk = 0;
         for (int i=0 ; i<polm[irec]->GetN() ; i++)
         {
            Double_t xr,yr,zr;
            polm[irec]->GetPoint(i,xr,yr,zr);
            if(zr<300) continue; //// count only the active layers
            int inacclayer = accept(xr,yr);
            acctrk += inacclayer;
           
            Double_t xg,yg,zg;
            polm_gen[outigen]->GetPoint(i,xg,yg,zg);
            h_dxrel->Fill((xr-xg)/xg);
         }
         accepttrk[irec] = (acctrk==nlayers);
         if(accepttrk[irec]) nacc++;
      }
 	 if(iev==0) WriteGeometry(polm,poll,accepttrk);
 	 fOut->cd();
 	 tOut->Fill();
   }
   fOut->cd();
   tOut->Write();
   h_dxrel->Write();
   h_dErel->Write();
   h_dPrel->Write();
   h_dPrel_Ruben->Write();
   fOut->Write();
   fOut->Close();
}