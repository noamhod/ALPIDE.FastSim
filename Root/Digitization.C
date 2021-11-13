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

TMapTSTH1D histos1;
TMapTSTH2D histos2;

TString storage =  gSystem->ExpandPathName("$STORAGEDIR");

KMCDetectorFwd* det = 0;
vector<double>* zlayer = new vector<double>;
double meMeV = 0.5109989461; //MeV
double meGeV = meMeV/1000.;

//// staves geometry
double Hstave = 1.5;  // cm
double Lstave = 50; //27.12;   // cm
double Rbeampipe = 2.413; // cm
double yUp = +Hstave/2.;
double yDn = -Hstave/2.;

//// dipole geometry
double xW = 33.0;
double yH = 10.8;
double z1 = 100;
double z2 = 202.9;

double zL1I = -999;
double zL1O = -999;

double zL2I = -999;
double zL2O = -999;

double zL3I = -999;
double zL3O = -999;

double zL4I = -999;
double zL4O = -999;

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

vector<TString> layernames = {"L1I","L1O", "L2I","L2O", "L3I","L3O", "L4I","L4O"};

void setParametersFromDet()
{
	cout << "====================================" << endl;
	cout << "============DEFINITIONS+============" << endl;
	cout << "====================================" << endl;
	KMCLayerFwd* layer_outer = det->GetLayer("L1O");
	KMCLayerFwd* layer_inner = det->GetLayer("L1I");
	
	Hstave = layer_outer->GetYMax()-layer_outer->GetYMin();
	Lstave = layer_outer->GetXMaxP()-layer_outer->GetXMinP();
	cout << "Hstave=" << Hstave << ", Lstave=" << Lstave << endl;
	
	xMinEI = layer_inner->GetXMinE();
	xMinEO = layer_outer->GetXMinE();
	xMinPI = layer_inner->GetXMinP();
	xMinPO = layer_outer->GetXMinP();
	cout << "xMinEI=" << xMinEI << ", xMinEO=" << xMinEO << ", xMinPI=" << xMinPI << ", xMinPO=" << xMinPO << endl;
	
	xMaxEI = layer_inner->GetXMaxE();
	xMaxEO = layer_outer->GetXMaxE();
	xMaxPI = layer_inner->GetXMaxP();
	xMaxPO = layer_outer->GetXMaxP();
	cout << "xMaxEI=" << xMaxEI << ", xMaxEO=" << xMaxEO << ", xMaxPI=" << xMaxPI << ", xMaxPO=" << xMaxPO << endl;
		
	yUp = +Hstave/2.;
	yDn = -Hstave/2.;
	cout << "yUp=" << yUp << ", yDn=" << yDn << endl;
	
	//// get the Bfield from the setup
	TVirtualMagField* fld = TGeoGlobalMagField::Instance()->GetField();
	MagField* fldm = (MagField*) fld;
	const double* BfieldObj = fldm->GetBVals(0);    // region 0
	const double* BfieldXminObj = fldm->GetXMin(); // region 0
	const double* BfieldXmaxObj = fldm->GetXMax(); // region 0
	const double* BfieldYminObj = fldm->GetYMin(); // region 0
	const double* BfieldYmaxObj = fldm->GetYMax(); // region 0
	const double* BfieldZminObj = fldm->GetZMin(); // region 0
	const double* BfieldZmaxObj = fldm->GetZMax(); // region 0
	double BfieldValTesla = BfieldObj[1]/10; /// the B dield is only in the y direction (0,B,0), hence the index 1
	double BfieldXmin = BfieldXminObj[0]; // region 0
   double BfieldXmax = BfieldXmaxObj[0]; // region 0
   double BfieldYmin = BfieldYminObj[0]; // region 0
   double BfieldYmax = BfieldYmaxObj[0]; // region 0
   double BfieldZmin = BfieldZminObj[0]; // region 0
   double BfieldZmax = BfieldZmaxObj[0]; // region 0
	cout << "BfieldValTesla=" << BfieldValTesla << ", BfieldXmin=" << BfieldXmin << ", BfieldXmax=" << BfieldXmax << ", BfieldYmin=" << BfieldYmin << ", BfieldYmax=" << BfieldYmax << ", BfieldZmin=" << BfieldZmin << ", BfieldZmax=" << BfieldZmax << endl;
	
	xW = BfieldXmax-BfieldXmin;
	yH = BfieldYmax-BfieldYmin;
	z1 = BfieldZmin;
	z2 = BfieldZmax;
	cout << "xW=" << xW << ", yH=" << yH << ", z1=" << z1 << ", z2=" << z2 << endl;
	
	zL1I = det->GetLayer("L1I")->GetZ();
	zL1O = det->GetLayer("L1O")->GetZ();
	cout << "zL1I=" << zL1I << ", zL1O=" << zL1O << endl;

	zL2I = det->GetLayer("L2I")->GetZ();
	zL2O = det->GetLayer("L2O")->GetZ();
	cout << "zL2I=" << zL2I << ", zL2O=" << zL2O << endl;

	zL3I = det->GetLayer("L3I")->GetZ();
	zL3O = det->GetLayer("L3O")->GetZ();
	cout << "zL3I=" << zL3I << ", zL3O=" << zL3O << endl;

	zL4I = det->GetLayer("L4I")->GetZ();
	zL4O = det->GetLayer("L4O")->GetZ();
	cout << "zL4I=" << zL4I << ", zL4O=" << zL4O << endl;
	
	zLastLayer  = zL4I;
	zFirstLayer = zL1O;
	cout << "zLastLayer=" << zLastLayer << ", zFirstLayer=" << zFirstLayer << endl;
	
	cout << "====================================" << endl;
	cout << "====================================" << endl;
}

int acceptcls(double x, double y, double z, double step=0.1)
{
	if( (abs(z-zL1I)>step && abs(z-zL1O)>step) && (abs(z-zL2I)>step && abs(z-zL2O)>step) && (abs(z-zL3I)>step && abs(z-zL3O)>step) && (abs(z-zL4I)>step && abs(z-zL4O)>step) ) return 0;
	if(x<xMinEO || x>xMaxPO) return 0;
	if(x>xMaxEI && x<xMinPI) return 0;
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
    TPolyMarker3D *polymarker = new TPolyMarker3D(zlayer->size());
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
   view_pl3d->SetRange(-80,-50,0, +80,+50,350);
   view_pl3d->ShowAxis();
   
   TCanvas* cnv_pm3d = new TCanvas("cnv_pm3d"+suff,"",500,500);
   TView* view_pm3d = TView::CreateView(1);
   view_pm3d->SetRange(-80,-50,0, +80,+50,350);
   view_pm3d->ShowAxis();
   
   TPolyLine3D* stave1PI = GetLayer("P","I",zL1I,kGreen+3);
   TPolyLine3D* stave1PO = GetLayer("P","O",zL1O,kGreen+3);
   TPolyLine3D* stave1EI = GetLayer("E","I",zL1I,kGreen+3);
   TPolyLine3D* stave1EO = GetLayer("E","O",zL1O,kGreen+3);
   TPolyLine3D* stave2PI = GetLayer("P","I",zL2I,kGreen+3);
   TPolyLine3D* stave2PO = GetLayer("P","O",zL2O,kGreen+3);
   TPolyLine3D* stave2EI = GetLayer("E","I",zL2I,kGreen+3);
   TPolyLine3D* stave2EO = GetLayer("E","O",zL2O,kGreen+3);
   TPolyLine3D* stave3PI = GetLayer("P","I",zL3I,kGreen+3);
   TPolyLine3D* stave3PO = GetLayer("P","O",zL3O,kGreen+3);
   TPolyLine3D* stave3EI = GetLayer("E","I",zL3I,kGreen+3);
   TPolyLine3D* stave3EO = GetLayer("E","O",zL3O,kGreen+3);
   TPolyLine3D* stave4PI = GetLayer("P","I",zL4I,kGreen+3);
   TPolyLine3D* stave4PO = GetLayer("P","O",zL4O,kGreen+3);
   TPolyLine3D* stave4EI = GetLayer("E","I",zL4I,kGreen+3);
   TPolyLine3D* stave4EO = GetLayer("E","O",zL4O,kGreen+3);
   TPolyLine3D* dipole  = GeDipole(kGray);
   
   cnv_pl3d->cd();
   dipole->Draw();
   stave1PI->Draw();
   stave1PO->Draw();
   stave1EI->Draw();
   stave1EO->Draw();
   stave2PI->Draw();
   stave2PO->Draw();
   stave2EI->Draw();
   stave2EO->Draw();
   stave3PI->Draw();
   stave3PO->Draw();
   stave3EI->Draw();
   stave3EO->Draw();
   stave4PI->Draw();
   stave4PO->Draw();
   stave4EI->Draw();
   stave4EO->Draw();
   
   cnv_pm3d->cd();
   dipole->Draw();
   stave1PI->Draw();
   stave1PO->Draw();
   stave1EI->Draw();
   stave1EO->Draw();
   stave2PI->Draw();
   stave2PO->Draw();
   stave2EI->Draw();
   stave2EO->Draw();
   stave3PI->Draw();
   stave3PO->Draw();
   stave3EI->Draw();
   stave3EO->Draw();
   stave4PI->Draw();
   stave4PO->Draw();
   stave4EI->Draw();
   stave4EO->Draw();
   
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
   stave1PI->Write();
   stave1PO->Write();
   stave1EI->Write();
   stave1EO->Write();
   stave2PI->Write();
   stave2PO->Write();
   stave2EI->Write();
   stave2EO->Write();
   stave3PI->Write();
   stave3PO->Write();
   stave3EI->Write();
   stave3EO->Write();
   stave4PI->Write();
   stave4PO->Write();
   stave4EI->Write();
   stave4EO->Write();
   leg->Write();
   // flines->Write();
   flines->Close();
}


bool accepttrk(vector<TVector3>& clusters, bool fullacc, double step=0.1, int nMinLayers=3)
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
void SetOutTree(TString fOutName)
{
   fOut = new TFile(fOutName,"RECREATE");
   tOut = new TTree("dig","dig");
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
}
void KillOutTree()
{
   fOut->cd();
   tOut->Write();
	for(TMapTSTH1D::iterator it=histos1.begin() ; it!=histos1.end() ; ++it) it->second->Write();
	for(TMapTSTH2D::iterator it=histos2.begin() ; it!=histos2.end() ; ++it) it->second->Write();
   fOut->Write();
   fOut->Close();
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


void AddCluster(int slvidx, int index_offset, TString process, TString LYR)
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
	/// !!! dirty fix !!! ///
	/////////////////////////
	unsigned int ncls = clusters_r[slvidx].size();
	double xprev = (ncls>0) ? clusters_r[slvidx][ncls-1].X() : x;
	if(abs(x-xprev)>2) return;
	if(x*xprev<0) return;
	
   clusters_xyz[slvidx]->SetNextPoint(x,y,z);
   clusters_r[slvidx].push_back( v );
	clusters_id[slvidx].push_back( layerid*index_offset+slvidx ); // assuming no chance to have >index_offset tracks
	clusters_layerid[slvidx].push_back( layerid );
	clusters_type[slvidx].push_back( (process.Contains("bkg")) ? 0 : 1 );
	
	cluster->Reset();
	// det->GetLayer(LYR)->GetMCCluster()->Kill();
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
	if(argc==3 and !((TString)argv[2]).Contains("-evnt=")) { printf("argc=3 but cannot parse %s\n",argv[2]); exit(-1); }
	if(argc==4 and !((TString)argv[3]).Contains("-seed=")) { printf("argc=4 but cannot parse %s\n",argv[3]); exit(-1); }
	//// assign inputs
	TString process = ((TString)argv[1]).ReplaceAll("-proc=",""); // mandatory
	int     evnt    = (argc>2) ? toint(((TString)argv[2]).ReplaceAll("-evnt=","")) : -1; // job id [optional]
	int     Seed    = (argc>3) ? toint(((TString)argv[3]).ReplaceAll("-seed=","")) : 12345; // seed [optional]
	//// print assigned inputs
	cout << "process=" << process << endl;
	cout << "evnt=" << evnt << endl;
	cout << "Seed=" << Seed << endl;
	
	
	TString eventid = (evnt<0) ? "" : FormatEventID(evnt);
	TString proc = process;
	proc.ReplaceAll("_bkg","");
	TString setup = "../setup/setupLUXE_"+proc+".txt";
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
	
	
	///////////////////////////
	setParametersFromDet(); ///
	///////////////////////////
	
   
   zlayer->push_back(0);     // IP (vertex)
   zlayer->push_back(z1);    // start of dipol
   zlayer->push_back(z2);    // end of dipol
   zlayer->push_back(zL1I); // L1 inner
   zlayer->push_back(zL1O); // L1 outer
   zlayer->push_back(zL2I); // L2 inner
   zlayer->push_back(zL2O); // L2 outer
   zlayer->push_back(zL3I); // L3 inner
   zlayer->push_back(zL3O); // L3 outer
   zlayer->push_back(zL4I); // L4 inner
   zlayer->push_back(zL4O); // L4 outer



   int outN = 100;
   int nMaxEventsPerFile = 100;
   int nFiles = 1;

   // TString process = "bppp";  /// trident or bppp or bppp_bkg or trident_bkg
	// if(process.Contains("trident")) resetToTridentGeometry();
	
	int index_offset = (process.Contains("bkg")) ? 10000 : 100000;

   /// get the particles from a ttree
   TFile* fIn = new TFile(storage+"/data/root/raw_"+process+".root","READ");
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
   gInterpreter->GenerateDictionary("vector<TLorentzVector>",       "vector");
   gInterpreter->GenerateDictionary("vector<TPolyMarker3D*>",       "vector");
   gInterpreter->GenerateDictionary("vector<TPolyLine3D*>",         "vector");
	gInterpreter->GenerateDictionary("vector<vector<int> >",         "vector");
	gInterpreter->GenerateDictionary("vector<vector<TVector3> >",    "vector");
	gInterpreter->GenerateDictionary("vector<vector<TPolyLine3D> >", "vector");
	gSystem->Exec("mkdir -p "+storage+"/data/root/dig");
	TString fOutName = storage+"/data/root/dig/dig_"+process+"_"+eventid+".root";
	SetOutTree(fOutName);
	
	// TString hname = "";
	// hname = "h2_z_vs_x"; histos2.insert( make_pair(hname, new TH2D(hname,";x [cm];z [cm];Tracks",1000,-100,+100, 2000,0,+400)) );
	// hname = "h2_z_vs_y"; histos2.insert( make_pair(hname, new TH2D(hname,";y [cm];z [cm];Tracks",1000,-100,+100, 2000,0,+400)) );
	
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
      // int pair = 1;
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
			vector<TVector3> vtmpTVector3d;
         int q = (pdgId->at(igen)==11) ? -1 : +1;
			ptmp.SetXYZM(px->at(igen), py->at(igen), pz->at(igen), meGeV);
			
			///////////////////////////////////
			/// look just on positrons for now
			// if(q<0) continue; /////////////////
			///////////////////////////////////
			
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
			if(!process.Contains("bkg")) trkpts_fullrange.push_back( TrackMarker3d(trutrk,0,zLastLayer+5,0.1,trkcol(ptmp.E()),false,true) );
			trkpts.push_back( TrackMarker3d(trutrk,0,zLastLayer+5,0.1,trkcol(ptmp.E())) );
			trklin.push_back( TrackLine3d(trutrk,zLastLayer+5,1,trkcol(ptmp.E())) );
			
			// for(Int_t n=0 ; n<trkpts_fullrange[slvidx]->GetN() ; ++n)
			// {
			// 	Double_t xp, yp, zp;
			// 	trkpts_fullrange[slvidx]->GetPoint(n,xp,yp,zp);
			// 	histos2["h2_z_vs_x"]->Fill(xp,zp);
			// 	histos2["h2_z_vs_y"]->Fill(yp,zp);
			// }
         
			clusters_id.push_back( vtmp );
			clusters_layerid.push_back( vtmp );
			clusters_type.push_back( vtmp );
			clusters_xyz.push_back( new TPolyMarker3D() );
			clusters_r.push_back( vtmpTVector3d );
			acc.push_back( 0 );

			
			// get the reconstructed propagated to the vertex
			for(unsigned int k=0 ; k<layernames.size() ; k++)
			{
				TString LYR = layernames[k];
				if(!det->GetLayer(LYR)->IsITS()) continue;
				AddCluster(slvidx,index_offset,process,LYR);
			}
			int nclusters = clusters_layerid[slvidx].size(); // same as layers hit by the track
			cout << "slvidx=" << slvidx << ", E=" << trkp4[slvidx].E() << ", Q=" << crg[slvidx] << " --> Nclusters=" << nclusters << endl;
			for(int j=0 ; j<nclusters ; ++j)
			{
			   Double_t x,y,z;
			   // clusters_xyz[slvidx]->GetPoint(j,x,y,z);
			   x = clusters_r[slvidx][j].X();
			   y = clusters_r[slvidx][j].Y();
			   z = clusters_r[slvidx][j].Z();
				int layerid = clusters_layerid[slvidx][j];
				cout << "point " << j << ", layerid=" << layerid << ", layername=" << det->GetLayer(layerid)->GetName() << ", xyz=(" << x << "," << y << "," << z << ")" << endl;
				// det->GetLayer(layerid)->Print();
			}
			
			
			
			/// check acceptance
			// acc[slvidx] = (accepttrk(clusters_xyz[slvidx],false) && acceptpts(trkpts[slvidx],false));
			// acc[slvidx] = (accepttrk(clusters_xyz[slvidx],false));
			acc[slvidx] = (accepttrk(clusters_r[slvidx],false));
			if(acc[slvidx]) nacc++;
			
			
			
      }
		if(iev==0) WriteGeometry(trkpts,trklin,process,acc,clusters_xyz,"_truth");
		if(iev%1==0) cout << "iev=" << iev << " --> ngen=" << ngen << ", nslv=" << nslv << ", nacc=" << nacc << endl;
      if(nslv!=ngen and !process.Contains("bkg")) cout << "Warning: nslv=" << nslv << ", ngen=" << ngen << " --> problem" << endl;
		
		/// fill the tree
      fOut->cd();
      tOut->Fill();
		
		int nevprocessed = iev+1; // iev starts from 0
		if(process.Contains("bkg") and fullloop)
		{
			if(nevprocessed%nMaxEventsPerFile==0)
			{
				cout << "Killing file #" << nFiles << " and redefining new one" << endl;
		   	KillOutTree();
				
				cout << "Renaming file #" << nFiles << " to include its index" << endl;
				RenameOutTree(fOutName,nFiles);
				nFiles++; // must propagate!
				
				// reset the file and tree as usual
				cout << "Resetting file #" << nFiles << " and tree" << endl;
				SetOutTree(fOutName);
			}
		}
			
   }
   KillOutTree();
	if(process.Contains("bkg") and eventid=="") RenameOutTree(fOutName,nFiles);
	
	return 0;
}
