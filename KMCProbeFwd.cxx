#include "AliLog.h"
#include "KMCProbeFwd.h"
#include "KMCLayerFwd.h"
#include <TGeoGlobalMagField.h>
#include <TRandom.h>

ClassImp(KMCProbeFwd)

Int_t    KMCProbeFwd::fgNITSLayers = 0;
Double_t KMCProbeFwd::fgMissingHitPenalty = 2.;

//_______________________________________________________________________
KMCProbeFwd::KMCProbeFwd() 
  :fWeight(1) 
  ,fMass(0.10566)
  ,fChi2(0)
  ,fChi2ITS(0)
  ,fHits(0)
  ,fFakes(0)
  ,fNHits(0)
  ,fNHitsITS(0)
  ,fNHitsMS(0)
  ,fNHitsTR(0)
  ,fNHitsITSFake(0)
  ,fInnLrCheck(fgNITSLayers)
  ,fTrack()
{
}

//_______________________________________________________________________
KMCProbeFwd::KMCProbeFwd(double *xyz, double *pxyz, Int_t sign, double errLoose) 
  :fWeight(1) 
  ,fMass(0.10566)
  ,fChi2(0)
  ,fChi2ITS(0)
  ,fHits(0)
  ,fFakes(0)
  ,fNHits(0)
  ,fNHitsITS(0)
  ,fNHitsMS(0)
  ,fNHitsTR(0)
  ,fNHitsITSFake(0)
  ,fInnLrCheck(fgNITSLayers)
  ,fTrack()
{
  // create track
  Init(xyz,pxyz,sign,errLoose);
  for(int i=kMaxITSLr;i--;) fClID[i]=-2; 
  if (AliLog::GetGlobalDebugLevel()>=2) {
    AliDebug(2,Form("XYZ: %+e %+e %+e PXYZ: %+e %+e %+e Q=%d",xyz[0],xyz[1],xyz[2],pxyz[0],pxyz[1],pxyz[2],sign));
    Print("etp");
  }
}

//_______________________________________________________________________
KMCProbeFwd::KMCProbeFwd(const KMCProbeFwd& src)
:  TObject(src)
  ,fWeight(src.fWeight)
  ,fMass(src.fMass)
  ,fChi2(src.fChi2)
  ,fChi2ITS(src.fChi2ITS)
  ,fHits(src.fHits)
  ,fFakes(src.fFakes)
  ,fNHits(src.fNHits)
  ,fNHitsITS(src.fNHitsITS)
  ,fNHitsMS(src.fNHitsMS)
  ,fNHitsTR(src.fNHitsTR)
  ,fNHitsITSFake(src.fNHitsITSFake)
  ,fInnLrCheck(src.fInnLrCheck)
  ,fTrack( src.fTrack )
{
  for (int i=kMaxITSLr;i--;) fClID[i]=-src.fClID[i]; 
}

//_______________________________________________________________________
KMCProbeFwd& KMCProbeFwd::operator=(const KMCProbeFwd& src)
{
  if (this==&src) return *this;
  this->TObject::operator=(src);
  fWeight = src.fWeight;
  fMass = src.fMass;
  fChi2 = src.fChi2;
  fChi2ITS = src.fChi2ITS;
  fHits = src.fHits;
  fFakes = src.fFakes;
  fNHits = src.fNHits;
  fNHitsITS = src.fNHitsITS;
  fNHitsMS = src.fNHitsMS;
  fNHitsTR = src.fNHitsTR;
  fNHitsITSFake = src.fNHitsITSFake;
  fInnLrCheck = src.fInnLrCheck;
  fTrack = src.fTrack;
  for (int i=kMaxITSLr;i--;) fClID[i] = src.fClID[i];
  return *this;
}
  
//_______________________________________________________________________
void KMCProbeFwd::Reset() 
{
  fWeight = 1.;
  fMass=0.14; 
  fChi2=0; 
  fChi2ITS=0; 
  fHits=fFakes=0;  
  fTrack.Reset();
  ResetCovariance();
  for (int i=kMaxITSLr;i--;) fClID[i]=-2; 
  fNHits = fNHitsITS = fNHitsMS = fNHitsTR = fNHitsITSFake;
  fTrack.Reset();
}
  
//_______________________________________________________________________
Bool_t KMCProbeFwd::Init(const double *xyz, const double *pxyz, Int_t sign, double errLoose)
{
  // Init with track position/momentum in usual Lab frame
  // If errLoose>0 then scale initially small errors by this amount
  double xyzL[3],pxyzL[3];
  //printf("SetLab %f %f %f\n",pxyz[0],pxyz[1],pxyz[2]);
  Lab2Trk(xyz,xyzL);
  Lab2Trk(pxyz,pxyzL);
  double cov[21] = {1.e-6,    // assign small errors first
		    0.   ,1.e-6, 
		    0.   ,0.   ,1.e-6, 
		    0.   ,0.   ,0.   ,1.e-4,
		    0.   ,0.   ,0.   ,0.    ,1.e-4, 
		    0.   ,0.   ,0.   ,0.    ,0.   ,1.e-3};
  //
  //printf("SetTrk %f %f %f\n",pxyzL[0],pxyzL[1],pxyzL[2]);
  fTrack.Set(xyzL,pxyzL,cov,sign);
  //  printf("BefROT "); fTrack.Print();
  if (errLoose>0) fTrack.ResetCovariance(errLoose);
  //  printf("AftRes ");   fTrack.Print();
  fTrack.Rotate(pxyz[2] > 0 ? 0.:TMath::Pi());
  //  printf("AftROT ");   fTrack.Print();
  return kTRUE;
}

//__________________________________________________________________________
Int_t KMCProbeFwd::Compare(const TObject* obj) const
{
  // compare to tracks
  const KMCProbeFwd* trc = (KMCProbeFwd*) obj;
  if (trc->IsKilled()) {
    if (IsKilled()) return 0;
    return -1;
  }
  else if (IsKilled()) return 1;
  double chi2a = GetNormChi2ITS(kTRUE);
  double chi2b = trc->GetNormChi2ITS(kTRUE);
  if (chi2a<chi2b) return  -1;
  if (chi2a>chi2b) return   1;
  return 0;
}

//_______________________________________________________________________
Bool_t KMCProbeFwd::PropagateToZBxByBz(double z, double maxDZ, Double_t xOverX0, Double_t xTimesRho, Bool_t modeMC)
{
  // propagate the track to position Z in uniform material with xOverX0 rad lgt and xTimesRho lgt*density
  double zCurr = GetZ();
  double dz = z - zCurr;
  if (TMath::Abs(dz)<kAlmost0) return kTRUE;
  int nz = TMath::Abs(dz)/maxDZ + 1;
  double zstep = dz/nz;
  double xyz[3],bxyz[3],bxyzFwd[3];
  // AliDebug(2,Form("from Z=%f to Z=%f, X/X0: %f X*rho:%f, max step:%f Mode:%d (%d steps)", GetZ(),z,xOverX0,xTimesRho,maxDZ,modeMC,nz));
  for (int iz=0;iz<nz;iz++) {
    GetXYZ(xyz);              // coordinates in Lab frame
    TGeoGlobalMagField::Instance()->Field(xyz,bxyz);
    Lab2Trk(bxyz,bxyzFwd);  // field in Fwd frame
    //    printf("z = %f Field: %f %f %f\n",zCurr,bxyzFwd[0],bxyzFwd[1],bxyzFwd[2]);
    zCurr += zstep;
    if (!PropagateToZBxByBz(zCurr,bxyzFwd)) {/*printf("Fail1\n");*/ return kFALSE;}
    //if (!PropagateToZBxByBz(zCurr,bxyz)) return kFALSE;
    if (TMath::Abs(xTimesRho)>1e-6 && 
	!CorrectForMeanMaterial(xOverX0/nz, xTimesRho/nz, modeMC)) {/*printf("Fail2\n");*/ return kFALSE;}
    //	!fTrack.CorrectForMeanMaterial(xOverX0/nz, xTimesRho/nz, fMass, modeMC)) return kFALSE;
    //    fTrack.Print();
  }
  return kTRUE;
}

//_______________________________________________________________________
Bool_t KMCProbeFwd::PropagateToDCA(KMCProbeFwd* partner)
{
  // propagate the track to position Z of closest approach to partner track
  //
  double xyz[3],bxyz[3],bxyzFwd[3];
  GetXYZ(xyz);              // coordinates in Lab frame
  TGeoGlobalMagField::Instance()->Field(xyz,bxyz);
  Lab2Trk(bxyz,bxyzFwd);  // field in Fwd frame
  double zthis=0, zpartner=0;
  Double_t dca=fTrack.GetDCA(&partner->fTrack,bxyzFwd[2],zthis,zpartner);

  if (!PropagateToZBxByBz(zthis) || !partner->PropagateToZBxByBz(zpartner)) return kFALSE;
  if (AliLog::GetGlobalDebugLevel()>=2) {
    printf("Prop to DCA at %f %f | DCA = %f\n",zthis,zpartner,dca);
    Print("etp");
    partner->Print("etp");
  }
  return kTRUE;
  //
}


//_______________________________________________________________________
Double_t KMCProbeFwd::GetSigmaP2() const
{
  // error^2 on total momentum, P = sqrt(1+tgl^2)/(1/pt)
  double pinv = fTrack.GetSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double tgl  = fTrack.GetTgl();
  double tglE  = fTrack.GetSigmaTgl2();
  double pinvE = fTrack.GetSigma1Pt2();
  double pinvtgE = fTrack.GetSigma1PtTgl();
  //
  double tp12 = TMath::Sqrt(1.+tgl*tgl);
  double dt =  tgl/pinv/tp12;  // dP/dtgl
  double dc = -tp12/pinv/pinv; // dP/dC
  double err2 = dt*dt*tglE +dc*dc*pinvE + dt*dc*pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
Double_t KMCProbeFwd::GetSigmaPX2() const
{
  // error^2 on Px in Lav frame (Z in the tracking frame, = tgl/(1/pt)
  double pinv = fTrack.GetSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  double tgl  = fTrack.GetTgl();
  double tglE  = fTrack.GetSigmaTgl2();
  double pinvE = fTrack.GetSigma1Pt2();
  double pinvtgE = fTrack.GetSigma1PtTgl();
  //
  double dt =  1./pinv;      // dP/dtgl
  double dc = -tgl/pinv/pinv; // dP/dC
  double err2 = dt*dt*tglE +dc*dc*pinvE + dt*dc*pinvtgE;
  return err2;
  //
}

//_______________________________________________________________________
Double_t KMCProbeFwd::GetSigmaPY2() const
{
  // error^2 on Py in Lab frame (Y in the tracking frame, = sinp/(1/pt)
  double pinv = fTrack.GetSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  Double_t cosAlp=TMath::Cos(fTrack.GetAlpha()), sinAlp=TMath::Sin(fTrack.GetAlpha());
  double snp  = fTrack.GetSnp();
  double csp  = TMath::Sqrt((1. - snp)*(1. + snp));
  double snpE  = fTrack.GetSigmaSnp2();
  double pinvE = fTrack.GetSigma1Pt2();
  double pinvsnpE = fTrack.GetSigma1PtSnp();
  //
  double ds =  (cosAlp-sinAlp*snp/csp)/pinv; // dP/dsnp
  double dc = -(snp*cosAlp+csp*sinAlp)/pinv/pinv; // dP/dC
  double err2 = ds*ds*snpE +dc*dc*pinvE + ds*dc*pinvsnpE;
  return err2;

  //
}

//_______________________________________________________________________
Double_t KMCProbeFwd::GetSigmaPZ2() const
{
  // error^2 on Pz in Lab frame (X in the tracking frame, = (sqrt(1-snp^2)*cosAlp-snp*sinAlp)/(1/pt)
  double pinv = fTrack.GetSigned1Pt();
  if (TMath::Abs(pinv)<kAlmost0) return 0;
  Double_t cosAlp=TMath::Cos(fTrack.GetAlpha()), sinAlp=TMath::Sin(fTrack.GetAlpha());
  double snp  = fTrack.GetSnp();
  double snpE  = fTrack.GetSigmaSnp2();
  double pinvE = fTrack.GetSigma1Pt2();
  double pinvsnpE = fTrack.GetSigma1PtSnp();
  //
  double csp = TMath::Sqrt( (1.-snp)*(1.+snp) );
  double ds = -snp/pinv*(cosAlp/csp + sinAlp);        // dP/dsnp
  double dc = -(csp*cosAlp-snp*sinAlp)/pinv/pinv;     // dP/dC
  double err2 = ds*ds*snpE +dc*dc*pinvE + ds*dc*pinvsnpE;
  return err2;
  //
}

//_______________________________________________________________________
void KMCProbeFwd::ResetCovariance(float err)
{
  // reset cov matrix
  double *trCov  = (double*)fTrack.GetCovariance();
  double *trPars = (double*)fTrack.GetParameter();
  const double kLargeErr2Coord = 50*50;
  const double kLargeErr2Dir = 0.6*0.6;
  const double kLargeErr2PtI = 2.;
  for (int ic=15;ic--;) trCov[ic] = 0.;
  trCov[kY2]   = trCov[kZ2]   = err<0 ? kLargeErr2Coord : err*err; 
  trCov[kSnp2] = trCov[kTgl2] = kLargeErr2Dir;
  trCov[kPtI2] = kLargeErr2PtI*trPars[kPtI]*trPars[kPtI];
  fTrack.CheckCovariance();
}

//_________________________________________________________
Bool_t KMCProbeFwd::ApplyMSEL(double x2X0, double xTimesRho)
{
  // simulate random modification of track params due to the MS
  // note: here we work directly in TrackPar internal coordinates
  if (x2X0<0 || xTimesRho>0) AliFatal(Form("Negative X/X0=%f or X*rho=%f",x2X0,xTimesRho));
  if (x2X0<1e-7) return kTRUE;
  //  printf("Before %e: ",x2X0); fTrack.Print();
  double alpha = fTrack.GetAlpha(); // store original alpha
  double snp = fTrack.GetSnp();
  double dip = fTrack.GetTgl();
  Double_t angle=TMath::Sqrt((1.+ dip*dip)/((1-snp)*(1.+snp)));
  x2X0 *= angle;
  //
  static double covCorr[15],covDum[21]={0};
  static double mom[3],pos[3];
  double *cov = (double*) fTrack.GetCovariance();
  memcpy(covCorr,cov,15*sizeof(double));
  fTrack.GetXYZ(pos);
  fTrack.GetPxPyPz(mom);
  double pt2 = mom[0]*mom[0]+mom[1]*mom[1];
  double pt = TMath::Sqrt(pt2);
  double ptot2 = pt2 + mom[2]*mom[2];
  double ptot  = TMath::Sqrt(ptot2);
  double beta = ptot/TMath::Sqrt(ptot2 + fMass*fMass);
  double sigth = TMath::Sqrt(x2X0)*0.014/(ptot*beta);
  //
  // a la geant
  double phiSC = gRandom->Rndm()*TMath::Pi();
  double thtSC = gRandom->Gaus(0,1.4142*sigth);
  //  printf("MS phi: %+.5f tht: %+.5f\n",phiSC,thtSC);
  double sn = TMath::Sin(thtSC);
  double dx = sn*TMath::Sin(phiSC);
  double dy = sn*TMath::Cos(phiSC);  
  double dz = TMath::Cos(thtSC);
  double v[3];
  //  printf("Before: %+.3e %+.3e %+.3e | MS: %+.3e %+.3e\n",mom[0],mom[1],mom[2],thtSC,phiSC);
  for (int i=3;i--;) mom[i] /= ptot;
  double vmm = TMath::Sqrt(mom[0]*mom[0]+mom[1]*mom[1]);
  if (!IsZero(pt)) {
    double pd1 = mom[0]/vmm;
    double pd2 = mom[1]/vmm;
    v[0] = pd1*mom[2]*dx - pd2*dy + mom[0]*dz;
    v[1] = pd2*mom[2]*dx + pd1*dy + mom[1]*dz;
    v[2] = -vmm*dx                + mom[2]*dz;
  }
  else {
    v[0] = dx;
    v[1] = dy;
    v[2] = dz*TMath::Sign(1.,mom[2]);
  }
  double nrm = TMath::Sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  //  printf("before :%+e %+e %+e  || %+e %+e %+e %+e\n",mom[0],mom[1],mom[2],  sigth, x2X0, pt, beta);
  //  trc->Print();
  // direction cosines -> p
  for (int i=3;i--;) mom[i] = ptot*v[i]/nrm;
  //  printf("After : %+.3e %+.3e %+.3e\n",mom[0],mom[1],mom[2]);
  fTrack.Set(pos,mom,covDum,fTrack.Charge());
  //
  if (!fTrack.Rotate(alpha)) return kFALSE;
  //
  if (xTimesRho<-1e-7) {
    // Apply ELoss
    double p = fTrack.GetP();
    Double_t bg = p/fMass;
    Double_t dEmean=TrackPar::BetheBlochSolid(bg)*xTimesRho; // mean loss
    const Double_t knst=0.07; // To be tuned.  
    Double_t sigmadE=knst*TMath::Sqrt(TMath::Abs(dEmean));  // fluctuation
    Double_t dE = 0;
    while( (dE=gRandom->Gaus(dEmean,sigmadE))>=0) {}
    double p2 = p*p;
    Double_t e=TMath::Sqrt(p2 + fMass*fMass);
    if ( TMath::Abs(dE) > 0.3*e ) return kFALSE; //30% energy loss is too much!
    double fc = (1.+ dE/p2*(dE + 2*e));
    if ( fc < 0. ) return kFALSE;
    fc = 1./TMath::Sqrt(fc);
    if (TMath::Abs(fTrack.Get1P()*fc)>100.) return kFALSE; //Do not track below 10 MeV/c
    //
    double *pars = (double*)fTrack.GetParameter();
    double errcor = (sigmadE*e/p2*pars[4]);
    pars[4] *= fc;
    covCorr[14] = errcor*errcor;
  }
  // 
  memcpy(cov,covCorr,15*sizeof(double));
  //printf("After %e: ",x2X0); fTrack.Print();
  return kTRUE;
  //
}

//_______________________________________________________________________
void KMCProbeFwd::Print(Option_t* opt) const
{
  printf("Killed:%d M:%.3f Chi2:%.3f(%.3f) Chi2ITS:%.3f(%.3f) | Wgh:%.2e\nHits: Tot:%2d ITS:%d ITSFake:%d (last check:%d)",
	 IsKilled(),fMass,GetChi2(),GetNormChi2(1),
	 GetChi2ITS(),GetNormChi2ITS(1),GetWeight(),
	 fNHits,fNHitsITS,fNHitsITSFake,fInnLrCheck);
  // hit pattern
  printf(" Pattern: |");
  for (int i=0;i<fgNITSLayers;i++) {
    if (!(fHits&(0x1<<i))) printf(".");
    else if ( fFakes&(0x1<<i) ) printf("-");
    else printf("+");
  }
  printf("|  ");
  TString opts = opt; opts.ToLower();
  if (opts.Contains("clid")) {
    for (int ilr=0;ilr<fgNITSLayers;ilr++) {
      printf("%4d|",fClID[ilr]);
    }
  }
  printf("\n");
  if (opts.Contains("etp")) fTrack.Print();
}

//_______________________________________
void KMCProbeFwd::AddHit(const KMCLayerFwd* lr, double chi2, Int_t clID) {
  // note: lr is active layer ID
  if (!lr) return;
  fNHits++;
  fChi2 += chi2;
  int lrID = lr->GetActiveID();
  if (lr->IsITS()) {
    fChi2ITS += chi2;
    SetWBit(fHits,lrID); 
    fNHitsITS++;
    if (clID>-1) {
      SetWBit(fFakes,lrID);
      fNHitsITSFake++;
    }
    fClID[lrID] = clID;
    //else ResetWBit(fFakes,lr);
  }
  else if (lr->IsMS()) fNHitsMS++;
  else if (lr->IsTrig()) fNHitsTR++;
  else {
    lr->Print();
    AliFatal("Invalid layer type");
  }
}

//____________________________________________
void KMCProbeFwd::ImposeKinematics(const double* xyzLab,const double* cosinesLab,
				   double en, double mass, int charge) 
{
  // RS: note: we assume p=e here to avoid problem with e+ e- interpreted as muon
  double p = en*en - mass*mass;
  if (p<=0) {
    printf("Anomalous kinematics: E:%e M:%e",en,mass);
    exit(1);
  }
  p = TMath::Sqrt(p);
  double pxyz[3] = {p*cosinesLab[0],p*cosinesLab[1],p*cosinesLab[2]};
  //  printf("Imp : %f %f %f | %f\n",pxyz[0],pxyz[0],pxyz[0], en);
  Init(xyzLab,pxyz,charge,1.e4);
  SetMass(mass);
  ResetCovariance();// reset cov.matrix
  //  printf("set at %f %f %f \n",xyzLab[0],xyzLab[1],xyzLab[2]);
}
