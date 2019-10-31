#include "KMCClusterFwd.h"
#include <TString.h>

ClassImp(KMCClusterFwd)

//_________________________________________________________________________
KMCClusterFwd::KMCClusterFwd(KMCClusterFwd &src) 
:  TObject(src)
  ,fX(src.fX)
  ,fY(src.fY)
  ,fZ(src.fZ)
{}

//__________________________________________________________________________
KMCClusterFwd& KMCClusterFwd::operator=(const KMCClusterFwd& src) 
{
  if (this!=&src) {
    TObject::operator=(src);
    fX = src.fX;
    fY = src.fY;
    fZ = src.fZ;
  }
  return *this;
}

//_________________________________________________________________________
void KMCClusterFwd::Print(Option_t *opt) const 
{
  TString opts = opt;
  opts.ToLower();
  if (opts.Contains("lc")) 
    printf("Tr#%4d Loc (%+.4e,%+.4e %+.4e) %s",GetTrID(),fX,fY,fZ,IsKilled()?"Killed":""); 
  else 
    printf("Tr#%4d Lab (%+.4e,%+.4e %+.4e) %s",GetTrID(),GetXLab(),GetYLab(),GetZLab(),IsKilled()?"Killed":""); 
  if (opts.Contains("nl")) return;
  printf("\n");
}

//_________________________________________________________________________
Bool_t KMCClusterFwd::IsEqual(const TObject* obj) const 
{
  // check if clusters are equal
  KMCClusterFwd* cl = (KMCClusterFwd*)obj;
  const double kTiny = 1e-12;
  if (TMath::Abs(GetXLab()-cl->GetXLab())>kTiny || 
      TMath::Abs(GetYLab()-cl->GetYLab())>kTiny) return kFALSE;
  return kTRUE;
}

//_________________________________________________________________________
Int_t KMCClusterFwd::Compare(const TObject* obj) const
{
  // compare 1st labx, then laby
  KMCClusterFwd* cl = (KMCClusterFwd*)obj;
  if (GetXLab() > cl->GetXLab()) return -1; // tracking -Z = labX
  if (GetXLab() < cl->GetXLab()) return  1;
  if (GetYLab() < cl->GetYLab()) return -1; // tracking Y = labY
  if (GetYLab() > cl->GetYLab()) return  1;
  return 0;
}
