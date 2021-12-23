#include "LuxeLabels.h"

#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TPad.h"
#include "TMarker.h"


void LUXELabel(Double_t x,Double_t y,const char* text,Color_t color){
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(52);
  l.SetTextColor(color);

  double delx = 0.85*0.115*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"LUXE");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }

}

