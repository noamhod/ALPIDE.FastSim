//
// LUXE Style based on ATLAS Style, based on a style file from BaBar
//

#include <iostream>

#include "LuxeStyle.h"

#include "TROOT.h"

void SetLuxeStyle ()
{
  static TStyle* luxeStyle = 0;
  std::cout << "\nApplying LUXE style settings...\n" << std::endl ;
  if ( luxeStyle==0 ) luxeStyle = LuxeStyle();
  gROOT->SetStyle("LUXE");
  gROOT->ForceStyle();
}

TStyle* LuxeStyle() 
{
  TStyle *luxeStyle = new TStyle("LUXE","Luxe style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  luxeStyle->SetFrameBorderMode(icol);
  luxeStyle->SetFrameFillColor(icol);
  luxeStyle->SetCanvasBorderMode(icol);
  luxeStyle->SetCanvasColor(icol);
  luxeStyle->SetPadBorderMode(icol);
  luxeStyle->SetPadColor(icol);
  luxeStyle->SetStatColor(icol);
  //luxeStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  luxeStyle->SetPaperSize(20,26);

  // set margin sizes
  luxeStyle->SetPadTopMargin(0.05);
  luxeStyle->SetPadRightMargin(0.05);
  luxeStyle->SetPadBottomMargin(0.16);
  luxeStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  luxeStyle->SetTitleXOffset(1.4);
  luxeStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05;
  luxeStyle->SetTextFont(font);

  luxeStyle->SetTextSize(tsize);
  luxeStyle->SetLabelFont(font,"x");
  luxeStyle->SetTitleFont(font,"x");
  luxeStyle->SetLabelFont(font,"y");
  luxeStyle->SetTitleFont(font,"y");
  luxeStyle->SetLabelFont(font,"z");
  luxeStyle->SetTitleFont(font,"z");
  
  luxeStyle->SetLabelSize(tsize,"x");
  luxeStyle->SetTitleSize(tsize,"x");
  luxeStyle->SetLabelSize(tsize,"y");
  luxeStyle->SetTitleSize(tsize,"y");
  luxeStyle->SetLabelSize(tsize,"z");
  luxeStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  luxeStyle->SetMarkerStyle(20);
  luxeStyle->SetMarkerSize(1.2);
  luxeStyle->SetHistLineWidth(2.);
  luxeStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars (as recommended in LUXE figure guidelines)
  //luxeStyle->SetErrorX(0.0001); // this prevents the E2 draw option from working, use X0 option instead
  // get rid of error bar caps
  luxeStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  luxeStyle->SetOptTitle(0);
  //luxeStyle->SetOptStat(1111);
  luxeStyle->SetOptStat(0);
  //luxeStyle->SetOptFit(1111);
  luxeStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  luxeStyle->SetPadTickX(1);
  luxeStyle->SetPadTickY(1);

  //Legend style
  luxeStyle->SetLegendBorderSize(0);
  luxeStyle->SetLegendFillColor(0);
  luxeStyle->SetLegendFont(42);
  luxeStyle->SetLegendTextSize(0.);

  return luxeStyle;

}

