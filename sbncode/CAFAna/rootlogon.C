// To set this as default, you need a .rootrc file in your home directory,
// containing the following line:
// Rint.Logon: /full/path/to/rootlogon.C

#ifndef ROOTLOGON_C
#define ROOTLOGON_C

#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"

void rootlogon()
{
  printf("Welcome to the ROOT of all evils \n");

  // Defaults to classic style, but that's OK, we can fix it
  TStyle* sbnStyle = new TStyle("sbnStyle", "SBN Style");

  // Centre title
  sbnStyle->SetTitleAlign(22);
  sbnStyle->SetTitleX(.5);
  sbnStyle->SetTitleY(.95);
  sbnStyle->SetTitleBorderSize(0);

  // No info box
  sbnStyle->SetOptStat(0);

  //set the background color to white
  sbnStyle->SetFillColor(10);
  sbnStyle->SetFrameFillColor(10);
  sbnStyle->SetCanvasColor(10);
  sbnStyle->SetPadColor(10);
  sbnStyle->SetTitleFillColor(0);
  sbnStyle->SetStatColor(10);

  // Don't put a colored frame around the plots
  sbnStyle->SetFrameBorderMode(0);
  sbnStyle->SetCanvasBorderMode(0);
  sbnStyle->SetPadBorderMode(0);

  // Set the default line color for a fit function to be red
  sbnStyle->SetFuncColor(kRed);

  // Marker settings
  //  sbnStyle->SetMarkerStyle(kFullCircle);

  // No border on legends
  sbnStyle->SetLegendBorderSize(0);

  // Axis titles
  sbnStyle->SetTitleSize(.055, "xyz");
  sbnStyle->SetTitleOffset(.8, "xyz");
  // More space for y-axis to avoid clashing with big numbers
  sbnStyle->SetTitleOffset(.9, "y");
  // This applies the same settings to the overall plot title
  sbnStyle->SetTitleSize(.055, "");
  sbnStyle->SetTitleOffset(.8, "");
  // Axis labels (numbering)
  sbnStyle->SetLabelSize(.04, "xyz");
  sbnStyle->SetLabelOffset(.005, "xyz");

  // Prevent ROOT from occasionally automatically zero-suppressing
  sbnStyle->SetHistMinimumZero();

  // Thicker lines
  sbnStyle->SetHistLineWidth(2);
  sbnStyle->SetFrameLineWidth(2);
  sbnStyle->SetFuncWidth(2);

  // Set the number of tick marks to show
  sbnStyle->SetNdivisions(506, "xyz");

  // Set the tick mark style
  sbnStyle->SetPadTickX(1);
  sbnStyle->SetPadTickY(1);

  // Fonts
  const int kSbnFont = 42;
  sbnStyle->SetStatFont(kSbnFont);
  sbnStyle->SetLabelFont(kSbnFont, "xyz");
  sbnStyle->SetTitleFont(kSbnFont, "xyz");
  sbnStyle->SetTitleFont(kSbnFont, ""); // Apply same setting to plot titles
  sbnStyle->SetTextFont(kSbnFont);
  sbnStyle->SetLegendFont(kSbnFont);

  // Get moodier colours for colz
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  sbnStyle->SetNumberContours(NCont);

  gROOT->SetStyle("sbnStyle");

  // Uncomment this line if you want to force all plots loaded from files
  // to use this same style
  //gROOT->ForceStyle();
}

// Put a "SBN Preliminary" tag in the corner
void Preliminary()
{
  TLatex* prelim = new TLatex(.9, .95, "Sbn Preliminary");
  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

// Put a "SBN Preliminary" tag on the right
void PreliminarySide()
{
  TLatex* prelim = new TLatex(.93, .9, "SBN Preliminary");
  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAngle(270);
  prelim->SetTextAlign(12);
  prelim->Draw();
}

// Put a "SBN Simulation" tag in the corner
void Simulation()
{
  TLatex* prelim = new TLatex(.9, .95, "SBN Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

// Put a "SBN Simulation" tag on the right
void SimulationSide()
{
  TLatex* prelim = new TLatex(.93, .9, "SBN Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAngle(270);
  prelim->SetTextAlign(12);
  prelim->Draw();
}

void CenterTitles(TH1* histo)
{
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  histo->GetZaxis()->CenterTitle();  
}

#endif
