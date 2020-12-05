#pragma once

#include "/icarus/app/users/dmendez/develop/srcs/sbncode/sbncode/CAFAna/Core/rootlogon.C"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "helper_nuesel_icarus.h"
#include "plotting_tools.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"


using namespace ana;


void FillWithDimColor(TH1* h, bool usealpha=false, float dim=0.8)
{
  if ( usealpha ){
    h->SetFillColorAlpha(h->GetLineColor(),dim);
    return;
  }
  TColor *color = gROOT->GetColor(h->GetLineColor());
  float R,G,B,hR,hG,hB,hHue,hSat,hVal;
  color->GetRGB(hR,hG,hB);
  color->RGB2HSV(hR,hG,hB,hHue,hSat,hVal);
  color->HSV2RGB(hHue,dim*hSat,hVal,R,G,B);
  h->SetFillColor(color->GetColor(R,G,B));
}

void plot_spectra_nuesel_icarus(
  std::string input = "nucosmics",
  bool crtveto = false,
  bool crtvars = false)
{

  std::string inFile = input+"_spectra.root";
  // std::string innue  = "nuecosmics_spectra.root";
  // std::string innus  = "nucosmics_spectra.root";
  // std::string incos  = "cosmics_spectra.root";

  std::string pot_tag = "6.6 #times 10^{20} POT";
  double POT = 6.6E20;
  double Livetime = 1;

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // std::vector<float> rem_vect;
  std::vector<float> sig_vect;
  std::vector<float> bkg_vect;

  // I want to make a plot for each var
  bool print_int = true;
  for(unsigned int iVar = 0; iVar < kNVar; ++iVar){
    for(unsigned int iSel = 0; iSel < kNSel; ++iSel){

      std::string thiscornertag = sels[iSel].suffix;
      std::string mysuffix = sels[iSel].suffix + "_" + plots[iVar].suffix;
      if(crtveto) mysuffix = mysuffix+"_veto";
      Spectrum *spec_nue   = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffix).release();
      Spectrum *spec_numu  = LoadFromFile<Spectrum>(inFile, "numu_"+mysuffix).release();
      Spectrum *spec_nc    = LoadFromFile<Spectrum>(inFile, "nunc_"+mysuffix).release();
      Spectrum *spec_total = LoadFromFile<Spectrum>(inFile, "total_"+mysuffix).release();
      TH1* hnue   = spec_nue->ToTH1(POT);
      TH1* hnumu  = spec_numu->ToTH1(POT);
      TH1* hnc    = spec_nc->ToTH1(POT);
      TH1* htotal = spec_total->ToTH1(POT);
      TH1* htotbkg  = (TH1*)htotal->Clone(); // total bkg
      htotbkg->Add(hnue,-1);
      TH1* hother = (TH1*)htotal->Clone(); // other bkg that's not numu or nc
      hother->Add(hnue, -1);
      hother->Add(hnumu, -1);
      hother->Add(hnc, -1);

/* // when separatedly adding in time cosmics
      Spectrum *spec_nue   = LoadFromFile<Spectrum>(innue, "nue_"+mysuffix).release();
      Spectrum *spec_nuenus= LoadFromFile<Spectrum>(innumu, "nue_"+mysuffix).release();
      Spectrum *spec_numu  = LoadFromFile<Spectrum>(innumu, "numu_"+mysuffix).release();
      Spectrum *spec_nc    = LoadFromFile<Spectrum>(innumu, "nunc_"+mysuffix).release();
      Spectrum *spec_other = LoadFromFile<Spectrum>(innumu, "other_"+mysuffix).release();
      Spectrum *spec_cos   = LoadFromFile<Spectrum>(incos,  "cosmics_"+mysuffix).release();
      Spectrum *spec_cosnus= LoadFromFile<Spectrum>(innumu,  "cosmics_"+mysuffix).release();
      TH1* hnue   = spec_nue->ToTH1(POT);     // nue+cosmics
      TH1* hnuenus= spec_nuenus->ToTH1(POT);  // nus+cosmics
      TH1* hnumu  = spec_numu->ToTH1(POT);    // nus+cosmics
      TH1* hnc    = spec_nc->ToTH1(POT);      // nus+cosmics
      TH1* hother = spec_other->ToTH1(POT);   // nus+cosmics
      TH1* hcosnus= spec_cosnus->ToTH1(POT);  // nus+cosmics
      TH1* htotal = spec_total->ToTH1(POT);   // nus+cosmics
      TH1* hcos   = spec_cos->ToTH1(POT);     // cosmics only
      TH1* htotbkg  = (TH1*)hcos->Clone();
      htotbkg->Add(hnumu,1);
      htotbkg->Add(hnc,1);
      htotbkg->Add(hother,1);
*/

      float inue    = hnue->Integral();
      float inumu   = hnumu->Integral();
      float inc     = hnc->Integral();
      float iother  = hother->Integral();
      float itotbkg = htotbkg->Integral();

      sig_vect.push_back(inue);
      bkg_vect.push_back(itotbkg);
      
      float pnue = 100 * inue / (inue + itotbkg);
      float pbkg = 100 * itotbkg / (inue + itotbkg);
      PimpHist(hnue, color_nue, line_nue, 2);
      PimpHist(hnumu, color_numu, line_numu, 2);
      PimpHist(hnc, color_nc, line_nc, 2);
      PimpHist(hother, color_other, line_other, 3);
      CenterTitles(hnue);
      CenterTitles(hnumu);
      CenterTitles(hnc);
      CenterTitles(hother);

      TCanvas *c = new TCanvas(plots[iVar].suffix.c_str(),plots[iVar].suffix.c_str(), 700, 500);
      hnue->GetYaxis()->SetRangeUser(0.0,1.4 * GetHistMax({hnue,hnumu,hnc,hother}));
      hnue->GetYaxis()->SetTitle("Slices");
      if(iVar==0){ // fake stacking
        hnumu->Add(hnc);
        hnc->Add(hnumu);
        hother->Add(hnc);
        FillWithDimColor(hnue, 0);
        FillWithDimColor(hnumu, 0);
        FillWithDimColor(hnc, 0);
        FillWithDimColor(hother, 0);
      }
      hother->Draw("hist same");
      hnc->Draw("hist same");
      hnumu->Draw("hist same");
      hnue->Draw("hist same");
      TPaveText *pText1 = new TPaveText(0.15, 0.78, 0.30, 0.85, "brNDC");
      TText *text1 = (pText1->AddText(pot_tag.c_str()));
      text1->SetTextSize(0.04);
      pText1->SetBorderSize(0);
      pText1->SetFillStyle(0);
      pText1->Draw();
      TPaveText *pText2 = new TPaveText(0.15, 0.65, 0.30, 0.75, "brNDC");
      TText *text2 = pText2->AddText(Form("sig: %2.f = %2.f %%", hnue->Integral(), pnue));
      text2->SetTextAlign(11);
      text2->SetTextSize(0.04);
      TText *text3 = pText2->AddText(Form("bkg: %2.f = %2.f %%", htotbkg->Integral(), pbkg));
      text3->SetTextAlign(11);
      text3->SetTextSize(0.04);
      pText2->SetBorderSize(0);
      pText2->SetFillStyle(0);
      pText2->Draw();
      // TLegend *l = new TLegend(0.60, 0.65, 0.85, 0.85, NULL,"brNDC");
      TLegend *l = new TLegend(0.70, 0.70, 0.85, 0.85, NULL,"brNDC");
      l->SetFillStyle(0);
      l->SetTextSize(0.035);
      // l->SetHeader(pot_tag.c_str());
      // l->SetHeader("Integral");
      // l->AddEntry(hnue, Form("#nu_{e} CC: %.2f", hnue->Integral()), "l");
      // l->AddEntry(hnumu, Form("#nu_{#mu} CC: %.2f", hnumu->Integral()), "l");
      // l->AddEntry(hnc, Form("NC: %.2f", hnc->Integral()), "l");
      // l->AddEntry(hother, Form("Other bkg: %.2f", hother->Integral()), "l");
      l->AddEntry(hnue, "#nu_{e} CC", "l");
      l->AddEntry(hnumu, "#nu_{#mu} CC", "l");
      l->AddEntry(hnc, "NC", "l");
      l->AddEntry(hother, "Other bkg", "l");
      l->Draw("");
      Simulation(true);
      CornerLabel(thiscornertag.c_str());

      c->Modified();
      c->Update();
      c->Print(("plots/"+input+"/"+mysuffix+"_eventsel_"+input+".pdf").c_str());

      if(c) c->Close();

    } // iSel simple

    print_int = false;
  } // iVar

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    std::cout << sels[iSel].suffix << ", " << "," << sig_vect[iSel] << "," << bkg_vect[iSel] << std::endl;
  }

}
