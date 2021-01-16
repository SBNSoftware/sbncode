///////////////////////////////////////////////////////////////////////
// Author: Diana Patricia Mendez                                     //
// Contact: dmendezme@bnl.gov                                        //
// Last edited: January 15 2021                                      //   
//                                                                   // 
// Plots the spectra produced by make_spectra.C                      //
// Setting effpur to true will plot a split canvas, with the spectra //
// at the top and efficiency and purity at the bottom.               //
///////////////////////////////////////////////////////////////////////

#pragma once

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "helper_nuesel_icarus.h"
#include "tools_nuesel_icarus.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"


using namespace ana;

void plot_spectra_nuesel_icarus(
  std::string input = "nucosmics",
  bool crtveto = false,
  bool crtvars = false,
  bool effpur  = false)
{

  std::string inFile = input+"_spectra.root";
  // std::string inFile_nue  = "nue_spectra.root";
  // std::string inFile_nus  = "nucosmics_spectra.root";
  // std::string inFile_cos  = "cosmics_spectra.root";

  std::string pot_tag = "6.6 #times 10^{20} POT";
  double POT = 6.6E20;
  double Livetime = 1;

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // std::vector<float> rem_vect;
  std::vector<float> sig_vect;
  std::vector<float> bkg_vect;

  // I want to make a plot for each var and cut
  for(unsigned int iVar = 0; iVar < kNVar; ++iVar){
    for(unsigned int iSel = 0; iSel < kNSel; ++iSel){

      std::string thiscornertag = sels[iSel].suffix;
      std::string mysuffix    = sels[iSel].suffix + "_" + plots[iVar].suffix;
      std::string mysuffixall = "nocut_" + plots[iVar].suffix;
      if(crtveto){
        mysuffix = mysuffix+"_veto";
        mysuffixall = mysuffixall+"_veto";
      }
      Spectrum *spec_allnue = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffixall).release();
      Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffix).release();
      Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile, "numu_"+mysuffix).release();
      Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile, "nunc_"+mysuffix).release();
      Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile, "cosmic_"+mysuffix).release();
      Spectrum *spec_total  = LoadFromFile<Spectrum>(inFile, "total_"+mysuffix).release();
      TH1* hallnue = spec_allnue->ToTH1(POT);
      TH1* hnue    = spec_nue->ToTH1(POT);
      TH1* hnumu   = spec_numu->ToTH1(POT);
      TH1* hnc     = spec_nc->ToTH1(POT);
      TH1* htotcos = spec_cos->ToTH1(POT);
      TH1* htotal  = spec_total->ToTH1(POT);
      TH1* htotbkg = (TH1*)htotal->Clone(); // total bkg
      htotbkg->Add(hnue,-1);
      TH1* hother = (TH1*)htotal->Clone(); // other bkg that's not numu or nc
      hother->Add(hnue, -1);
      hother->Add(hnumu, -1);
      hother->Add(hnc, -1);
      hother->Add(htotcos, -1);

/* // when separatedly adding in time cosmics
      Spectrum *spec_allnue = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffixall).release(); // needed for purity
      Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffix).release();
      Spectrum *spec_nuenus = LoadFromFile<Spectrum>(inFile_numu, "nue_"+mysuffix).release();
      Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile_numu, "numu_"+mysuffix).release();
      Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile_numu, "nunc_"+mysuffix).release();
      Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile_cos,  "cosmics_"+mysuffix).release();
      Spectrum *spec_cosnus = LoadFromFile<Spectrum>(inFile_numu,  "cosmics_"+mysuffix).release();
      TH1* hallnue = spec_allsig->ToTH1(POT);  // from nue
      TH1* hnue    = spec_nue->ToTH1(POT);     // from nue
      TH1* hnuenus = spec_nuenus->ToTH1(POT);  // from nus+cosmics
      TH1* hnumu   = spec_numu->ToTH1(POT);    // from nus+cosmics
      TH1* hnc     = spec_nc->ToTH1(POT);      // from nus+cosmics
      TH1* hcosnus = spec_cosnus->ToTH1(POT);  // from nus+cosmics
      TH1* htotal  = spec_total->ToTH1(POT);   // from nus+cosmics
      TH1* hother  = (TH1*)htotal->Clone();    // from nus+cosmics
      hother->Add(hcosnus,-1);
      hother->Add(hnuenus,-1);
      hother->Add(hnumu,-1);
      hother->Add(hnc,-1);
      TH1* hcos= spec_cos->ToTH1(POT);    // from cosmics only
      TH1* htotcos = (TH1*)hcos->Clone(); // from cosmics only
      TH1* htotcos->Add(hcosnus,1);       // add out of time cosmics
      TH1* htotbkg = (TH1*)htotal->Clone();
      htotbkg->Add(hnuenus,-1);
      htotbkg->Add(hcos,1);
*/

      float inue    = hnue->Integral();
      float inumu   = hnumu->Integral();
      float inc     = hnc->Integral();
      float itoscos = htotcos->Integral();
      float iother  = hother->Integral();
      float itotbkg = htotbkg->Integral();

      sig_vect.push_back(inue);
      bkg_vect.pussh_back(itotbkg);
      
      float pnue = 100 * inue / (inue + itotbkg);
      float pbkg = 100 * itotbkg / (inue + itotbkg);
      PimpHist(hnue, color_nue, line_nue, 2);
      PimpHist(hnumu, color_numu, line_numu, 2);
      PimpHist(hnc, color_nc, line_nc, 2);
      PimpHist(htotcos, color_cos, line_cos, 2);
      PimpHist(hother, color_other, line_other, 3);
      CenterTitles(hnue);
      CenterTitles(hnumu);
      CenterTitles(hnc);
      CenterTitles(htotcos);
      CenterTitles(hother);

      TCanvas *cEvents = new TCanvas(plots[iVar].suffix.c_str(),plots[iVar].suffix.c_str(), 700, 500);
      hnue->GetYaxis()->SetRangeUser(0.0,1.4 * GetHistMax({hnue,hnumu,hnc,htotcos,hother}));
      hnue->GetYaxis()->SetTitle("Slices");
      if(iVar==0){ // fake stacking
        hnumu->Add(hnc);
        hnc->Add(hnumu);
        htotcos->Add(hnc);
        hother->Add(htotcos);
        FillWithDimColor(hnue, 0);
        FillWithDimColor(hnumu, 0);
        FillWithDimColor(hnc, 0);
        FillWithDimColor(htotcos, 0);
        FillWithDimColor(hother, 0);
      }
      hother->Draw("hist same");
      htotcos->Draw("hist same");
      hnc->Draw("hist same");
      hnumu->Draw("hist same");
      hnue->Draw("hist same");
      DrawSigBkgIntText(hnue, htotbkg);
      DrawComponentsLegend(hnue, hnumu, hnc, htotcos, hother);
      Simulation(true);
      CornerLabel(thiscornertag.c_str());

      cEvents->Modified();
      cEvents->Update();
      cEvents->Print(("plots/"+input+"/"+mysuffix+"_eventsel_"+input+".pdf").c_str());
      if(cEvents) cEvents->Close();


      if(iSel!=0 && effpur){
        TCanvas *c1;
        TPad *padEvents, *padEffPur;
        SplitCanvas2(c1, padEvents, padEffPur);
        padEvents->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        hother->Draw("hist same");
        htotcos->Draw("hist same");
        hnc->Draw("hist same");
        hnumu->Draw("hist same");
        hnue->Draw("hist same");
        DrawSigBkgIntText(hnue, htotbkg);
        DrawComponentsLegend(hnue, hnumu, hnc, htotcos, hother);
        Simulation(true);
        CornerLabel(thiscornertag.c_str());
        padEffPur->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        TGraph* gEff = SelEFForPURvsX(hnue, htotbkg, hallnue, true);
        TGraph* gPur = SelEFForPURvsX(hnue, htotbkg, hallnue, false);
        gEff->Draw("l same");
        gPur>Draw("l same");
        DrawEffPurLegend(gEff, "Efficiency", gPur, "Purity");
        cEventsEffPur->Modified();
        cEventsEffPur->Update();
        cEventsEffPur->Print(("plots/"+input+"/"+mysuffix+"_eventsel_"+input+"__effpur.pdf").c_str());
        if(cEventsEffPurr) cEventsEffPur->Close();

      } // end if effpur

    } // iSel simple
  } // iVar

}
