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

#include "../../../CAFAna/rootlogon.C"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "helper_nuesel_icarus.h"
#include "tools_nuesel_icarus.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"


using namespace ana;

void plot_spectra_nuesel_icarus(
  std::string input = "nucosmics",
  bool combo = true,
  bool logscale = true,
  bool effpur  = false,
  bool crtveto = false,
  bool crtvars = false)
{

  std::string inDir  = "/icarus/data/users/dmendez/SBNAna/ananue/files/Jan2021/";
  std::string inFile = inDir + input + "_spectra_hadded1_slice.root";
  std::string inFile_nue  = inDir + "nue_spectra_hadded1_slice.root";
  std::string inFile_nus  = inDir + "nucosmics_spectra_hadded1_slice.root";
  std::string inFile_cos  = inDir + "cosmics_spectra_hadded1_slice.root";
  std::string outDir  = "/icarus/data/users/dmendez/SBNAna/ananue/plots/Jan2021/";

  std::string ext_tag = (logscale ? "_logscale" : "");
  std::string pot_tag = "6.6 #times 10^{20} POT";
  double POT = 6.6E20;
  double Livetime = 1;

  // select all slices or a single slice per spill
  std::vector<PlotDef> plots = plots_slice;
  std::vector<SelDef> types  = types_slice;
  std::vector<SelDef> sels   = sels_slice;
  // // Not automatically working as these are different structures
  // // Move from using these structures or make class
  // // or remember to change by hand ...
  // if(selspill){
  //   plots = plots_spill;
  //   types = types_spill;
  //   sels  = sels_spill;
  // }

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  // I want to make a plot for each var and cut
  for(unsigned int iVar = 0; iVar < kNVar; ++iVar){
    for(unsigned int iSel = 0; iSel < kNSel; ++iSel){

      std::string thiscornertag = sels[iSel].suffix;
      std::string mysuffix    = sels[iSel].suffix + "_" + plots[iVar].suffix;
      std::string mysuffixall = "nocut_" + plots[iVar].suffix;
      if(crtveto){
        mysuffix = sels[iSel].suffix + "_veto_" + plots[iVar].suffix;
        mysuffixall = "nocut_veto_" + plots[iVar].suffix;
      }

      TH1D* hallnue;
      TH1D* hnue;
      TH1D* hnuenus;
      TH1D* hnumu;
      TH1D* hnc;
      TH1D* hcosnus;
      TH1D* htotal;
      TH1D* hother;
      TH1D* hcos;
      TH1D* htotcos;
      TH1D* htotbkg;
      if(combo){ // make plots with all the samples
        Spectrum *spec_allnue = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffixall).release(); // needed for purity
        Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffix).release();
        Spectrum *spec_nuenus = LoadFromFile<Spectrum>(inFile_nus, "nue_"+mysuffix).release();
        Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile_nus, "numu_"+mysuffix).release();
        Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile_nus, "nunc_"+mysuffix).release();
        Spectrum *spec_total  = LoadFromFile<Spectrum>(inFile_nus, "total_"+mysuffix).release();
        Spectrum *spec_cosnus = LoadFromFile<Spectrum>(inFile_nus, "cosmic_"+mysuffix).release();
        Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile_cos, "cosmic_"+mysuffix).release();
        hallnue = spec_allnue->ToTH1(POT);  // from nue
        hnue    = spec_nue->ToTH1(POT);     // from nue
        hnuenus = spec_nuenus->ToTH1(POT);  // from nus+cosmics
        hnumu   = spec_numu->ToTH1(POT);    // from nus+cosmics
        hnc     = spec_nc->ToTH1(POT);      // from nus+cosmics
        hcosnus = spec_cosnus->ToTH1(POT);  // from nus+cosmics
        htotal  = spec_total->ToTH1(POT);   // from nus+cosmics
        hother  = (TH1D*)htotal->Clone();    // from nus+cosmics
        hcos    = spec_cos->ToTH1(POT);    // from cosmics only
        htotcos = (TH1D*)hcos->Clone(); // from cosmics only
        htotbkg = (TH1D*)htotal->Clone();
        hother->Add(hcosnus,-1);
        hother->Add(hnuenus,-1);
        hother->Add(hnumu,-1);
        hother->Add(hnc,-1);
        htotcos->Add(hcosnus,1); // add out of time cosmics
        htotbkg->Add(hnuenus,-1);
        htotbkg->Add(hcos,1);
      }
      else{
        Spectrum *spec_allnue = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffixall).release(); // needed for purity
        Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffix).release();
        Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile, "numu_"+mysuffix).release();
        Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile, "nunc_"+mysuffix).release();
        Spectrum *spec_total  = LoadFromFile<Spectrum>(inFile, "total_"+mysuffix).release();
        Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile, "cosmic_"+mysuffix).release();
        hallnue = spec_allnue->ToTH1(POT);  // from nue
        hnue    = spec_nue->ToTH1(POT);     // from nue
        hnumu   = spec_numu->ToTH1(POT);    // from nus+cosmics
        hnc     = spec_nc->ToTH1(POT);      // from nus+cosmics
        hcos    = spec_cos->ToTH1(POT);    // from cosmics only
        htotal  = spec_total->ToTH1(POT);   // from nus+cosmics
        hother  = (TH1D*)htotal->Clone();    // from nus+cosmics
        htotcos = (TH1D*)hcos->Clone(); // from cosmics only
        htotbkg = (TH1D*)htotal->Clone();
        hother->Add(hcos,-1);
        hother->Add(hnue,-1);
        hother->Add(hnumu,-1);
        hother->Add(hnc,-1);
        htotbkg->Add(hnue,-1);
      }

      // Make efficiency and purity graphs
      TGraph* gEff = SelEFForPURvsX(hnue, htotbkg, hallnue, true);
      TGraph* gPur = SelEFForPURvsX(hnue, htotbkg, hallnue, false);
      gEff->SetTitle("");
      gPur->SetTitle("");
      TH2* axEffPur = new TH2F("", ";GeV;", hnue->GetSize() - 2, hnue->GetXaxis()->GetXmin(), hnue->GetXaxis()->GetXmax(), 30, 0.1, 1.1);
      axEffPur->GetXaxis()->SetTitle(gPur->GetXaxis()->GetTitle());
      axEffPur->GetXaxis()->CenterTitle();
      axEffPur->GetYaxis()->SetTitle("");
      gEff->SetLineWidth(3); gEff->SetLineColor(color_eff);
      gPur->SetLineWidth(3); gPur->SetLineColor(color_pur);

      float iallnue = hallnue->Integral();
      float inue    = hnue->Integral();
      float inumu   = hnumu->Integral();
      float inc     = hnc->Integral();
      float itoscos = htotcos->Integral();
      float iother  = hother->Integral();
      float itotbkg = htotbkg->Integral();

      float efficiency = inue / iallnue;
      float purity     = inue / (inue + itotbkg);      
      float pnue = 100. * inue / (inue + itotbkg);
      float pbkg = 100. * itotbkg / (inue + itotbkg);
      // if(iSel==0){
      //   count_sig=inue;
      //   count_bkg=itotbkg;
      // }

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
      if(logscale){
        gPad->SetLogy();
        hother->GetYaxis()->SetRangeUser(1.,1.3 * GetHistMax({hnue,hnumu,hnc,htotcos,hother}));
      }
      else{
        hother->GetYaxis()->SetRangeUser(0.,1.3 * GetHistMax({hnue,hnumu,hnc,htotcos,hother}));        
      }
      hother->GetYaxis()->SetTitle("Slices");
      hother->Draw("hist same");
      htotcos->Draw("hist same");
      hnc->Draw("hist same");
      hnumu->Draw("hist same");
      hnue->Draw("hist same");
      hnue->Draw("hist same");
      DrawSigBkgIntText(hnue, htotbkg, 0.04);
      DrawComponentsLegend(hnue, hnumu, hnc, htotcos, hother);
      Simulation(true);
      CornerLabel(thiscornertag.c_str());

      // cEvents->Modified();
      cEvents->Update();
      std::string eventsPlotName = outDir+mysuffix+"_eventsel_"+(combo ? "combo" : input)+ ext_tag + ".pdf";
      cEvents->Print(eventsPlotName.c_str());
      if(cEvents) cEvents->Close();


      if(iSel!=0 && effpur){
        TCanvas *cEventsEffPur;
        TPad *padEvents, *padEffPur;
        SplitCanvas2(cEventsEffPur,  padEvents, padEffPur);
        padEvents->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogy(0);
        hother->GetXaxis()->SetLabelSize(0); // remove labels
        hother->GetYaxis()->SetRangeUser(0.0,1.3 * GetHistMax({hnue,hnumu,hnc,htotcos,hother}));
        hother->Draw("hist same");
        htotcos->Draw("hist same");
        hnc->Draw("hist same");
        hnumu->Draw("hist same");
        hnue->Draw("hist same");
        hnue->Draw("hist same");
        DrawSigBkgIntText(hnue, htotbkg, 0.03);
        DrawComponentsLegend(hnue, hnumu, hnc, htotcos, hother);
        Simulation(true);
        // CornerLabel(thiscornertag.c_str());
        cEventsEffPur->Update();

        padEffPur->cd();
        gPad->SetTickx(1);
        gPad->SetTicky(1);
        gPad->SetLogy(0);
        axEffPur->Draw();
        gEff->Draw("same L");
        gPur->Draw("same L");
        // DrawEffPurLegend(gEff, "Efficiency", gPur, "Purity");
        DrawIntEffPurLegend(efficiency, gEff, "Efficiency", purity, gPur, "Purity");
        // cEventsEffPur->Modified();
        cEventsEffPur->Update();
        gEff->Draw("same L");
        gPur->Draw("same L");
        gPad->RedrawAxis();
        // std::string eventsEffPurPlotName = outDir+mysuffix+"_eventsel_"+(combo ? "combo" : input)+"__effpur.pdf";
        std::string eventsEffPurPlotName = outDir+mysuffix+"_eventsel_"+(combo ? "combo" : input)+ext_tag+"__effpur.pdf";
        cEventsEffPur->Print(eventsEffPurPlotName.c_str());
        if(cEventsEffPur) cEventsEffPur->Close();

      } // end if effpur

    } // iSel simple
  } // iVar

}
