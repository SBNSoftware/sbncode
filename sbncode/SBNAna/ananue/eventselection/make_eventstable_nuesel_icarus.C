///////////////////////////////////////////////////////////////////////
// Author: Diana Patricia Mendez                                     //
// Contact: dmendezme@bnl.gov                                        //
// Last edited: January 15 2021                                      //   
//                                                                   // 
// Makes summary event count tables from the spectra produced by     //
// make_spectra.C                                                    //
// Setting N1table to true will make a table of the N1 cuts only     //
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

void plot_eventstable_nuesel_icarus(bool N1table=false)
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
  const int limitN1 = 13; // N1 cuts start at 13 position in the sels object
  
  // std::vector<float> rem_vect;
  std::vector<float> events_nue;
  std::vector<float> events_numu;
  std::vector<float> events_nc;
  std::vector<float> events_cos;
  std::vector<float> events_other;
  std::vector<float> events_totbkg;
  std::vector<float> events_total;

  // I want to make a plot for each var
  bool print_int = true;
  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){

    std::string thiscornertag = sels[iSel].suffix;
    std::string mysuffix    = sels[iSel].suffix + "_count";
    if(crtveto){
      mysuffix = mysuffix+"_veto";
    }
    Spectrum *spec_nue   = LoadFromFile<Spectrum>(inFile, "nue_"+mysuffix).release();
    Spectrum *spec_numu  = LoadFromFile<Spectrum>(inFile, "numu_"+mysuffix).release();
    Spectrum *spec_nc    = LoadFromFile<Spectrum>(inFile, "nunc_"+mysuffix).release();
    Spectrum *spec_cos   = LoadFromFile<Spectrum>(inFile, "cosmic_"+mysuffix).release();
    Spectrum *spec_total = LoadFromFile<Spectrum>(inFile, "total_"+mysuffix).release();
    TH1* hnue    = spec_nue->ToTH1(POT);
    TH1* hnumu   = spec_numu->ToTH1(POT);
    TH1* hnc     = spec_nc->ToTH1(POT);
    TH1* htotcos = spec_cos->ToTH1(POT);
    TH1* htotal  = spec_total->ToTH1(POT);
    TH1* htotbkg  = (TH1*)htotal->Clone(); // total bkg
    htotbkg->Add(hnue,-1);
    TH1* hother = (TH1*)htotal->Clone(); // other bkg that's not numu or nc
    hother->Add(hnue, -1);
    hother->Add(hnumu, -1);
    hother->Add(hnc, -1);
    hother->Add(htotcos, -1);

/* // when separatedly adding in time cosmics
      Spectrum *spec_allsig = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffixall).release(); // needed for purity
      Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffix).release();
      Spectrum *spec_nuenus = LoadFromFile<Spectrum>(inFile_numu, "nue_"+mysuffix).release();
      Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile_numu, "numu_"+mysuffix).release();
      Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile_numu, "nunc_"+mysuffix).release();
      Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile_cos,  "cosmics_"+mysuffix).release();
      Spectrum *spec_cosnus = LoadFromFile<Spectrum>(inFile_numu,  "cosmics_"+mysuffix).release();
      TH1* hallsig = spec_allsig->ToTH1(POT);  // from nue
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
    float icos    = htotcos->Integral();
    float iother  = hother->Integral();
    float itotbkg = htotbkg->Integral();
    float itotal  = htotal->Integral();

    events_nue.push_back(inue);
    events_numu.push_back(inumu);
    events_nc.push_back(inc);
    events_cos.push_back(icos);
    events_other.push_back(iother);
    events_totbkg.push_back(itotbkg);
    events_total.push_back(itotal);
      
  } // iSel simple

  int initLoop = 0;
  int endLoop = limitN1;
  if(N1table){
    initLoop = limitN1;
    endLoop = kNSel;
  }

  printTableHeader();
  if(N1table){
    // Always print the no cut result for reference
    printEventsLine(cutname, events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel]);
  }
  for(int iSel = initLoop; iSel < endLoop; ++iSel){
    std::string cutname = sels[iSel].label;
    printEventsLine(cutname, events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel]);
  }
  printTableFooter();

} // end make_eventstable
