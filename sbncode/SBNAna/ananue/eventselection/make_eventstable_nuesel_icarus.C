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

void make_eventstable_nuesel_icarus(bool N1table=false,
  bool crtveto = false)
{

  std::string inDir  = "/icarus/data/users/dmendez/SBNAna/ananue/files/Jan2021/";
  std::string inFile_nue  = inDir + "nue_spectra_hadded1_slice.root";
  std::string inFile_nus  = inDir + "nucosmics_spectra_hadded1_slice.root";
  std::string inFile_cos  = inDir + "cosmics_spectra_hadded1_slice.root";
  std::string outDir  = "/icarus/data/users/dmendez/SBNAna/ananue/tables/Jan2021/";

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

    Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffix).release();
    Spectrum *spec_nuenus = LoadFromFile<Spectrum>(inFile_nus, "nue_"+mysuffix).release();
    Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile_nus, "numu_"+mysuffix).release();
    Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile_nus, "nunc_"+mysuffix).release();
    Spectrum *spec_total  = LoadFromFile<Spectrum>(inFile_nus, "total_"+mysuffix).release();
    Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile_cos, "cosmic_"+mysuffix).release();
    Spectrum *spec_cosnus = LoadFromFile<Spectrum>(inFile_nus, "cosmic_"+mysuffix).release();
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
    htotcos->Add(hcosnus,1); // add out of time cosmics
    TH1* htotbkg = (TH1*)htotal->Clone();
    htotbkg->Add(hnuenus,-1);
    htotbkg->Add(hcos,1);

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
    printEventsLine("No cut", events_nue[0], events_numu[0], events_nc[0], events_cos[0], events_other[0]);
  }
  for(int iSel = initLoop; iSel < endLoop; ++iSel){
    std::string cutname = sels[iSel].label;
    printEventsLine(cutname, events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel]);
  }
  printTableFooter();

} // end make_eventstable
