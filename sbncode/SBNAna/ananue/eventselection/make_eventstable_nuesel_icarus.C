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

void make_eventstable_nuesel_icarus(bool crtveto = true)
{

  std::string inDir  = "/icarus/data/users/dmendez/SBNAna/ananue/files/Jan2021/";
  std::string inFile_nue  = inDir + "nue_spectra_slice.root";
  std::string inFile_nus  = inDir + "nucosmics_spectra_slice.root";
  std::string inFile_cos  = inDir + "cosmics_spectra_slice.root";
  std::string outDir  = "/icarus/data/users/dmendez/SBNAna/ananue/tables/Jan2021/";

  std::string pot_tag = "6.6 #times 10^{20} POT";
  double POT = 6.6E20;

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

  const unsigned int kNSel = sels.size();
  
  std::vector<std::string> names_cut;
  std::vector<float> events_nue;
  std::vector<float> events_numu;
  std::vector<float> events_nc;
  std::vector<float> events_cos;
  std::vector<float> events_other;
  std::vector<float> events_totbkg;
  std::vector<float> events_total;
  std::vector<float> eff;
  std::vector<float> pur;

  unsigned int endLoop = kNSel;
  if(crtveto) endLoop = kNSel+2; // add one for crtveto and one for the N-1 crtveto
  for(unsigned int iSel = 0; iSel < endLoop; iSel++){ 

    // std::string mysuffixall;
    std::string mysuffix;
    std::string cutname;
    if(iSel<kNSel){
      // mysuffixall = "nocut_count";
      mysuffix      = sels[iSel].suffix + "_count";
      cutname       = sels[iSel].label;
    }
    if(crtveto){
      if(iSel==(limitN1-1)){ // everything id, and crtveto
        mysuffix = sels[iSel].suffix + "_veto_count"; // this should equal everything_veto_cout
        cutname  = "Everything (incl. CRT veto)";
      }
      if(iSel==(endLoop-2)){
        mysuffix = "nocut_veto_count"; // get the crtveto result only
        cutname  = "CRT veto";
      }
      if(iSel==(endLoop-1)){
        mysuffix = "everything_count"; // get the N-1 crtveto result, which would be everything
        cutname  = "N1 CRT veto";
      }
    }
    names_cut.push_back(cutname);
    
    std::cout << "cutname: " << cutname << ", mysuffix: " << mysuffix << std::endl;
    Spectrum *spec_allnue = LoadFromFile<Spectrum>(inFile_nue, "nue_nocut_count").release();
    Spectrum *spec_nue    = LoadFromFile<Spectrum>(inFile_nue, "nue_"+mysuffix).release();
    Spectrum *spec_nuenus = LoadFromFile<Spectrum>(inFile_nus, "nue_"+mysuffix).release();
    Spectrum *spec_numu   = LoadFromFile<Spectrum>(inFile_nus, "numu_"+mysuffix).release();
    Spectrum *spec_nc     = LoadFromFile<Spectrum>(inFile_nus, "nunc_"+mysuffix).release();
    Spectrum *spec_total  = LoadFromFile<Spectrum>(inFile_nus, "total_"+mysuffix).release();
    Spectrum *spec_cosnus = LoadFromFile<Spectrum>(inFile_nus, "cosmic_"+mysuffix).release();
    Spectrum *spec_cos    = LoadFromFile<Spectrum>(inFile_cos, "cosmic_"+mysuffix).release();
    TH1* hallnue = spec_allnue->ToTH1(POT);     // from nue
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

    float iallnue = hallnue->Integral();
    float inue    = hnue->Integral();
    float inumu   = hnumu->Integral();
    float inc     = hnc->Integral();
    float icos    = htotcos->Integral();
    float iother  = hother->Integral();
    float itotbkg = htotbkg->Integral();
    float itotal  = htotal->Integral();
    float ieff    = inue / iallnue;
    float ipur    = inue / (inue + itotbkg);      

    events_nue.push_back(inue);
    events_numu.push_back(inumu);
    events_nc.push_back(inc);
    events_cos.push_back(icos);
    events_other.push_back(iother);
    events_totbkg.push_back(itotbkg);
    events_total.push_back(itotal);
    if(iSel==0){
      eff.push_back(0.);
      pur.push_back(0.);      
    }
    else{
      eff.push_back(ieff);
      pur.push_back(ipur);
    }
    
  } // iSel simple

  // print single cuts only table
  printTableHeader();
  for(unsigned int iSel=0; iSel<limitN1; iSel++){
    if(iSel==(limitN1-1) && crtveto){
      printEventsLine(names_cut[endLoop-2], events_nue[endLoop-2], events_numu[endLoop-2], events_nc[endLoop-2], events_cos[endLoop-2], events_other[endLoop-2], eff[endLoop-2], pur[endLoop-2]);
      printEventsLine(names_cut[iSel], events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel], eff[iSel], pur[iSel]);
    }
    else
      printEventsLine(names_cut[iSel], events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel], eff[iSel], pur[iSel]);
  }
  printTableFooter();

  // print n-1 cuts table
  printTableHeader();
  printEventsLine(names_cut[0], events_nue[0], events_numu[0], events_nc[0], events_cos[0], events_other[0], eff[0], pur[0]); // print nocut for reference
  printEventsLine(names_cut[limitN1-1], events_nue[limitN1-1], events_numu[limitN1-1], events_nc[limitN1-1], events_cos[limitN1-1], events_other[limitN1-1], eff[limitN1-1], eff[limitN1-1]); // print everything for reference
  for(unsigned int iSel=limitN1; iSel<kNSel; iSel++){
    printEventsLine(names_cut[iSel], events_nue[iSel], events_numu[iSel], events_nc[iSel], events_cos[iSel], events_other[iSel], eff[iSel], pur[iSel]);
  }
  if(crtveto)
    printEventsLine(names_cut[endLoop-1], events_nue[endLoop-1], events_numu[endLoop-1], events_nc[endLoop-1], events_cos[endLoop-1], events_other[endLoop-1], eff[endLoop-1], pur[endLoop-1]);
  printTableFooter();

  names_cut.clear();
  events_nue.clear();
  events_numu.clear();
  events_nc.clear();
  events_cos.clear();
  events_other.clear();
  events_totbkg.clear();
  events_total.clear();

} // end make_eventstable
