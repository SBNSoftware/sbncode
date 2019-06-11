#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TColor.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

using namespace selection;

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/breakdown/";
  std::string plots_location  = "../Output_Selection_Tool/plots/breakdown/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap numu_signal_map  = GeneralAnalysisHelper::GetNuMuTopologyMap();
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  /*
   *
   *        STATISTICAL BREAKDOWN OF EVENTS
   *        
   *        Total
   *
   *        Reconstructed vertex contained
   *
   *        All reconstructed tracks contained
   *
   *        For each of:
   *          True vertex contained 
   *          True vertex not contained
   *
   *            For each of: 
   *              Selected
   *              Signal
   *              True
   *              
   *                nu mu
   *                nu mu CC
   *                nu mu CC 0pi
   *                nu mu CC 1pi
   *
   *            Percentage of true particles with less than 10 hits associated
   *              Protons
   *              Muons
   *              Pions
   *        
   */

  unsigned int total                         = 0;
  unsigned int contained_vertex              = 0;
  unsigned int contained_tracks              = 0;
  unsigned int min_one                       = 0;
  unsigned int true_proton_less_than_10_hits = 0;   
  unsigned int true_muon_less_than_10_hits   = 0;   
  unsigned int true_pion_less_than_10_hits   = 0;
  unsigned int total_true_proton             = 0;
  unsigned int total_true_muon               = 0;
  unsigned int total_true_pion               = 0;
  unsigned int true_contained                = 0;
  unsigned int true_outside                  = 0;
  unsigned int sel_contained_numu            = 0;
  unsigned int sel_contained_numu_cc         = 0;
  unsigned int sel_contained_numu_cc_0pi     = 0;
  unsigned int sel_contained_numu_cc_1pi     = 0;
  unsigned int sel_outside_numu              = 0;
  unsigned int sel_outside_numu_cc           = 0;
  unsigned int sel_outside_numu_cc_0pi       = 0;
  unsigned int sel_outside_numu_cc_1pi       = 0;
  unsigned int sig_contained_numu            = 0;
  unsigned int sig_contained_numu_cc         = 0;
  unsigned int sig_contained_numu_cc_0pi     = 0;
  unsigned int sig_contained_numu_cc_1pi     = 0;
  unsigned int sig_outside_numu              = 0;
  unsigned int sig_outside_numu_cc           = 0;
  unsigned int sig_outside_numu_cc_0pi       = 0;
  unsigned int sig_outside_numu_cc_1pi       = 0;
  unsigned int tru_contained_numu            = 0;
  unsigned int tru_contained_numu_cc         = 0;
  unsigned int tru_contained_numu_cc_0pi     = 0;
  unsigned int tru_contained_numu_cc_1pi     = 0;
  unsigned int tru_outside_numu              = 0;
  unsigned int tru_outside_numu_cc           = 0;
  unsigned int tru_outside_numu_cc_0pi       = 0;
  unsigned int tru_outside_numu_cc_1pi       = 0;

  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 50;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    //if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    //name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    name = "/home/rhiannon/Samples/LocalSamples/analysis/200219_neutrino_only/selection/"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    
    //std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;
    EventSelectionTool::GetTimeLeft(start,total_files,i);

    /*
    TFile f(file_name);
    TTree *cut = (TTree*) f.Get("cut_tree");

    TBranch *b_total = cut->GetBranch("c_total");
    TBranch *b_contained = cut->GetBranch("c_contained");
    TBranch *b_tracks = cut->GetBranch("c_contained_tracks");
    TBranch *b_minone = cut->GetBranch("c_min_one");
    
    cut->GetEntry(0);

    total            += b_total->GetLeaf("c_total")->GetValue();
    contained_vertex += b_contained->GetLeaf("c_contained")->GetValue();
    contained_tracks += b_tracks->GetLeaf("c_contained_tracks")->GetValue();
    min_one          += b_minone->GetLeaf("c_min_one")->GetValue();
    */

  }
  std::cout << std::endl;
  /*
   *
   *          OF THE REMAINING EVENTS
   *
   *
   */

  for(const Event &e : events){

    ParticleList reco = e.GetRecoParticleList();
    ParticleList mc   = e.GetMCParticleList();

    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;

    if(e.IsSBNDTrueFiducial()){
      // Topology based contained and external (true) studies
      true_contained++;
      // Particle type hit-based studies
      for(const Particle &p : mc){
        if(p.GetPdgCode() == 13) total_true_muon++;
        if(p.GetPdgCode() == 211 || p.GetPdgCode() == -211) total_true_pion++;
        if(p.GetPdgCode() == 2212) total_true_proton++;
        if(p.GetNumberOfHits() <= 10){
          if(p.GetPdgCode() == 13) true_muon_less_than_10_hits++;
          if(p.GetPdgCode() == 211 || p.GetPdgCode() == -211) true_pion_less_than_10_hits++;
          if(p.GetPdgCode() == 2212) true_proton_less_than_10_hits++;
        }
      }
      // Selected and signal
      if(e.CheckRecoTopology(numu_signal_map)) {
        sel_contained_numu++; 
        if(e.CheckMCTopology(numu_signal_map)) sig_contained_numu++;
      }
      if(e.CheckRecoTopology(cc_signal_map)) {
        sel_contained_numu_cc++; 
        if(e.CheckMCTopology(cc_signal_map)) sig_contained_numu_cc++;
      }
      if(e.CheckRecoTopology(cc0pi_signal_map)) {
        sel_contained_numu_cc_0pi++; 
        if(e.CheckMCTopology(cc0pi_signal_map)) sig_contained_numu_cc_0pi++;
      }
      if(e.CheckRecoTopology(cc1pi_signal_map)) {
        sel_contained_numu_cc_1pi++;
        if(e.CheckMCTopology(cc1pi_signal_map)) sig_contained_numu_cc_1pi++;
      }
      // True
      if(e.CheckMCTopology(numu_signal_map))  tru_contained_numu++; 
      if(e.CheckMCTopology(cc_signal_map))    tru_contained_numu_cc++; 
      if(e.CheckMCTopology(cc0pi_signal_map)) tru_contained_numu_cc_0pi++; 
      if(e.CheckMCTopology(cc1pi_signal_map)) tru_contained_numu_cc_1pi++;
    }
    else{
      true_outside++;
      // Selected and signal
      if(e.CheckRecoTopology(numu_signal_map)) {
        sel_outside_numu++; 
        if(e.CheckMCTopology(numu_signal_map)) sig_outside_numu++;
      }
      if(e.CheckRecoTopology(cc_signal_map)) {
        sel_outside_numu_cc++; 
        if(e.CheckMCTopology(cc_signal_map)) sig_outside_numu_cc++;
      }
      if(e.CheckRecoTopology(cc0pi_signal_map)) {
        sel_outside_numu_cc_0pi++; 
        if(e.CheckMCTopology(cc0pi_signal_map)) sig_outside_numu_cc_0pi++;
      }
      if(e.CheckRecoTopology(cc1pi_signal_map)) {
        sel_outside_numu_cc_1pi++;
        if(e.CheckMCTopology(cc1pi_signal_map)) sig_outside_numu_cc_1pi++;
      }
      // True
      if(e.CheckMCTopology(numu_signal_map))  tru_outside_numu++; 
      if(e.CheckMCTopology(cc_signal_map))    tru_outside_numu_cc++; 
      if(e.CheckMCTopology(cc0pi_signal_map)) tru_outside_numu_cc_0pi++; 
      if(e.CheckMCTopology(cc1pi_signal_map)) tru_outside_numu_cc_1pi++;
    }
  }

  /*
   *
   *      OUTPUTS
   *
   */

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"statistical_breakdown.txt");


  file << "===========================================================" << std::endl;
  file << " Reco vertex contained          : " << events.size() << std::endl;
  file << " Reco and true vertex contained : " << true_contained << std::endl; 
  file << " Reco contained, true outside   : " << true_outside << std::endl; 
  file << "===========================================================" << std::endl;
  file << " True contained " << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Truth " << std::endl;
  file << "     NuMu                       : " << tru_contained_numu << std::endl;
  file << "     NuMu CC Inc.               : " << tru_contained_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << tru_contained_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << tru_contained_numu_cc_1pi << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Selected " << std::endl;
  file << "     NuMu                       : " << sel_contained_numu << std::endl;
  file << "     NuMu CC Inc.               : " << sel_contained_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << sel_contained_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << sel_contained_numu_cc_1pi << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Signal " << std::endl;
  file << "     NuMu                       : " << sig_contained_numu << std::endl;
  file << "     NuMu CC Inc.               : " << sig_contained_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << sig_contained_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << sig_contained_numu_cc_1pi << std::endl;
  file << "===========================================================" << std::endl;
  file << " True outside " << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Truth " << std::endl;
  file << "     NuMu                       : " << tru_outside_numu << std::endl;
  file << "     NuMu CC Inc.               : " << tru_outside_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << tru_outside_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << tru_outside_numu_cc_1pi << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Selected " << std::endl;
  file << "     NuMu                       : " << sel_outside_numu << std::endl;
  file << "     NuMu CC Inc.               : " << sel_outside_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << sel_outside_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << sel_outside_numu_cc_1pi << std::endl;
  file << "-----------------------------------------------------------" << std::endl;
  file << "   Signal " << std::endl;
  file << "     NuMu                       : " << sig_outside_numu << std::endl;
  file << "     NuMu CC Inc.               : " << sig_outside_numu_cc << std::endl;
  file << "     NuMu CC 0Pi                : " << sig_outside_numu_cc_0pi << std::endl;
  file << "     NuMu CC 1Pi                : " << sig_outside_numu_cc_1pi << std::endl;
  file << "===========================================================" << std::endl;
  file << " Protons with < 10 hits         : " << 100*(true_proton_less_than_10_hits/double(total_true_proton)) << std::endl;
  file << " Muons with < 10 hits           : " << 100*(true_muon_less_than_10_hits/double(total_true_muon)) << std::endl;
  file << " Pions with < 10 hits           : " << 100*(true_pion_less_than_10_hits/double(total_true_pion)) << std::endl;
  file << "===========================================================" << std::endl;

  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
