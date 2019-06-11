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
  std::cout << " Start local time and date: " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total; ++i ){

    //if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    
    //std::cout << "Loaded file " << std::setw(4) << i << '\r' << flush;
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
 
  /*
   *
   *      COUNTERS
   *
   */
  unsigned int true_cc0pi     = 0;
  unsigned int signal_cc0pi   = 0;
  unsigned int selected_cc0pi = 0;

  unsigned int nc_true         = 0;
  unsigned int cc_true         = 0;
  unsigned int nc0pi_true      = 0;
  unsigned int nc1pi_true      = 0;
  unsigned int nc2pi_true      = 0;
  unsigned int cc1pi_true      = 0;
  unsigned int cc2pi_true      = 0;

  /*
   *
   *      TOPOLOGY MAPS
   *
   */
  
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc2pi_signal_map = GeneralAnalysisHelper::GetCC2PiTopologyMap();
  TopologyMap nc0pi_signal_map = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_signal_map = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nc2pi_signal_map = GeneralAnalysisHelper::GetNC2PiTopologyMap();

  for(const Event &e : events){
    // Only look at events with 1 escaping track
    if(!e.IsSBNDTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) != 1) continue;
    if(e.CheckRecoTopology(cc0pi_signal_map)) {
      selected_cc0pi++;
      if(e.CheckMCTopology(cc0pi_signal_map)) signal_cc0pi++;
      if(e.CheckMCTopology(nc_signal_map)) nc_true++;
      if(e.CheckMCTopology(cc_signal_map)) cc_true++;
      if(e.CheckMCTopology(nc0pi_signal_map)) nc0pi_true++;
      if(e.CheckMCTopology(nc1pi_signal_map)) nc1pi_true++;
      if(e.CheckMCTopology(nc2pi_signal_map)) nc2pi_true++;
      if(e.CheckMCTopology(cc1pi_signal_map)) cc1pi_true++;
      if(e.CheckMCTopology(cc2pi_signal_map)) cc2pi_true++;
    }
    if(e.CheckMCTopology(cc0pi_signal_map)) true_cc0pi++;
  }

  /*
   *
   *      OUTPUTS
   *
   */

  // Files to hold particle statistics
  ofstream file;
  file.open(stats_location+"cc0pi_performance.txt");


  file << "============================="        << std::endl;
  file << " True CC0Pi     : " << true_cc0pi     << std::endl;
  file << " Selected CC0Pi : " << selected_cc0pi << std::endl;
  file << " Signal CC0Pi   : " << signal_cc0pi   << std::endl;
  file << "-----------------------------"        << std::endl;
  file << " BG NC0Pi       : " << nc0pi_true     << std::endl;
  file << " BG NC1Pi       : " << nc1pi_true     << std::endl;
  file << " BG NC2Pi       : " << nc2pi_true     << std::endl;
  file << " BG CC1pi       : " << cc1pi_true     << std::endl;
  file << " BG CC2pi       : " << cc2pi_true     << std::endl;
  file << " BG CC Other    : " << cc_true - signal_cc0pi - cc1pi_true - cc2pi_true << std::endl;
  file << " BG NC Other    : " << nc_true - nc0pi_true - nc1pi_true - nc2pi_true << std::endl;
  file << " CC + NC Inc.   : " << cc_true + nc_true << std::endl;
  file << "============================="        << std::endl;

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
