#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include "../include/Plane.h"
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

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<int> &exceptions);

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/mcp0_9/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 991;
  double pot = 0.; 

  std::vector<int> exceptions;
  exceptions.clear();

  // Read in txt file of list of empty input directories
  std::fstream exception_file("/home/rhiannon/Samples/LocalSamples/analysis/mcp0.9_neutrino_with_subrun/selection/exceptions.txt");
  std::string s_exc;
  while (std::getline(exception_file, s_exc)) {
    int i_exc;
    std::istringstream ss_exc(s_exc);
    ss_exc >> i_exc;
    exceptions.push_back(i_exc); 
    ss_exc.str(std::string());
    s_exc.clear();
  }

  std::cout << " Skipping files ending in : " << std::endl;
  for(const int & ex : exceptions)
    std::cout << " - " << ex << " - ";
  std::cout << std::endl;

  LoadAllEvents(events, total_files, start, pot, exceptions);

  // COUNTERS
  unsigned int nue_true  = 0; 
  
  unsigned int ncinc_true  = 0; 
  unsigned int ncinc_sig   = 0; 
  unsigned int ncinc_sel   = 0; 
  
  unsigned int nc0pi_true  = 0; 
  unsigned int nc0pi_sig   = 0; 
  unsigned int nc0pi_sel   = 0; 
  
  unsigned int nc1pi_true  = 0; 
  unsigned int nc1pi_sig   = 0; 
  unsigned int nc1pi_sel   = 0; 
  
  unsigned int ccinc_cc0pi = 0;
  unsigned int ccinc_cc1pi = 0;
  unsigned int ccinc_ccoth = 0;
  unsigned int ccinc_nc0pi = 0;
  unsigned int ccinc_nc1pi = 0;
  unsigned int ccinc_ncoth = 0;
  unsigned int ccinc_nue   = 0;
  unsigned int ccinc_true  = 0; 
  unsigned int ccinc_sig   = 0; 
  unsigned int ccinc_sel   = 0; 
  
  unsigned int cc0pi_cc0pi = 0;
  unsigned int cc0pi_cc1pi = 0;
  unsigned int cc0pi_ccoth = 0;
  unsigned int cc0pi_nc0pi = 0;
  unsigned int cc0pi_nc1pi = 0;
  unsigned int cc0pi_ncoth = 0;
  unsigned int cc0pi_nue   = 0;
  unsigned int cc0pi_true  = 0; 
  unsigned int cc0pi_sig   = 0; 
  unsigned int cc0pi_sel   = 0; 
  
  unsigned int cc1pi_cc0pi = 0;
  unsigned int cc1pi_cc1pi = 0;
  unsigned int cc1pi_ccoth = 0;
  unsigned int cc1pi_nc0pi = 0;
  unsigned int cc1pi_nc1pi = 0;
  unsigned int cc1pi_ncoth = 0;
  unsigned int cc1pi_nue   = 0;
  unsigned int cc1pi_true  = 0; 
  unsigned int cc1pi_sig   = 0; 
  unsigned int cc1pi_sel   = 0; 

  unsigned int cc0pi2p_cc0pi = 0;
  unsigned int cc0pi2p_cc1pi = 0;
  unsigned int cc0pi2p_ccoth = 0;
  unsigned int cc0pi2p_nc0pi = 0;
  unsigned int cc0pi2p_nc1pi = 0;
  unsigned int cc0pi2p_ncoth = 0;
  unsigned int cc0pi2p_nue   = 0;
  unsigned int cc0pi2p_true  = 0; 
  unsigned int cc0pi2p_sig   = 0; 
  unsigned int cc0pi2p_sel   = 0; 
  
  unsigned int all_tracks_contained   = 0;
  unsigned int max_one_escaping_track = 0;

  TopologyMap nc_map      = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc_map      = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi_map   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_map   = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc0pi2p_map = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap nc0pi_map   = GeneralAnalysisHelper::GetNC0PiTopologyMap();
  TopologyMap nc1pi_map   = GeneralAnalysisHelper::GetNC1PiTopologyMap();
  TopologyMap nue_map     = GeneralAnalysisHelper::GetNuETopologyMap();

  std::vector< TopologyMap > maps({cc_map, cc0pi_map, cc1pi_map, cc0pi2p_map, nc_map, nc0pi_map, nc1pi_map, nue_map});  

  // First, ensure all tracks are contained
  for(const Event &e : events){
 
    // Check the true vertex is in the fiducial volume
    if(e.IsSBNDTrueFiducial()){
//      if(GeneralAnalysisHelper::NumberEscapingTracks(e) != 0) continue;
      if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;
      max_one_escaping_track++;
  
      if(e.CheckRecoTopology(maps[0])){
        if(e.CheckMCTopology(maps[1]))      ccinc_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) ccinc_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) ccinc_ccoth++;
        else if(e.CheckMCTopology(maps[5])) ccinc_nc0pi++;
        else if(e.CheckMCTopology(maps[6])) ccinc_nc1pi++;
        else if(e.CheckMCTopology(maps[4])) ccinc_ncoth++;
        else if(e.CheckMCTopology(maps[7])) ccinc_nue++;
      }
      if(e.CheckRecoTopology(maps[1])){
        if(e.CheckMCTopology(maps[1]))      cc0pi_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) cc0pi_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) cc0pi_ccoth++;
        else if(e.CheckMCTopology(maps[5])) cc0pi_nc0pi++;
        else if(e.CheckMCTopology(maps[6])) cc0pi_nc1pi++;
        else if(e.CheckMCTopology(maps[4])) cc0pi_ncoth++;
        else if(e.CheckMCTopology(maps[7])) cc0pi_nue++;
        if(e.CheckRecoTopology(maps[3])){
          if(e.CheckMCTopology(maps[1]))      cc0pi2p_cc0pi++;
          else if(e.CheckMCTopology(maps[2])) cc0pi2p_cc1pi++;
          else if(e.CheckMCTopology(maps[0])) cc0pi2p_ccoth++;
          else if(e.CheckMCTopology(maps[5])) cc0pi2p_nc0pi++;
          else if(e.CheckMCTopology(maps[6])) cc0pi2p_nc1pi++;
          else if(e.CheckMCTopology(maps[4])) cc0pi2p_ncoth++;
        else if(e.CheckMCTopology(maps[7]))   cc0pi2p_nue++;
        }
      }
      if(e.CheckRecoTopology(maps[2])){
        if(e.CheckMCTopology(maps[1]))      cc1pi_cc0pi++;
        else if(e.CheckMCTopology(maps[2])) cc1pi_cc1pi++;
        else if(e.CheckMCTopology(maps[0])) cc1pi_ccoth++;
        else if(e.CheckMCTopology(maps[5])) cc1pi_nc0pi++;
        else if(e.CheckMCTopology(maps[6])) cc1pi_nc1pi++;
        else if(e.CheckMCTopology(maps[4])) cc1pi_ncoth++;
        else if(e.CheckMCTopology(maps[7])) cc1pi_nue++;
      }
      // Overall efficiencies 
      if(e.CheckMCTopology(maps[0])){
        ccinc_true++;
        if(e.CheckRecoTopology(maps[0])) ccinc_sig++;
      }
      if(e.CheckMCTopology(maps[1])){
        cc0pi_true++;
        if(e.CheckRecoTopology(maps[1])) cc0pi_sig++;
      }
      if(e.CheckMCTopology(maps[2])) {
        cc1pi_true++;
        if(e.CheckRecoTopology(maps[2])) cc1pi_sig++;
      }
      if(e.CheckMCTopology(maps[3])){
        cc0pi2p_true++;
        if(e.CheckRecoTopology(maps[3])) cc0pi2p_sig++;
      }
      if(e.CheckMCTopology(maps[4])){
        ncinc_true++;
        if(e.CheckRecoTopology(maps[4])) ncinc_sig++;
      }
      if(e.CheckMCTopology(maps[5])){
        nc0pi_true++;
        if(e.CheckRecoTopology(maps[5])) nc0pi_sig++;
      }
      if(e.CheckMCTopology(maps[6])){
        nc1pi_true++;
        if(e.CheckRecoTopology(maps[6])) nc1pi_sig++;
      }
      if(e.CheckMCTopology(maps[7])) nue_true++;

      // Overall purities
      if(e.CheckRecoTopology(maps[0])) ccinc_sel++;
      if(e.CheckRecoTopology(maps[1])) cc0pi_sel++;
      if(e.CheckRecoTopology(maps[2])) cc1pi_sel++;
      if(e.CheckRecoTopology(maps[3])) cc0pi2p_sel++;
      if(e.CheckRecoTopology(maps[4])) ncinc_sel++;
      if(e.CheckRecoTopology(maps[5])) nc0pi_sel++;
      if(e.CheckRecoTopology(maps[6])) nc1pi_sel++;
    }
  }

  // Files to hold particle statistics
  ofstream file;
  
  file.open(stats_location+"topology_breakdown.txt");

  file << "===============================================================================" << std::endl;
  //file << " Total number of events with all tracks contained : " << all_tracks_contained << std::endl;
  file << " Total number of events with maximum one escaping track : " << max_one_escaping_track << std::endl;
  file << std::setw(12) << "True \\ Reco" << "||" <<  std::setw(10) << " CC Inc. " << std::setw(10) << " CC 0Pi " << std::setw(10) << " CC 0Pi 2P" << std::setw(10) << " CC 1Pi " << std::endl;
  file << std::setw(12) << " CC 0Pi "     << "||" << std::setw(10) << ccinc_cc0pi << std::setw(10) << cc0pi_cc0pi << std::setw(10) << cc0pi2p_cc0pi << std::setw(10) << cc1pi_cc0pi << std::endl;  
  file << std::setw(12) << " CC 1Pi "     << "||" << std::setw(10) << ccinc_cc1pi << std::setw(10) << cc0pi_cc1pi << std::setw(10) << cc0pi2p_cc1pi << std::setw(10) << cc1pi_cc1pi << std::endl;  
  file << std::setw(12) << " CC Other "   << "||" << std::setw(10) << ccinc_ccoth << std::setw(10) << cc0pi_ccoth << std::setw(10) << cc0pi2p_ccoth << std::setw(10) << cc1pi_ccoth << std::endl;  
  file << std::setw(12) << " NC 0Pi "     << "||" << std::setw(10) << ccinc_nc0pi << std::setw(10) << cc0pi_nc0pi << std::setw(10) << cc0pi2p_nc0pi << std::setw(10) << cc1pi_nc0pi << std::endl;  
  file << std::setw(12) << " NC 1Pi "     << "||" << std::setw(10) << ccinc_nc1pi << std::setw(10) << cc0pi_nc1pi << std::setw(10) << cc0pi2p_nc1pi << std::setw(10) << cc1pi_nc1pi << std::endl;  
  file << std::setw(12) << " NC Other "   << "||" << std::setw(10) << ccinc_ncoth << std::setw(10) << cc0pi_ncoth << std::setw(10) << cc0pi2p_ncoth << std::setw(10) << cc1pi_ncoth << std::endl;  
  file << std::setw(12) << " Nu E    "    << "||" << std::setw(10) << ccinc_nue   << std::setw(10) << cc0pi_nue   << std::setw(10) << cc0pi2p_nue   << std::setw(10) << cc1pi_nue << std::endl;  
  file << "===============================================================================" << std::endl;
  file << " CC Inc.    true       : " << ccinc_true << std::endl; 
  file << " CC 0Pi     true       : " << cc0pi_true << std::endl; 
  file << " CC 1Pi     true       : " << cc1pi_true << std::endl; 
  file << " NC Inc.    true       : " << ncinc_true << std::endl; 
  file << " NC 0Pi     true       : " << nc0pi_true << std::endl; 
  file << " NC 1Pi     true       : " << nc1pi_true << std::endl; 
  file << " Nu E       true       : " << nue_true   << std::endl; 
  file << "===============================================================================" << std::endl;
  file << " CC Inc.    efficiency : " << ccinc_sig/double(ccinc_true) << std::endl; 
  file << " CC 0Pi     efficiency : " << cc0pi_sig/double(cc0pi_true) << std::endl; 
  file << " CC 0Pi 2P  efficiency : " << cc0pi2p_sig/double(cc0pi2p_true) << std::endl; 
  file << " CC 1Pi     efficiency : " << cc1pi_sig/double(cc1pi_true) << std::endl; 
  file << " NC Inc.    efficiency : " << ncinc_sig/double(ncinc_true) << std::endl; 
  file << " NC 0Pi     efficiency : " << nc0pi_sig/double(nc0pi_true) << std::endl; 
  file << " NC 1Pi     efficiency : " << nc1pi_sig/double(nc1pi_true) << std::endl; 
  file << "===============================================================================" << std::endl;
  file << " CC Inc.   purity      : " << ccinc_sig/double(ccinc_sel)  << std::endl; 
  file << " CC 0Pi    purity      : " << cc0pi_sig/double(cc0pi_sel)  << std::endl; 
  file << " CC 0Pi 2P purity      : " << cc0pi2p_sig/double(cc0pi2p_sel)  << std::endl; 
  file << " CC 1Pi    purity      : " << cc1pi_sig/double(cc1pi_sel)  << std::endl; 
  file << " NC Inc.   purity      : " << ncinc_sig/double(ncinc_sel)  << std::endl; 
  file << " NC 0Pi    purity      : " << nc0pi_sig/double(nc0pi_sel)  << std::endl; 
  file << " NC 1Pi    purity      : " << nc1pi_sig/double(nc1pi_sel)  << std::endl; 
  file << "===============================================================================" << std::endl;

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<int> &exceptions) {
  double total_pot = 0;
  std::vector<int>::iterator it;
  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    it = std::find(exceptions.begin(), exceptions.end(),i);
    if(it != exceptions.end()) continue;
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
//    name = "/home/rhiannon/Samples/LocalSamples/analysis/test/output_file.root";
      name = "/home/rhiannon/Samples/LocalSamples/analysis/mcp0.9_neutrino_with_subrun/selection/"+std::to_string(i)+"/output_file.root";
//    name = "/home/rhiannon/Samples/LocalSamples/analysis/200219_neutrino_only/selection/"+std::to_string(i)+"/output_file.root";
//    name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    double temp_pot = 0.;
    EventSelectionTool::LoadEventList(file_name, events, i, temp_pot);
    EventSelectionTool::GetTimeLeft(start_time,total_files,i);
    total_pot += temp_pot;
  }
  std::cout << std::endl;
  pot = total_pot;
} // LoadAllEvents
