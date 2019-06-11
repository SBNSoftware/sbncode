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

using namespace selection;

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << endl;
 
  // Output file location
  std::string stats_location = "../Output_Selection_Tool/statistics/mcp0_9/";
  std::string plots_location  = "../Output_Selection_Tool/plots/mcp0_9";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
 
  int start = static_cast<int>(time(NULL));
  unsigned int total = 991;

  // Load the events into the event list
  for( unsigned int i = 0; i < total; ++i ){

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
    EventSelectionTool::GetTimeLeft(start,total,i);
  }

  std::cout << std::endl;
  
  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Files to hold particle statistics
  ofstream all_file;
  all_file.open(stats_location+"smalldistancecut_particle_stats.txt");

  ofstream mis_id_file;
  mis_id_file.open(stats_location+"smalldistancecut_mis_identification_stats.txt");

  GeneralAnalysisHelper::FillGeneralParticleStatisticsFile(events, all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, nc_signal_map, "NC Inclusive",  all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc_signal_map, "CC Inclusive",  all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc0pi_signal_map, "CC 0 Pi",    all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, cc1pi_signal_map, "CC 1 Pi+/-", all_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleStatisticsFile(events, ccpi0_signal_map, "CC Pi0",     all_file);

  GeneralAnalysisHelper::FillGeneralParticleMisIdStatisticsFile(events, mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, nc_signal_map, "NC Inclusive",  mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc_signal_map, "CC Inclusive",  mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc0pi_signal_map, "CC 0 Pi",    mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, cc1pi_signal_map, "CC 1 Pi+/-", mis_id_file);
  GeneralAnalysisHelper::FillTopologyBasedParticleMisIdStatisticsFile(events, ccpi0_signal_map, "CC Pi0",     mis_id_file);

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
