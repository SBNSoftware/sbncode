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
 
  // Output file location
  // I recommend making a directory close to the Selection_Tool directory
  // to store any output files you might produce
  // For instance:
  std::string stats_location = "../Output_Selection_Tool/statistics/";
  std::string plots_location = "../Output_Selection_Tool/plots/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
 
  // So that we can estimate how much time is left while we fill the event list
  int start = static_cast<int>(time(NULL));

  // Load the events into the event list
  // ===========================
  //            FYI
  //  Running over all 500 files 
  //  takes ~15 minutes with a 
  //  good connection
  //
  //  I usually test on ~50
  // ===========================
  unsigned int total_files = 500;
  for( unsigned int i = 0; i < total_files; ++i ){
    // These are bad directories in this case, so continue past them
    //if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7) continue;

    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/pnfs/sbnd/persistent/users/rsjones/analysis_sample/120918_ana_files/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    //Load the event list with the contents of the current file
    EventSelectionTool::LoadEventList(file_name, events, i);
    
    // See how much time is left in the loading
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;

  // TopologyMap is a pre-defined map which holds a topology
  // See include/Event.h for the definition
  TopologyMap cc0pi_topology = GeneralAnalysisHelper::GetCC0PiTopologyMap();

  // COUNTERS
  // For a simple check of 
  //    How many cc0pi events we see in truth
  //    How many we reconstruct
  //    How many cc0pi reconstructed events are also true cc0pi events
  //
  unsigned int number_of_true_cc0pi_events                   = 0;
  unsigned int number_of_reconstructed_cc0pi_events          = 0;
  unsigned int number_of_true_and_reconstructed_cc0pi_events = 0;

  // Loop over all the events we have loaded
  for(const Event &e : events){

    // Make sure the true vertex is contained and that only 1 track escapes
    if(!e.IsSBNDTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) != 1) continue;
    
    // Find out if the current event was a true cc0pi final state
    if(e.CheckMCTopology(cc0pi_topology)) {
      number_of_true_cc0pi_events++;
      // Now check if it was also a reconstructed cc0pi final state
      if(e.CheckRecoTopology(cc0pi_topology))
        number_of_true_and_reconstructed_cc0pi_events++;
    }
    // Check if it was a reconstructed cc0pi final state
    if(e.CheckRecoTopology(cc0pi_topology))
      number_of_reconstructed_cc0pi_events++;

    // Get the list of reconstructed particle from the current event
    ParticleList reconstructed_particles = e.GetRecoParticleList();
    
    // Loop over the particles
    for(const Particle &p : reconstructed_particles){
      // Get the particle's PDG code
      int pdg = p.GetPdgCode();
    }
  }

  // File to hold particle statistics
  ofstream file;
  file.open(stats_location+"example_cc0pi_efficiency_and_purity.txt");

  file << "================================================================" << std::endl;
  file << " CC 0Pi  efficiency : " << number_of_true_and_reconstructed_cc0pi_events/double(number_of_true_cc0pi_events)          << std::endl; 
  file << " CC 0Pi  purity     : " << number_of_true_and_reconstructed_cc0pi_events/double(number_of_reconstructed_cc0pi_events) << std::endl; 
  file << "================================================================" << std::endl;

  file.close();

  return 0;

} // MainTest()
