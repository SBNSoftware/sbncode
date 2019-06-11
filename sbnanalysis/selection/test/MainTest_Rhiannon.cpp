#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <numeric>
#include <time.h>
#include "TVector3.h"

using namespace selection;

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start: Local time and date:  " << asctime(timeinfo) << std::endl;
  std::cout << " 100 sample " << std::endl;

  std::string filename = "/home/jtv/Desktop/Event_Selection_Tool/analysis_trees_2.root";

  EventSelectionTool::EventList events;
  EventSelectionTool::LoadEventList(filename, events);

  // Counters
  unsigned int correctly_reconstructed = 0;
  unsigned int true_topology = 0;
  unsigned int reco_topology = 0;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    TopologyMap signal_map_all;
   
    std::vector< int > cc_0pi_mu;
    std::vector< int > cc_0pi_pi;
    
    cc_0pi_mu.push_back( 13 );
    cc_0pi_pi.push_back( 211 );
    cc_0pi_pi.push_back(-211 );
    cc_0pi_pi.push_back( 111 );
    
    signal_map_all.insert( std::make_pair( cc_0pi_mu, 1 ) );
    signal_map_all.insert( std::make_pair( cc_0pi_pi, 0 ) );

    if(e.CheckMCTopology(signal_map_all) && e.CheckRecoTopology(signal_map_all)) correctly_reconstructed++;
    if(e.CheckMCTopology(signal_map_all)) true_topology++;
    if(e.CheckRecoTopology(signal_map_all)) reco_topology++;
  }

  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Fraction of correctly reconstructed CC 0Pi events : " << correctly_reconstructed/double(true_topology) << std::endl;
  std::cout << " Fraction of mis-identified CC 0Pi events          : " << (reco_topology - correctly_reconstructed)/double(reco_topology) << std::endl; 

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;

  return 0;

} // MainTest()
