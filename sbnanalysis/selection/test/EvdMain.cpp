#include "../include/EventSelectionTool.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/Event.h"
#include <iostream>
#include <iomanip>
#include <fstream>
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
  std::cout << " Start: Local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Running all files " << std::endl;
  
  std::string evd_location = "../Output_Selection_Tool/evd/";
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
    
  }
  std::cout << std::endl;
  
  TopologyMap cc0pi   = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc0pi2p = GeneralAnalysisHelper::GetCC0Pi2PTopologyMap();
  TopologyMap cc0pi3p = GeneralAnalysisHelper::GetCC0Pi3PTopologyMap();
  TopologyMap cc0pi5p = GeneralAnalysisHelper::GetCC0Pi5PTopologyMap();
 
  // Initialise the file to hold file and event ids for different topologies 
  ofstream file;
  file.open(evd_location+"event_display_ids.txt");
  file << std::endl;
  file << "-------------------------------------------------------------------------------------------------" << std::endl;
  file << std::endl;
  file << std::setw(16) << " Type " << std::setw(16) << "MC, Reco, Signal" << std::setw(8) << " File " << std::setw(8) << " Event " << std::endl;
  file << std::endl;
  file << "-------------------------------------------------------------------------------------------------" << std::endl;
  file << std::endl;

  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);

    if(e.CheckMCTopology(cc0pi2p) && e.IsSBNDTrueFiducial()) {
      file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " MC " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;

      if(e.CheckRecoTopology(cc0pi2p)) {
        file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;

        bool found_one = false;
        TVector3 dir1, dir2;
        // Find a hammer event
        for(const Particle &p : e.GetRecoParticleList()){
          if(p.GetPdgCode() == 2212){
            if(!found_one){
              dir1 = (1. / (p.GetEnd() - p.GetVertex()).Mag())*(p.GetEnd() - p.GetVertex());
              found_one = true;
            }
            else{
              dir2 = (1. / (p.GetEnd() - p.GetVertex()).Mag())*(p.GetEnd() - p.GetVertex());
              break;
            }
          }
        }
        // Find the angle between the two
        double cosangle = (1. / (dir1.Mag() * dir2.Mag())) * dir1.Dot(dir2);
        if(cosangle >= -1 && cosangle < -0.95) {
          file      << std::setw(16) <<" HAMMER " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
          std::cout << std::setw(16) <<" HAMMER " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
        }
      }
    }
    if(e.CheckRecoTopology(cc0pi2p) && e.IsSBNDTrueFiducial())
      file << std::setw(16) <<" CC 0Pi 2P " << std::setw(16) << " Reco " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;

    // CC0Pi Signal
    if(e.CheckMCTopology(cc0pi) && e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(cc0pi))
        file << std::setw(16) << " CC 0Pi " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
    }
    // CC0Pi3 Signal
    if(e.CheckMCTopology(cc0pi3p) && e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(cc0pi3p))
        file << std::setw(16) << " CC 0Pi 3P " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
    }
    // CC0Pi3 Signal
    if(e.CheckMCTopology(cc0pi5p) && e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(cc0pi5p)){
        file      << std::setw(16) << " CC 0Pi 5P " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
        std::cout << std::setw(16) << " CC 0Pi 5P " << std::setw(16) << " Signal " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
      }
    }
    // Muon-pion mixup 
    if(e.CheckMCTopology(cc0pi) && e.IsSBNDTrueFiducial() && !e.CheckRecoTopology(cc0pi)){
      // Get the true muon and find out if it was reconstructed as a pion
      for(const Particle &mc : e.GetMCParticleList()){
        if(mc.GetPdgCode() == 13){
          for(const Particle &re : e.GetRecoParticleList()){
            if(!re.GetFromRecoTrack()) continue;
            if(mc.GetMCId() == re.GetMCParticleIdHits() && (re.GetPdgCode() == 211 || re.GetPdgCode() == -211))
              file << std::setw(16) << " CC 0Pi, reco pi" << std::setw(16) << " MC " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
          }
        }
      }
    }
    if(e.CheckRecoTopology(cc0pi) && e.IsSBNDTrueFiducial() && !e.CheckMCTopology(cc0pi)){
      // Get the true muon and find out if it was reconstructed as a pion
      for(const Particle &re : e.GetRecoParticleList()){
        if(!re.GetFromRecoTrack()) continue;
        if(re.GetPdgCode() == 13){
          for(const Particle &mc : e.GetMCParticleList()){
            if(mc.GetMCId() == re.GetMCParticleIdHits() && (mc.GetPdgCode() == 211 || mc.GetPdgCode() == -211))
              file << std::setw(16) << " 1pi, reco mu" << std::setw(16) << " Reco " << std::setw(8) << e.GetFileId() << std::setw(8) << e.GetId() << std::endl;
          }
        }
      }
    }
  }
 
  file.close();

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
