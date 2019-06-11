#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include "../include/Plane.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include <stdexcept>
#include "TROOT.h"
#include "TMath.h"
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
  std::string plots_location = "../Output_Selection_Tool/plots/mcp0_9/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  TopologyMap cc0pi_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap ccinc_map = GeneralAnalysisHelper::GetCCIncTopologyMap();
 
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 991;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    //name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    name = "/home/rhiannon/Samples/LocalSamples/analysis/200219_neutrino_only/selection/"+std::to_string(i)+"/output_file.root";
    //name = "/hepstore/rjones/Samples/FNAL/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
 
  TH2D *h_muon_momentum_angle_ccinc = new TH2D("h_muon_momentum_angle_ccinc","Selected muon momentum vs. angle", 40,0,2, 50, -1, 1);
  TH2D *h_muon_momentum_angle_cc0pi = new TH2D("h_muon_momentum_angle_cc0pi","Selected muon momentum vs. angle", 40,0,2, 50, -1, 1);
  
  for(const Event &e : events){
    //Counter for event-based track counting
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e) || !e.IsSBNDTrueFiducial()) continue;
    if(e.CheckRecoTopology(ccinc_map)){
      for(const Particle &p : e.GetRecoParticleList()){
        if(p.GetPdgCode() == 13){
          h_muon_momentum_angle_ccinc->Fill(p.GetModulusMomentum(),p.GetCosTheta());
        } 
      }
    }
    if(e.CheckRecoTopology(cc0pi_map)){
      for(const Particle &p : e.GetRecoParticleList()){
        if(p.GetPdgCode() == 13){
          h_muon_momentum_angle_cc0pi->Fill(p.GetModulusMomentum(),p.GetCosTheta());
        } 
      }
    }
  }

  /*
   *
   *      OUTPUTS
   *
   */
  // Length vs dist 2D plot
  TCanvas *c = new TCanvas("c","",800,600);

  gStyle->SetPalette(57);
  h_muon_momentum_angle_ccinc->SetContour(250);
  h_muon_momentum_angle_ccinc->GetXaxis()->SetTitle("Muon momentum [GeV/c]");
  h_muon_momentum_angle_ccinc->GetYaxis()->SetTitle("cos(#theta)");
  h_muon_momentum_angle_ccinc->SetStats(0);
  h_muon_momentum_angle_ccinc->Draw("colz");

  c->SaveAs((plots_location+"muon_momentum_angle_ccinc.png").c_str());
  c->SaveAs((plots_location+"muon_momentum_angle_ccinc.root").c_str());

  c->Clear();
  
  h_muon_momentum_angle_cc0pi->SetContour(250);
  h_muon_momentum_angle_cc0pi->GetXaxis()->SetTitle("Muon momentum [GeV/c]");
  h_muon_momentum_angle_cc0pi->GetYaxis()->SetTitle("cos(#theta)");
  h_muon_momentum_angle_cc0pi->SetStats(0);
  h_muon_momentum_angle_cc0pi->Draw("colz");

  c->SaveAs((plots_location+"muon_momentum_angle_cc0pi.png").c_str());
  c->SaveAs((plots_location+"muon_momentum_angle_cc0pi.root").c_str());

  c->Clear();
  
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
