#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
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
#include "TROOT.h"
#include "TAxis.h"

using namespace selection;

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<unsigned int> &exceptions);

int MainTest(){
  
  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string file_location  = "../Output_Selection_Tool/files/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  // Maps
  TopologyMap ccinc = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap cc0pi = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap cc2pi = GeneralAnalysisHelper::GetCC2PiTopologyMap();
  TopologyMap ccpi0 = GeneralAnalysisHelper::GetCCPi0TopologyMap();

  // Variables to get
  bool iscc, isnc;  // Charged or neutral current event current
  int nu_pdg, mode; // Neutrino pdg code and scattering mode
  double enu_true, enu_reco, qsqr; // Neutrino kinematics, truth and reco
  double mu_momentum, mu_cos_z; // Muon kinematics, reco
  int nkaons, npip, npim, npi0; // Particle counting, truth
  // Enumeration (all CC):
  //   -1 = undefined
  //    0 = 0pi
  //    1 = 1chpi
  //    2 = 2chpi
  //    3 = 1pi0
  //    4 = other
  int true_topology = -1;
  int reco_topology = -1;
  
  double pot; // Subrun information

  TTree *t_run    = new TTree("valor_tree"," Tree to hold variables needed for the VALOR analysis");
  TTree *t_subrun = new TTree("subrun_tree"," Tree to hold subrun variables, such as POT");
  t_run->SetDirectory(0);
  t_subrun->SetDirectory(0);
  
  t_run->Branch("iscc",          &iscc,          "iscc/O");
  t_run->Branch("isnc",          &isnc,          "isnc/O");
  t_run->Branch("nu_pdg",        &nu_pdg,        "nu_pdg/I");
  t_run->Branch("enu_true",      &enu_true,      "enu_true/D");
  t_run->Branch("enu_reco",      &enu_reco,      "enu_reco/D");
  t_run->Branch("mu_momentum",   &mu_momentum,   "mu_momentum/D");
  t_run->Branch("mu_cos_z",      &mu_cos_z,      "mu_cos_z/D");
  t_run->Branch("qsqr",          &qsqr,          "qsqr/D");
  t_run->Branch("mode",          &mode,          "mode/I");
  t_run->Branch("nkaons",        &nkaons,        "nkaons/I");
  t_run->Branch("npip",          &npip,          "npip/I");
  t_run->Branch("npim",          &npim,          "npim/I");
  t_run->Branch("npi0",          &npi0,          "npi0/I");
  t_run->Branch("true_topology", &true_topology, "true_topology/I");
  t_run->Branch("reco_topology", &reco_topology, "reco_topology/I");

  t_subrun->Branch("pot",        &pot,           "pot/D");
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 991;
  std::vector<unsigned int> exceptions;
  exceptions.clear();

  // Read in txt file of list of empty input directories
  std::fstream exception_file("exceptions.txt");
  if(!exception_file)
    std::cout << " File not open! " << std::endl;
  std::string s_exc;
  while (std::getline(exception_file, s_exc)) {
    unsigned int i_exc;
    std::istringstream ss_exc(s_exc);
    ss_exc >> i_exc;
    exceptions.push_back(i_exc);
    ss_exc.str(std::string());
    s_exc.clear();
  }

  std::cout << " Skipping files in directory : " << std::endl;
  for(const unsigned int & ex : exceptions){
    if(ex < total_files)
      std::cout << " - " << ex << " - ";
  }
  std::cout << std::endl;

  LoadAllEvents(events, total_files, start, pot, exceptions);
  

  // Counter to quantify how many events have strange particle energies in them
  int bad_events = 0;

  // Loop over events and perform vertexing study
  for(const Event &e : events){
    bool bad_event = false;

    // TPC track criteria
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;
    iscc     = e.GetIsCC();
    isnc     = !e.GetIsCC();
    nu_pdg   = e.GetNeutrinoPdgCode();
    enu_true = e.GetTrueNuEnergy();
    qsqr     = e.GetTrueNuQ2();
    mode     = e.GetPhysicalProcess();

    // Start the counters
    nkaons = 0;
    npip   = 0;
    npim   = 0;
    npi0   = 0;
    enu_reco = 0.;
    
    // Particles 
    ParticleList mc   = e.GetMCParticleList();
    ParticleList reco = e.GetRecoParticleList();

    // Neutrino vertex is within the ficudial border
    if(e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(ccinc)){
        if(!bad_event){
          for(const Particle &p : reco){
            if(p.GetKineticEnergy() <= 0.) continue;
            if(p.GetKineticEnergy() > 10.) {
              bad_event = true;
              bad_events++;
              std::cerr << " Badly defined energy of the particle, skipping event " << std::endl;
              break;
            }
            if(p.GetPdgCode() == 13){
              mu_momentum = p.GetModulusMomentum();
              mu_cos_z    = p.GetCosTheta();
            } 
            enu_reco += (p.GetKineticEnergy() + p.GetMass());
          }
          for(const Particle &p : mc){
            if(p.GetPdgCode() ==  311 || p.GetPdgCode() == -321 || p.GetPdgCode() == 321) nkaons++;
            if(p.GetPdgCode() ==  211) npip++;
            if(p.GetPdgCode() == -211) npim++;
            if(p.GetPdgCode() ==  111) npi0++; 
          } // Particles
          if(e.CheckRecoTopology(cc0pi))
            reco_topology = 0;
          else if(e.CheckRecoTopology(cc1pi))
            reco_topology = 1;
          else if(e.CheckRecoTopology(cc2pi))
            reco_topology = 2;
          else if(e.CheckRecoTopology(ccpi0))
            reco_topology = 3;
          else
            reco_topology = 4;
          if(e.CheckMCTopology(cc0pi))
            true_topology = 0;
          else if(e.CheckMCTopology(cc1pi))
            true_topology = 1;
          else if(e.CheckMCTopology(cc2pi))
            true_topology = 2;
          else if(e.CheckMCTopology(ccpi0))
            true_topology = 3;
          else
            true_topology = 4;
          t_run->Fill();
        } // Bad event
      } // Topology
    } // Fiducial
  } // Events
  std::cout << " Total events with particles with bad energy : " << bad_events << std::endl;
  
  // Print the total pot from all the samples
  std::cout << " Total POT in the samples is: " << pot << std::endl;
  t_subrun->Fill();

  // Output TFile
  TFile f((file_location+"test_ccinc_selection.root").c_str(), "RECREATE");

  t_run->Write();
  t_subrun->Write();

  f.Write();
  f.Close();

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

} // MainTest

void LoadAllEvents(EventSelectionTool::EventList &events, const unsigned int &total_files, const int &start_time, double &pot, std::vector<unsigned int> &exceptions) {
  double total_pot = 0;
  std::vector<unsigned int>::iterator it;
  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    it = std::find(exceptions.begin(), exceptions.end(),i);
    if(it != exceptions.end()) continue;
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
//    name = "/home/rhiannon/Samples/LocalSamples/analysis/test/output_file.root";
      name = "/pnfs/sbnd/persistent/users/rsjones/mcp0.9_neutrino_with_subrun/selection/"+std::to_string(i)+"/output_file.root";
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
