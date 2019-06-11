#include "../include/CC0piAnalysisHelper.h"
#include "../include/CC1piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
using namespace selection;

int MainTest(){

  time_t rawtime;
  struct tm * timeinfo;
  time (&rawtime);
  timeinfo = localtime (&rawtime);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " Start: Local time and date:  " << asctime(timeinfo)         << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  std::vector< double > CountMC    ( 5 ) ;
  std::vector< double > CountTReco ( 5 ) ;
  std::vector< double > CountReco  ( 5 ) ;

  std::vector< vector< double> > Count_MC_Topology    ( 5, vector< double > (5) );
  std::vector< vector< double> > Count_TReco_Topology ( 5, vector< double > (5) );
  std::vector< vector< double> > Count_Reco_Topology  ( 5, vector< double > (5) );

  TH1 * h_mu_MC_L       = new TH1D ( "h_mu_MC_L"      , " TRUE: mu track length      "    , 20, 0., 400.  );
  TH1 * h_mu_TReco_L    = new TH1D ( "h_mu_TReco_L"   , " SIGNAL: mu track length    "    , 20, 0., 400.  );
  TH1 * h_mu_Reco_L     = new TH1D ( "h_mu_Reco_L"    , " SELECTED: mu track length  "    , 20, 0., 400.  );
  
  TH1 * h_pi_MC_L       = new TH1D ( "h_pi_MC_L"      , " TRUE: pi track length      "    , 20, 0., 200.  );
  TH1 * h_pi_TReco_L    = new TH1D ( "h_pi_TReco_L"   , " SIGNAL: pi track length    "    , 20, 0., 200.  );
  TH1 * h_pi_Reco_L     = new TH1D ( "h_pi_Reco_L"    , " SELECTED: pi track length  "    , 20, 0., 200.  );
  
  TH1 * h_pi0_MC_L      = new TH1D ( "h_pi0_MC_L"     , " TRUE: pi0 track length     "    , 20, 0., 200.  );
  TH1 * h_pi0_TReco_L   = new TH1D ( "h_pi0_TReco_L"  , " SIGNAL: pi0 track length   "    , 20, 0., 200.  );
  TH1 * h_pi0_Reco_L    = new TH1D ( "h_pi0_Reco_L"   , " SELECTED: pi0 track length "    , 20, 0., 200.  );

  TH1 * h_mu_MC_T       = new TH1D ( "h_mu_MC_T"      , " TRUE: mu track angle       "    , 20, -1.,  1.  );
  TH1 * h_mu_TReco_T    = new TH1D ( "h_mu_TReco_T"   , " SIGNAL: mu track angle     "    , 20, -1.,  1.  );
  TH1 * h_mu_Reco_T     = new TH1D ( "h_mu_Reco_T"    , " SELECTED: mu track angle   "    , 20, -1.,  1.  );
  
  TH1 * h_pi_MC_T       = new TH1D ( "h_pi_MC_T"      , " TRUE: pi track angle       "    , 20, -1.,  1.  );
  TH1 * h_pi_TReco_T    = new TH1D ( "h_pi_TReco_T"   , " SIGNAL: pi track angle     "    , 20, -1.,  1.  );
  TH1 * h_pi_Reco_T     = new TH1D ( "h_pi_Reco_T"    , " SELECTED: pi track angle   "    , 20, -1.,  1.  );

  TH1 * h_pi0_MC_T      = new TH1D ( "h_pi0_MC_T"     , " TRUE: pi0 track angle      "    , 20, -1.,  1.  );
  TH1 * h_pi0_TReco_T   = new TH1D ( "h_pi0_TReco_T"  , " SIGNAL: pi0 track angle    "    , 20, -1.,  1.  );
  TH1 * h_pi0_Reco_T    = new TH1D ( "h_pi0_Reco_T"   , " SELECTED: pi0 track angle  "    , 20, -1.,  1.  );

  TH1 * h_MC_E_mu       = new TH1D ( "h_MC_E_mu"      , " TRUE: E for mu             "    , 100,  0.,  3.  );   
  TH1 * h_TReco_E_mu    = new TH1D ( "h_TReco_E_mu"   , " SIGNAL: E for mu           "    , 100,  0.,  3.  );  
  TH1 * h_Reco_E_mu     = new TH1D ( "h_Reco_E_mu"    , " SELECTED: E for mu         "    , 100,  0.,  3.  );   
  
  TH1 * h_MC_E_pi       = new TH1D ( "h_MC_E_pi"      , " TRUE: E for pi             "    , 100,  0.,  3.  );   
  TH1 * h_TReco_E_pi    = new TH1D ( "h_TReco_E_pi"   , " SIGNAL: E for pi           "    , 100,  0.,  3.  );  
  TH1 * h_Reco_E_pi     = new TH1D ( "h_Reco_E_pi"    , " SELECTED: E for pi         "    , 100,  0.,  3.  ); 

  TH1 * h_MC_Delta      = new TH1D ( "h_MC_Delta"     ," TRUE: Delta mass            "    , 20 , 1.,  1.5 );   
  TH1 * h_TReco_Delta   = new TH1D ( "h_TReco_Delta"  , " SIGNAL: Delta mass         "    , 20, 0.,  1.5 );  
  TH1 * h_Reco_Delta    = new TH1D ( "h_Reco_Delta"   , " SELECTED: Delta mass       "    , 20 , 1.,  1.5 ); 

  THStack * hS_mu_Reco_L  = new THStack ( "hS_mu_Reco_L"    , "Reco: mu track length "  );
  TH1 * h_mu_Reco_cc0pi_L = new TH1D ( "h_mu_Reco_cc0pi_L"    , "Reco: mu track length "    , 20, 0., 400.  );
  TH1 * h_mu_Reco_cc1pi_L = new TH1D ( "h_mu_Reco_cc1pi_L"    , "Reco: mu track length "    , 20, 0., 400.  );
  TH1 * h_mu_Reco_ccpi0_L = new TH1D ( "h_mu_Reco_ccpi0_L"    , "Reco: mu track length "    , 20, 0., 400.  );
  TH1 * h_mu_Reco_nc_L    = new TH1D ( "h_mu_Reco_nc_L"       , "Reco: mu track length "    , 20, 0., 400.  );

  THStack * hS_mu_Reco_T  = new THStack ( "hS_mu_Reco_T"    , "Reco: mu track angle "    );
  TH1 * h_mu_Reco_cc0pi_T = new TH1D( "h_mu_Reco_cc0pi_T","sthn", 20,-1.,1);
  TH1 * h_mu_Reco_cc1pi_T = new TH1D( "h_mu_Reco_cc1pi_T","sthn", 20,-1.,1);
  TH1 * h_mu_Reco_ccpi0_T = new TH1D( "h_mu_Reco_ccpi0_T","sthn", 20,-1.,1);
  TH1 * h_mu_Reco_nc_T    = new TH1D( "h_mu_Reco_nc_T","sthn"   , 20,-1.,1);
  
  THStack * hS_mu_Reco_E  = new THStack ( "hS_mu_Reco_E"    , "Reco: mu track angle "    );
  TH1 * h_mu_Reco_cc0pi_E = new TH1D( "h_mu_Reco_cc0pi_E","sthn", 20,0.,3);
  TH1 * h_mu_Reco_cc1pi_E = new TH1D( "h_mu_Reco_cc1pi_E","sthn", 20,0.,3);
  TH1 * h_mu_Reco_ccpi0_E = new TH1D( "h_mu_Reco_ccpi0_E","sthn", 20,0.,3);
  TH1 * h_mu_Reco_nc_E    = new TH1D( "h_mu_Reco_nc_E","sthn"   , 20,0.,3);

  THStack * hS_pi_Reco_L   = new THStack ( "hS_pi_Reco_L"   , "Reco: pi track length "  );
  TH1 * h_pi_Reco_cc0pi_L  = new TH1D ( "h_pi_Reco_cc0pi_L" , "Reco: pi track length "    , 20, 0., 200.  );
  TH1 * h_pi_Reco_cc1pi_L  = new TH1D ( "h_pi_Reco_cc1pi_L" , "Reco: pi track length "    , 20, 0., 200.  );
  TH1 * h_pi_Reco_ccpi0_L  = new TH1D ( "h_pi_Reco_ccpi0_L" , "Reco: pi track length "    , 20, 0., 200.  );
  TH1 * h_pi_Reco_nc_L     = new TH1D ( "h_pi_Reco_nc_L"    , "Reco: pi track length "    , 20, 0., 200.  );

  THStack * hS_pi_Reco_T  = new THStack ( "hS_pi_Reco_T"    , "Reco: pi track angle "    );
  TH1 * h_pi_Reco_cc0pi_T = new TH1D( "h_pi_Reco_cc0pi_T","sthn", 20,-1.,1);
  TH1 * h_pi_Reco_cc1pi_T = new TH1D( "h_pi_Reco_cc1pi_T","sthn", 20,-1.,1);
  TH1 * h_pi_Reco_ccpi0_T = new TH1D( "h_pi_Reco_ccpi0_T","sthn", 20,-1.,1);
  TH1 * h_pi_Reco_nc_T    = new TH1D( "h_pi_Reco_nc_T","sthn"   , 20,-1.,1);
  
  THStack * hS_pi_Reco_E  = new THStack ( "hS_pi_Reco_E"    , "Reco: pi track angle "    );
  TH1 * h_pi_Reco_cc0pi_E = new TH1D( "h_pi_Reco_cc0pi_E","sthn", 20,0.,3);
  TH1 * h_pi_Reco_cc1pi_E = new TH1D( "h_pi_Reco_cc1pi_E","sthn", 20,0.,3);
  TH1 * h_pi_Reco_ccpi0_E = new TH1D( "h_pi_Reco_ccpi0_E","sthn", 20,0.,3);
  TH1 * h_pi_Reco_nc_E    = new TH1D( "h_pi_Reco_nc_E","sthn"   , 20,0.,3);

  THStack * hS_pi0_Reco_L   = new THStack ( "hS_pi0u_Reco_L"  , "Reco: pi0 track length "  );
  TH1 * h_pi0_Reco_cc0pi_L  = new TH1D ( "h_pi0_Reco_cc0pi_L" , "Reco: pi0 track length "    , 20, 0., 200.  );
  TH1 * h_pi0_Reco_cc1pi_L  = new TH1D ( "h_pi0_Reco_cc1pi_L" , "Reco: pi0 track length "    , 20, 0., 200.  );
  TH1 * h_pi0_Reco_ccpi0_L  = new TH1D ( "h_pi0_Reco_ccpi0_L" , "Reco: pi0 track length "    , 20, 0., 200.  );
  TH1 * h_pi0_Reco_nc_L     = new TH1D ( "h_pi0_Reco_nc_L"    , "Reco: pi0 track length "    , 20, 0., 200.  );

  THStack * hS_pi0_Reco_T  = new THStack ( "hS_pi0_Reco_T"    , "Reco: pi0 track angle "    );
  TH1 * h_pi0_Reco_cc0pi_T = new TH1D( "h_pi0_Reco_cc0pi_T" ,"sthn", 20,-1.,1);
  TH1 * h_pi0_Reco_cc1pi_T = new TH1D( "h_pi0_Reco_cc1pi_T" ,"sthn", 20,-1.,1);
  TH1 * h_pi0_Reco_ccpi0_T = new TH1D( "h_pi0u_Reco_ccpi0_T","sthn", 20,-1.,1);
  TH1 * h_pi0_Reco_nc_T    = new TH1D( "h_pi0  _Reco_nc_T"  ,"sthn", 20,-1.,1);

  TH1F *h_MC_energy_nu_cc0pi    = new TH1F("h_MC_energy_nu_cc0pi", "TRUE : #nu_{#mu} CC 0#pi neutrino energy"     , 50, 0, 3);
  TH1F *h_Reco_energy_nu_cc0pi  = new TH1F("h_Reco_energy_nu_cc0pi", "SELECTED: #nu_{#mu} CC 0#pi neutrino energy", 50, 0, 3);
  TH1F *h_TReco_energy_nu_cc0pi = new TH1F("h_TReco_energy_nu_cc0pi", "SIGNAL: #nu_{#mu} CC 0#pi neutrino energy" , 50, 0, 3);

  TH1F *h_MC_energy_nu_cc1pi    = new TH1F("h_MC_energy_nu_cc1pi", "TRUE : #nu_{#mu} CC 1#pi neutrino energy"     , 50, 0, 3);
  TH1F *h_Reco_energy_nu_cc1pi  = new TH1F("h_Reco_energy_nu_cc1pi", "SELECTED: #nu_{#mu} CC 1#pi neutrino energy", 50, 0, 3);
  TH1F *h_TReco_energy_nu_cc1pi = new TH1F("h_TReco_energy_nu_cc1pi", "SIGNAL: #nu_{#mu} CC 1#pi neutrino energy" , 50, 0, 3);

  TH1F *h_MC_energy_nu_cc1pi_M2    = new TH1F("h_MC_energy_nu_cc1pi_M2", "TRUE : #nu_{#mu} CC 1#pi neutrino energy"     , 50, 0, 3);
  TH1F *h_Reco_energy_nu_cc1pi_M2  = new TH1F("h_Reco_energy_nu_cc1pi_M2", "SELECTED: #nu_{#mu} CC 1#pi neutrino energy", 50, 0, 3);
  TH1F *h_TReco_energy_nu_cc1pi_M2 = new TH1F("h_TReco_energy_nu_cc1pi_M2", "SIGNAL: #nu_{#mu} CC 1#pi neutrino energy" , 50, 0, 3);

  //
  TH2 * h_E_pi_MC_TReco = new TH2D ( "h_E_pi_MC_TReco", " TE Signal vs True          "    , 20,  0., 3., 20, -1., 1. );
  TH2 * h_Enu_MC_TReco = new TH2D ( "h_Enu_MC_TReco", " TE Signal vs True          "    , 20,  0., 3., 20, -1., 1. );
  TH2 * h_Q_MC_TReco  = new TH2D ( "h_Q_MC_TReco", " TE Signal vs True          "    , 20,  -1., 1., 20, -1., 1. );
  TH1 * h_R_E_mu        = new TH1D ( "h_MC_E_mu_R"    , " T((Signal - True )/True)   "    , 20,  -1., 1. );   
  TH1 * h_R_E_pi        = new TH1D ( "h_MC_E_pi_R"    , " T((Signal - True )/True)   "    , 20,  -1., 1. ); 
  //
  
  ofstream efile, lfile, afile ;
  efile.open( "event_information.txt" ) ;
  lfile.open( "length_information_per_event.txt" ) ;
  afile.open( "angular_information_per_event.txt" );
  float RecoEPi = 0;
  double count_mis_E_Reco=0, count_mis_E_TReco=0;
  //Initialise event list and the topology maps
 
  EventSelectionTool::EventList events;
  for( unsigned int i = 0; i < 500; ++i ){

    // Get the filename for each 2D histogram
    std::stringstream ss;
    ss.clear();
    
    std::string name;
    name.clear();
    
    char file_name[1024];
    
    ss << "/hepstore/rjones/Samples/FNAL/sbn_workshop_0318_new/4883618_" << i <<"/output_file.root";
    name = ss.str();
            
    strcpy( file_name, name.c_str() );
      
    EventSelectionTool::LoadEventList(file_name, events);
    }
  
  time_t rawtime_afterload;
  struct tm * timeinfo_afterload;
  time (&rawtime_afterload);
  timeinfo_afterload = localtime (&rawtime_afterload);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " After loading events: Local time and date:  " << asctime(timeinfo_afterload) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  // Histograms directory :
  std::string filepath = "../Output_Selection_Tool/plots/cc1pi_test/";

  for(unsigned int i = 0; i < events.size(); ++i){
    //---------------------------- Do analysis---------------------------------------------
    Event &e(events[i]);
    TopologyMap topology = GeneralAnalysisHelper::GetCC1PiTopologyMap();

    // ----------------------------Save Event Information----------------------------------
    GeneralAnalysisHelper::EventInformationParticles( e, "../Output_Selection_Tool/statistics/event_information.txt" , i );
    //GeneralAnalysisHelper::EventProperties( e, topology, "../Output_Selection_Tool/statistics/event_properties",  i );
    
    // ----------------------------Eficiency calculation values----------------------------

    GeneralAnalysisHelper::TopologyStatistics(e, GeneralAnalysisHelper::GetNCTopologyMap(),    CountMC[0], CountTReco[0], CountReco[0]);
    GeneralAnalysisHelper::TopologyStatistics(e, GeneralAnalysisHelper::GetCCIncTopologyMap(), CountMC[1], CountTReco[1], CountReco[1]);
    GeneralAnalysisHelper::TopologyStatistics(e, GeneralAnalysisHelper::GetCC0PiTopologyMap(), CountMC[2], CountTReco[2], CountReco[2]);
    GeneralAnalysisHelper::TopologyStatistics(e, GeneralAnalysisHelper::GetCC0PiTopologyMap(), CountMC[3], CountTReco[3], CountReco[3]);
    GeneralAnalysisHelper::TopologyStatistics(e, GeneralAnalysisHelper::GetCCPi0TopologyMap(), CountMC[4], CountTReco[4], CountReco[4]);
    
    //--------------------- TOPOLOGY MIS IDENTIFICATION : TOPOLOGY MATRIX -----------------
    GeneralAnalysisHelper::TopologyMatrix(e, Count_MC_Topology, Count_TReco_Topology, Count_Reco_Topology);
   
    /*************************************************************************************
     *    
     *    Commented out while we decide how to implement the new (fixed) functions
     *
     *************************************************************************************
     *
    //--------------------- Events vs Length, Angle and Kinetic Energy   ------------------
    if( e.CheckMCTopology( topology ) == 1 ) { 

      if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     && e.GetMCLengthWithPdg( 13 )  != 0 )  h_mu_MC_L -> Fill( e.GetMCLengthWithPdg( 13 )  );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetMCLengthWithPdg( 211 ) != 0 )  h_pi_MC_L -> Fill( e.GetMCLengthWithPdg( 211 ) );
      if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() && e.GetMCLengthWithPdg( 111 ) != 0 )  h_pi0_MC_L-> Fill( e.GetMCLengthWithPdg( 111 ) );
      
      if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     && e.GetMCCosThetaWithPdg( 13  ) !=0 ) h_mu_MC_T -> Fill( e.GetMCCosThetaWithPdg( 13 ) );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetMCCosThetaWithPdg( 211 ) !=0 ) h_pi_MC_T -> Fill( e.GetMCCosThetaWithPdg( 211 ) );
      if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() && e.GetMCCosThetaWithPdg( 111 ) !=0 ) h_pi0_MC_T -> Fill( e.GetMCCosThetaWithPdg( 111 ) );
      
      if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     ) h_MC_E_mu  -> Fill( e.GetMCKineticEnergyWithPdg( 13 ) );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) h_MC_E_pi  -> Fill( e.GetMCKineticEnergyWithPdg( 211 ) );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) h_MC_Delta -> Fill( e.GetMCDeltaEnergy() );

      if( e.CheckRecoTopology( topology ) == 1 ) {
        if( e.CheckRecoTopology( GeneralAnalysisHelper::GetCCPi0TopologyMap() ) == 1 && e.GetRecoCC0piNeutrinoEnergy()<0 ){}else if(e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 && e.GetRecoCC1piNeutrinoEnergy()<0 ){std::cout<<"cucu";} 
        else{
          if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     && e.GetRecoLengthWithPdg( 13   ) != 0 ) h_mu_TReco_L -> Fill( e.GetRecoLengthWithPdg( 13 )   );
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoLengthWithPdg( 211  ) != 0 && e.GetRecoCC1piNeutrinoEnergy()>0 ) h_pi_TReco_L -> Fill( e.GetRecoLengthWithPdg( 211 )  );
          if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() && e.GetRecoLengthWithPdg( 111  ) != 0 ) h_pi0_TReco_L-> Fill( e.GetRecoLengthWithPdg( 111 )  );
          
          if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     && e.GetRecoCosThetaWithPdg( 13  ) != 0 ) h_mu_TReco_T  -> Fill( e.GetRecoCosThetaWithPdg( 13 )  );
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoCosThetaWithPdg( 211 ) != 0 && e.GetRecoCC1piNeutrinoEnergy()>0 ) h_pi_TReco_T  -> Fill( e.GetRecoCosThetaWithPdg( 211 ) );
          if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() && e.GetRecoCosThetaWithPdg( 111 ) != 0 ) h_pi0_TReco_T -> Fill( e.GetRecoCosThetaWithPdg( 111 ) );
          
          if( topology != GeneralAnalysisHelper::GetNCTopologyMap()     && e.GetRecoKineticEnergyWithPdg( 13 ) !=0 ) h_TReco_E_mu-> Fill( e.GetRecoKineticEnergyWithPdg( 13 ) );
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoKineticEnergyWithPdg( 211 )!=0 && e.GetRecoCC1piNeutrinoEnergy()>0 ) h_TReco_E_pi-> Fill( e.GetRecoKineticEnergyWithPdg( 211 ) );
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoKineticEnergyWithPdg( 13 )!=0 && e.GetMCKineticEnergyWithPdg( 13 )!=0 && e.GetRecoCC1piNeutrinoEnergy()>0){
            h_R_E_mu->Fill( ( e.GetRecoKineticEnergyWithPdg( 13 )- e.GetMCKineticEnergyWithPdg( 211 )) / e.GetMCKineticEnergyWithPdg( 13 ) );	 }
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoKineticEnergyWithPdg( 211 )!=0 && e.GetMCKineticEnergyWithPdg( 211 )!=0 && e.GetRecoCC1piNeutrinoEnergy()>0 ){
            h_R_E_pi->Fill( ( e.GetRecoKineticEnergyWithPdg( 211 )- e.GetMCKineticEnergyWithPdg( 211 )) / e.GetMCKineticEnergyWithPdg( 211 ) );}
          if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && e.GetRecoKineticEnergyWithPdg( 211 )!=0  && e.GetRecoCC1piNeutrinoEnergy()>0 ){
              h_Q_MC_TReco->Fill( e.GetMCQ2WithPdg_cc1pi( 13) ,  (e.GetRecoQ2WithPdg_cc1pi( 13)-e.GetMCQ2WithPdg_cc1pi( 13) )/  e.GetMCQ2WithPdg_cc1pi( 13) );}
        }
      }
    }
    if( e.CheckRecoTopology( topology ) == 1) { 
      if( (e.CheckRecoTopology( GeneralAnalysisHelper::GetCCPi0TopologyMap() ) == 1 && e.GetRecoCC0piNeutrinoEnergy()<0) || (e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 && e.GetRecoCC1piNeutrinoEnergy()<0) ){} 
      else{
        if( topology != GeneralAnalysisHelper::GetNCTopologyMap()      && e.GetRecoLengthWithPdg( 13   ) != 0 ) h_mu_Reco_L -> Fill( e.GetRecoLengthWithPdg( 13   ) );
        if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap()  && e.GetRecoLengthWithPdg( 211  ) != 0 ) h_pi_Reco_L -> Fill( e.GetRecoLengthWithPdg( 211  ) );
        if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap()  && e.GetRecoLengthWithPdg( 111  ) != 0 ) h_pi0_Reco_L -> Fill( e.GetRecoLengthWithPdg( 111  ) );
        
        if( topology != GeneralAnalysisHelper::GetNCTopologyMap()      && e.GetRecoCosThetaWithPdg( 13   ) != 0 )h_mu_Reco_T -> Fill( e.GetRecoCosThetaWithPdg( 13   ) );
        if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap()  && e.GetRecoCosThetaWithPdg( 211  ) != 0 ) h_pi_Reco_T -> Fill( e.GetRecoCosThetaWithPdg( 211  ) );
        if( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap()  && e.GetRecoCosThetaWithPdg( 111  ) != 0 ) h_pi0_Reco_T -> Fill( e.GetRecoCosThetaWithPdg( 111  ) );
        
        if( topology != GeneralAnalysisHelper::GetNCTopologyMap()      && e.GetRecoKineticEnergyWithPdg( 13   ) != 0 )h_Reco_E_mu-> Fill( e.GetRecoKineticEnergyWithPdg( 13 ) );
        if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap()  && e.GetRecoKineticEnergyWithPdg( 211  ) != 0 ) h_Reco_E_pi-> Fill( e.GetRecoKineticEnergyWithPdg( 211 ) );
        if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap()  && e.GetRecoKineticEnergyWithPdg( 211 )!=0 && e.GetMCKineticEnergyWithPdg( 211 )!=0  ){
          h_R_E_pi->Fill( ( e.GetRecoKineticEnergyWithPdg( 211 )- e.GetMCKineticEnergyWithPdg( 211 ))/e.GetMCKineticEnergyWithPdg(211) );}
      } 
    }

    // ------------------------------------------ Background study ----------------------------------------------------
    if( e.CheckRecoTopology( topology ) == 1 ) { 
      if( topology != GeneralAnalysisHelper::GetNCTopologyMap() ){
	      if( (e.CheckRecoTopology( GeneralAnalysisHelper::GetCCPi0TopologyMap() ) == 1 && e.GetRecoCC0piNeutrinoEnergy()<0) || (e.CheckRecoTopology( GeneralAnalysisHelper::GetCC1PiTopologyMap() ) == 1 && e.GetRecoCC1piNeutrinoEnergy()<0) ){} 
        else{
	        if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap()) == 1) {                                                              
            if( e.GetRecoLengthWithPdg( 13 )!=0        ) h_mu_Reco_cc0pi_L -> Fill( e.GetRecoLengthWithPdg( 13 )        );       
            if( e.GetRecoCosThetaWithPdg( 13 ) != 0    ) h_mu_Reco_cc0pi_T -> Fill( e.GetRecoCosThetaWithPdg( 13 )      );
            if( e.GetRecoKineticEnergyWithPdg( 13 )!=0 ) h_mu_Reco_cc0pi_E -> Fill( e.GetRecoKineticEnergyWithPdg( 13 ) ); }		
            else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) == 1) {                                                              
            if( e.GetRecoLengthWithPdg( 13 )!=0        ) h_mu_Reco_ccpi0_L -> Fill( e.GetRecoLengthWithPdg( 13 )        );       
            if( e.GetRecoCosThetaWithPdg( 13 ) != 0    ) h_mu_Reco_ccpi0_T -> Fill( e.GetRecoCosThetaWithPdg( 13 )      ); 
            if( e.GetRecoKineticEnergyWithPdg( 13 )!=0 ) h_mu_Reco_ccpi0_E -> Fill( e.GetRecoKineticEnergyWithPdg( 13 ) ); }     
            else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1) {                                                              
            if( e.GetRecoLengthWithPdg( 13 )!=0        ) h_mu_Reco_cc1pi_L -> Fill( e.GetRecoLengthWithPdg( 13 )        );       
            if( e.GetRecoCosThetaWithPdg( 13 ) != 0    ) h_mu_Reco_cc1pi_T -> Fill( e.GetRecoCosThetaWithPdg( 13 )      ); 
            if( e.GetRecoKineticEnergyWithPdg( 13 )!=0 ) h_mu_Reco_cc1pi_E -> Fill( e.GetRecoKineticEnergyWithPdg( 13 ) ); }
          else {                                                              
            if( e.GetRecoLengthWithPdg( 13 ) !=0       ) h_mu_Reco_nc_L -> Fill( e.GetRecoLengthWithPdg( 13 )           );       
            if( e.GetRecoCosThetaWithPdg( 13 ) != 0    ) h_mu_Reco_nc_T -> Fill( e.GetRecoCosThetaWithPdg( 13 )         ); 
            if ( e.GetRecoKineticEnergyWithPdg( 13 )!=0) h_mu_Reco_nc_E-> Fill( e.GetRecoKineticEnergyWithPdg( 13 )     ); }   
        } 
      }
      if ( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ){
        if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap()) == 1) {                                                              
	        if( e.GetRecoLengthWithPdg( 211 )!=0        ) h_pi_Reco_cc0pi_L -> Fill( e.GetRecoLengthWithPdg( 211 )        );       
          if( e.GetRecoCosThetaWithPdg( 211 ) != 0    ) h_pi_Reco_cc0pi_T -> Fill( e.GetRecoCosThetaWithPdg( 211 )      );
	        if( e.GetRecoKineticEnergyWithPdg( 211 )!=0 ) h_pi_Reco_cc0pi_E -> Fill( e.GetRecoKineticEnergyWithPdg( 211 ) ); }		
	      else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) == 1) {                                                              
	        if( e.GetRecoLengthWithPdg( 211 )!=0        ) h_pi_Reco_ccpi0_L -> Fill( e.GetRecoLengthWithPdg( 211 )        );       
          if( e.GetRecoCosThetaWithPdg( 211 ) != 0    ) h_pi_Reco_ccpi0_T -> Fill( e.GetRecoCosThetaWithPdg( 211 )      ); 
	        if( e.GetRecoKineticEnergyWithPdg( 211 )!=0 ) h_pi_Reco_ccpi0_E -> Fill( e.GetRecoKineticEnergyWithPdg( 211 ) ); }     
	      else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1) {                                                              
	        if( e.GetRecoLengthWithPdg( 211 )!=0        ) h_pi_Reco_cc1pi_L -> Fill( e.GetRecoLengthWithPdg( 211 )        );       
          if( e.GetRecoCosThetaWithPdg( 211 ) != 0    ) h_pi_Reco_cc1pi_T -> Fill( e.GetRecoCosThetaWithPdg( 211 )      ); 
	        if( e.GetRecoKineticEnergyWithPdg( 211 )!=0 ) h_pi_Reco_cc1pi_E -> Fill( e.GetRecoKineticEnergyWithPdg( 211 ) ); }     
	      else {           
	        if( e.GetRecoLengthWithPdg( 211 ) !=0       ) h_pi_Reco_nc_L -> Fill( e.GetRecoLengthWithPdg( 211 )           );       
          if( e.GetRecoCosThetaWithPdg( 211 ) != 0    ) h_pi_Reco_nc_T -> Fill( e.GetRecoCosThetaWithPdg( 211 )         ); 
	        if( e.GetRecoKineticEnergyWithPdg( 211 ) !=0) h_pi_Reco_nc_E -> Fill( e.GetRecoKineticEnergyWithPdg( 211 )    ); }   
        
      }
      if ( topology == GeneralAnalysisHelper::GetCCPi0TopologyMap() ){
	    if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC0PiTopologyMap()) == 1) {                                                              
	      if( e.GetRecoLengthWithPdg( 111 )!=0        ) h_pi0_Reco_cc0pi_L -> Fill( e.GetRecoLengthWithPdg( 111 )        );       
        if( e.GetRecoCosThetaWithPdg( 111 ) != 0    ) h_pi0_Reco_cc0pi_T -> Fill( e.GetRecoCosThetaWithPdg( 111 )      ); }
	    else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCCPi0TopologyMap()) == 1) {                                                              
	      if( e.GetRecoLengthWithPdg( 111 )!=0        ) h_pi0_Reco_ccpi0_L -> Fill( e.GetRecoLengthWithPdg( 111 )        );       
        if( e.GetRecoCosThetaWithPdg( 111 ) != 0    ) h_pi0_Reco_ccpi0_T -> Fill( e.GetRecoCosThetaWithPdg( 111 )      ); } 
	    else if( e.CheckMCTopology(GeneralAnalysisHelper::GetCC1PiTopologyMap()) == 1) {                                                              
	      if( e.GetRecoLengthWithPdg( 111 )!=0        ) h_pi0_Reco_cc1pi_L -> Fill( e.GetRecoLengthWithPdg( 111 )        );       
        if( e.GetRecoCosThetaWithPdg( 111 ) != 0    ) h_pi0_Reco_cc1pi_T -> Fill( e.GetRecoCosThetaWithPdg( 111 )      ); }     
	    else {                                                              
	      if( e.GetRecoLengthWithPdg( 111 ) !=0       ) h_pi0_Reco_nc_L -> Fill( e.GetRecoLengthWithPdg( 111 )           );       
        if( e.GetRecoCosThetaWithPdg( 111 ) != 0    ) h_pi0_Reco_nc_T -> Fill( e.GetRecoCosThetaWithPdg( 111 )         ); }   
      }
    }
    */
    // ------------------------------------------    Efficiency     ----------------------------------------------------
    if( i == events.size()-1 ) std::cout << "Efficiency for the selected topology" << GeneralAnalysisHelper::Efficiency( CountMC, CountTReco, CountReco, topology ) << std::endl;
    // ------------------------------------------  Neutrino Energy  ----------------------------------------------------
    if( e.CheckMCTopology( topology ) == 1 ) { 
      if( topology == GeneralAnalysisHelper::GetCC0PiTopologyMap() ) h_MC_energy_nu_cc0pi -> Fill( e.GetTrueNuEnergy() );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) h_MC_energy_nu_cc1pi -> Fill( e.GetTrueNuEnergy() ); 
      if( e.CheckRecoTopology( topology ) == 1 ) { 
	      if( topology == GeneralAnalysisHelper::GetCC0PiTopologyMap() && CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e)>0 ) h_TReco_energy_nu_cc0pi -> Fill( CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e) );
	      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e)>0 ) h_TReco_energy_nu_cc1pi -> Fill( CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e) );
	      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e)>0 ) h_Enu_MC_TReco          -> Fill( CC1piAnalysisHelper::GetMCCC1piNeutrinoEnergy(e), (CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e)-CC1piAnalysisHelper::GetMCCC1piNeutrinoEnergy(e))/CC1piAnalysisHelper::GetMCCC1piNeutrinoEnergy(e) );
	    }
    }
    if( e.CheckRecoTopology( topology ) == 1 ) {
      if( topology == GeneralAnalysisHelper::GetCC0PiTopologyMap() && CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e)>0 ) h_Reco_energy_nu_cc0pi -> Fill( CC0piAnalysisHelper::GetRecoCC0piNeutrinoEnergy(e) );
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e)>0 ) h_Reco_energy_nu_cc1pi -> Fill( CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergy(e) );
    }

    // ------------------------------------------  Neutrino Energy method 2  ----------------------------------------------------
    if( e.CheckMCTopology( topology ) == 1 ) { 
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() ) h_MC_energy_nu_cc1pi_M2 -> Fill( e.GetTrueNuEnergy() );  
      if( e.CheckRecoTopology( topology ) == 1 ) { 
	      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergyMethod2(e)>0 ) h_TReco_energy_nu_cc1pi_M2 -> Fill( CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergyMethod2(e) );
	    }
    }
    if( e.CheckRecoTopology( topology ) == 1 ) {
      if( topology == GeneralAnalysisHelper::GetCC1PiTopologyMap() && CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergyMethod2(e)>0 ) h_Reco_energy_nu_cc1pi_M2 -> Fill( CC1piAnalysisHelper::GetRecoCC1piNeutrinoEnergyMethod2(e) );
    }
  } //endfor
  
  std::cout<<"END LOOP"<<std::endl;
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;

  // Plot histograms - Exemple CC1pi+/-
  TCanvas * c = new TCanvas( "c_mu_MC_L", "c_mu_MC_L", 600, 600 );
  TLegend *leg = new TLegend(0.9,0.7,0.48,0.9) ;
 
  
  //----------------HISTOGRAM 1----------------------- 
  h_mu_Reco_L->GetYaxis()->SetTitle("Events [#]");
  h_mu_Reco_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_Reco_L->SetLineColor(kRed);
  h_mu_Reco_L->SetStats(0);
  h_mu_Reco_L->GetYaxis()->SetTitleSize(20);
  h_mu_Reco_L->GetYaxis()->SetTitleFont(43);
  h_mu_Reco_L->GetYaxis()->SetLabelFont(43);
  h_mu_Reco_L->GetYaxis()->SetLabelSize(15);
  h_mu_Reco_L->GetXaxis()->SetTitleSize(20);
  h_mu_Reco_L->GetXaxis()->SetTitleFont(43);
  h_mu_Reco_L->GetXaxis()->SetLabelFont(43);
  h_mu_Reco_L->GetXaxis()->SetLabelSize(15);
  h_mu_Reco_L->GetXaxis()->SetRangeUser(0.,500.);
  h_mu_Reco_L->GetYaxis()->SetRangeUser(0.,200.);
  h_mu_Reco_L->SetTitle("") ;
  h_mu_Reco_L-> Draw();

  h_mu_MC_L->GetYaxis()->SetTitle("Event Rate [a.u]");
  h_mu_MC_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_MC_L->SetLineColor(kBlack);
  h_mu_MC_L-> SetStats(0);
  h_mu_MC_L->SetLineWidth(1);
  h_mu_MC_L->SetTitle("") ;
  leg -> AddEntry( h_mu_MC_L , " True: <L>=113.5" );
  h_mu_MC_L-> Draw("same");
  
  leg -> AddEntry( h_mu_Reco_L , " Selected: <L>=116.0 " ); 
  h_mu_TReco_L->GetYaxis()->SetTitle("Eff[%]");
  h_mu_TReco_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_TReco_L->SetFillColor( 8 );
  h_mu_TReco_L->SetFillStyle(3144);
  h_mu_TReco_L-> SetStats(0);
  leg -> AddEntry( h_mu_TReco_L , " Signal: <L>=116.4 " );
  h_mu_TReco_L->GetXaxis()->SetRangeUser(0.,500.);
  h_mu_TReco_L->GetYaxis()->SetRangeUser(0.,200.);
  h_mu_TReco_L->SetTitle("") ;
  h_mu_TReco_L-> Draw("same");
  TLatex *latex= new TLatex();
  leg -> Draw();
  c->Print( (filepath+"h_mu_ALL_L.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_ALL_L.root").c_str() );
  
  c->Clear();
  leg->Clear();
  //-----------------------------HISTOGRAM 2--------------------------------

  h_pi_MC_L->GetYaxis()->SetTitle("Events [#]");
  h_pi_MC_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_MC_L->SetLineColor(kBlack);
  h_pi_MC_L-> SetStats(0);
  h_pi_MC_L->SetLineWidth(1);
  h_pi_MC_L->SetTitle("") ;
  leg -> AddEntry( h_pi_MC_L , " True: <L>=28.1 " );
  h_pi_MC_L->GetYaxis()->SetRangeUser(0.,350.);
  h_pi_MC_L-> Draw();
  
  h_pi_Reco_L->GetYaxis()->SetTitle("Events [#]");
  h_pi_Reco_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_Reco_L->SetLineColor(kRed);
  h_pi_Reco_L-> SetStats(0);                                                                         
  h_pi_Reco_L->GetYaxis()->SetTitleSize(20);
  h_pi_Reco_L->GetYaxis()->SetTitleFont(43);
  h_pi_Reco_L->GetYaxis()->SetLabelFont(43);
  h_pi_Reco_L->GetYaxis()->SetLabelSize(15);
  h_pi_Reco_L->GetXaxis()->SetTitleSize(20);
  h_pi_Reco_L->GetXaxis()->SetTitleFont(43);
  h_pi_Reco_L->GetXaxis()->SetLabelFont(43);
  h_pi_Reco_L->GetXaxis()->SetLabelSize(15);
  h_pi_Reco_L->GetXaxis()->SetRangeUser(0.,350.);
  h_pi_Reco_L->SetTitle("");
  h_pi_Reco_L->GetYaxis()->SetRangeUser(0.,350.);
  h_pi_Reco_L-> Draw("same");

  leg -> AddEntry( h_pi_Reco_L , " Selected: <L>=24.4 " );
  h_pi_TReco_L->GetYaxis()->SetTitle("Events [#]");
  h_pi_TReco_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_TReco_L->SetFillColor( 8 );
  h_pi_TReco_L->SetFillStyle(3144);
  h_pi_TReco_L-> SetStats(0);
  leg -> AddEntry( h_pi_TReco_L , " Signal: <L>=30.7 " );
  h_pi_TReco_L->GetXaxis()->SetRangeUser(0.,500.);
  h_pi_TReco_L->GetYaxis()->SetRangeUser(0.,350.);
  h_pi_TReco_L->SetTitle("") ;
  h_pi_TReco_L-> Draw("same");
  leg -> Draw();

  c->Print( (filepath+"h_pi_ALL_L.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_ALL_L.root").c_str() );
  c->Clear();
  leg->Clear();

  //------------------------------HISTO 3
  h_mu_Reco_T->GetYaxis()->SetTitle("Events [#]");
  h_mu_Reco_T->GetXaxis()->SetTitle("cos(theta)");
  h_mu_Reco_T->SetLineColor(kRed);
  h_mu_Reco_T-> SetStats(0);                                                             
  h_mu_Reco_T->GetYaxis()->SetTitleSize(20);
  h_mu_Reco_T->GetYaxis()->SetTitleFont(43);
  h_mu_Reco_T->GetYaxis()->SetLabelFont(43);
  h_mu_Reco_T->GetYaxis()->SetLabelSize(15);
  h_mu_Reco_T->GetXaxis()->SetTitleSize(20);
  h_mu_Reco_T->GetXaxis()->SetTitleFont(43);
  h_mu_Reco_T->GetXaxis()->SetLabelFont(43);
  h_mu_Reco_T->GetXaxis()->SetLabelSize(15);
  h_mu_Reco_T->SetTitle("");
  h_mu_Reco_T->GetYaxis()->SetRangeUser(0.,250.);
  h_mu_Reco_T-> Draw();
  h_mu_MC_T->GetYaxis()->SetTitle("Events [#]");
  h_mu_MC_T->GetXaxis()->SetTitle("cos(theta)");
  h_mu_MC_T->SetLineColor(kBlack);
  h_mu_MC_T-> SetStats(0);
  h_mu_MC_T->SetLineWidth(1);
  h_mu_MC_T->SetTitle("") ;
  leg -> AddEntry( h_mu_MC_T , " True " );
  h_mu_MC_T-> Draw("same");
 
  leg -> AddEntry( h_mu_Reco_T , " Selected " );
  h_mu_TReco_T->GetYaxis()->SetTitle("Events [#]");
  h_mu_TReco_T->GetXaxis()->SetTitle("cos(\theta)");
  h_mu_TReco_T->SetFillColor( 8 );
  h_mu_TReco_T->SetFillStyle(3144);
  h_mu_TReco_T-> SetStats(0);
  leg -> AddEntry( h_mu_TReco_T , " Signal " );
  h_mu_TReco_T->SetTitle("") ;
  h_mu_TReco_T-> Draw("same");
  leg -> Draw();
  c->Print( (filepath+"h_mu_ALL_T.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_ALL_T.root").c_str() );
  c->Clear();
  leg->Clear();

  //-------------------HISTO 4
  h_pi_MC_T->GetYaxis()->SetTitle("Events [#]");
  h_pi_MC_T->GetXaxis()->SetTitle("cos(theta)");
  h_pi_MC_T->GetYaxis()->SetRangeUser(0.,250.);
  h_pi_MC_T->SetLineColor(kBlack);
  h_pi_MC_T-> SetStats(0);
  h_pi_MC_T->GetYaxis()->SetRangeUser(0.,120.);
  h_pi_MC_T->SetLineWidth(1);
  h_pi_MC_T->SetTitle("") ;
  leg -> AddEntry( h_pi_MC_T , " True " );
  h_pi_MC_T-> Draw();
  
  h_pi_Reco_T->GetYaxis()->SetTitle("Events [#]");
  h_pi_Reco_T->GetXaxis()->SetTitle("cos(theta)");
  h_pi_Reco_T->SetLineColor(kRed);
  h_pi_Reco_T-> SetStats(0);
  h_pi_Reco_T->GetYaxis()->SetTitleSize(20);
  h_pi_Reco_T->GetYaxis()->SetTitleFont(43);
  h_pi_Reco_T->GetYaxis()->SetLabelFont(43);
  h_pi_Reco_T->GetYaxis()->SetLabelSize(15);
  h_pi_Reco_T->GetXaxis()->SetTitleSize(20);
  h_pi_Reco_T->GetXaxis()->SetTitleFont(43);
  h_pi_Reco_T->GetXaxis()->SetLabelFont(43);
  h_pi_Reco_T->GetXaxis()->SetLabelSize(15);
  h_pi_Reco_T->SetTitle("");
  h_pi_Reco_T-> Draw("same");
  leg -> AddEntry( h_pi_Reco_T , " Selected" );

  h_pi_TReco_T->GetXaxis()->SetTitle("cos(theta)");
  h_pi_TReco_T->SetFillColor( 8 );
  h_pi_TReco_T->SetFillStyle(3144);
  h_pi_TReco_T-> SetStats(0);
  leg -> AddEntry( h_pi_TReco_T , " Signal " );
  h_pi_TReco_T->SetTitle("") ;
  h_pi_TReco_T-> Draw("same");
  leg -> Draw();
  c->Print( (filepath+"h_pi_ALL_T.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_ALL_T.root").c_str() );
  c->Clear();
  leg->Clear();
  //---------------------------HISTO 5  
  //---------------HISTO ENERGY
  h_Reco_E_mu->SetLineColor(kRed);
  h_Reco_E_mu-> SetStats(0);
  h_Reco_E_mu->GetYaxis()->SetRangeUser(0.,350.);
  h_Reco_E_mu->SetTitle("");
  h_Reco_E_mu->GetXaxis()->SetTitle("T/GeV");
  h_Reco_E_mu->Draw( );
  h_MC_E_mu->GetYaxis()->SetTitle("Events [#]");
  h_MC_E_mu->SetLineColor(kBlack);
  h_MC_E_mu-> SetStats(0);
  h_MC_E_mu->GetXaxis()->SetRangeUser(0.,3.);
  h_MC_E_mu->SetLineWidth(1);
  h_MC_E_mu->SetTitle("") ;
  h_MC_E_mu->GetXaxis()->SetTitle("T/GeV");
  leg -> AddEntry( h_MC_E_mu , " True: <T>=0.48 " );
  h_MC_E_mu->GetYaxis()->SetRangeUser(0.,350.);
  h_MC_E_mu->Draw("same");
  leg -> AddEntry( h_Reco_E_mu , " Selected: <T>=0.34 " );
  h_TReco_E_mu->SetFillColor( 8 );
  h_TReco_E_mu->SetFillStyle(3144);
  h_TReco_E_mu-> SetStats(0);
  h_TReco_E_mu->GetYaxis()->SetTitleSize(20);
  h_TReco_E_mu->GetYaxis()->SetTitleFont(43);
  h_TReco_E_mu->GetYaxis()->SetLabelFont(43);
  h_TReco_E_mu->GetYaxis()->SetLabelSize(15);
  h_TReco_E_mu->GetXaxis()->SetTitleSize(20);
  h_TReco_E_mu->GetXaxis()->SetTitleFont(43);
  h_TReco_E_mu->GetXaxis()->SetLabelFont(43);
  h_TReco_E_mu->GetXaxis()->SetLabelSize(15);
  h_TReco_E_mu->GetYaxis()->SetRangeUser(0.,350.);
  h_TReco_E_mu->SetTitle("");
  h_TReco_E_mu->GetXaxis()->SetTitle("T/GeV");
  h_TReco_E_mu->GetYaxis()->SetTitle("Event Rate [a.u]");
  h_TReco_E_mu->Draw("same");
  leg -> AddEntry( h_TReco_E_mu , " Signal: <T>=0.38 " );
  leg->Draw();
  c->Print( (filepath+"h_ALL_E_mu.pdf").c_str() );
  c->SaveAs( (filepath+"h_ALL_E_mu.root").c_str() );
  c->Clear();
  leg->Clear();
  // ------------------_HISTOGRAM TENERGY pi
  h_Reco_E_pi->SetLineColor(kRed);
  h_Reco_E_pi-> SetStats(0);
  h_Reco_E_pi->SetTitle("");
  h_Reco_E_pi->GetXaxis()->SetTitle("T/GeV");
  h_Reco_E_pi->Draw();
  h_MC_E_pi->GetYaxis()->SetTitle("Events [#]");
  h_MC_E_pi->GetXaxis()->SetTitle("T/GeV");
  h_MC_E_pi->SetLineColor(kBlack);
  h_MC_E_pi-> SetStats(0);
  h_MC_E_pi->SetLineWidth(1);
  h_MC_E_pi->GetXaxis()->SetRangeUser(0.,3.);
  h_MC_E_pi->SetTitle("") ;
  leg -> AddEntry( h_MC_E_pi , " True: <T>=0.22 " );
  h_MC_E_pi->Draw("same");
  leg -> AddEntry( h_Reco_E_pi , " Selected: <T>=0.34 " );
  h_TReco_E_pi->SetFillColor( 8 );
  h_TReco_E_pi->SetFillStyle(3144);
  h_TReco_E_pi-> SetStats(0);
  h_TReco_E_pi->GetYaxis()->SetTitleSize(20);
  h_TReco_E_pi->GetYaxis()->SetTitleFont(43);
  h_TReco_E_pi->GetYaxis()->SetLabelFont(43);
  h_TReco_E_pi->GetYaxis()->SetLabelSize(15);
  h_TReco_E_pi->GetXaxis()->SetTitleSize(20);
  h_TReco_E_pi->GetXaxis()->SetTitleFont(43);
  h_TReco_E_pi->GetXaxis()->SetLabelFont(43);
  h_TReco_E_pi->GetXaxis()->SetLabelSize(15);
  h_TReco_E_pi->SetTitle("");
  h_TReco_E_pi->GetXaxis()->SetTitle("T/GeV");
  h_TReco_E_pi->Draw("same");
  leg->Draw();
  leg -> AddEntry( h_TReco_E_pi , " Signal: <T>=0.38 " );
  c->Print( (filepath+"h_ALL_E_pi.pdf").c_str() );
  c->SaveAs( (filepath+"h_ALL_E_pi.root").c_str() );
  c->Clear();
  leg->Clear();
  // ------------------_HISTOGRAM TENERGY mu
  h_R_E_mu->SetTitle("");
  h_R_E_mu->GetYaxis()->SetTitle("Events [#]");
  h_R_E_mu->GetXaxis()->SetTitle("T#_#mu((signal-true)/true)");
  h_R_E_mu->Draw();
  
  c->Print( (filepath+"h_mu_ratio.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_ratio.root").c_str() );
  c->Clear();
  leg->Clear();

  // ------------------_HISTOGRAM TENERGY pi
  h_R_E_pi->SetTitle("");
  h_R_E_pi->GetYaxis()->SetTitle("Events [#]");
  h_R_E_pi->GetXaxis()->SetTitle("T#_#pi((signal-true)/true)");
  h_R_E_pi->Draw();
  
  c->Print( (filepath+"h_pi_ratio.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_ratio.root").c_str() );
  c->Clear();
  leg->Clear();

  // ------------------_HISTOGRAM TENERGY pi
  h_Enu_MC_TReco->SetTitle("");
  h_Enu_MC_TReco->GetXaxis()->SetTitle("E#_#nu [GeV]");
  h_Enu_MC_TReco->GetYaxis()->SetTitle("E#_#nu((signal-true)/true)");
  h_Enu_MC_TReco->Draw("COL");
  
  c->Print( (filepath+"h_Enu_MC_TReco.pdf").c_str() );
  c->SaveAs( (filepath+"h_Enu_MC_TReco.root").c_str() );
  c->Clear();
  leg->Clear();
  // ------------------_HISTOGRAM Q2  cc1pi
  h_Q_MC_TReco->SetTitle("");
  h_Q_MC_TReco->GetXaxis()->SetTitle("Q^2 [GeV^2/c^4]");
  h_Q_MC_TReco->GetYaxis()->SetTitle("Q^2((signal-true)/true)");
  h_Q_MC_TReco->Draw("COL");
  
  c->Print( (filepath+"h_Q_MC_TReco.pdf").c_str() );
  c->SaveAs( (filepath+"h_Q_MC_TReco.root").c_str() );
  c->Clear();
  leg->Clear();


  //----------------------------DELTA MASS HISTOGRAM
  TF1 * bw               = new TF1  ( "Bw"              , "TMath::BreitWigner(x, 1.232, 0.117 )", 1., 1.5 );
  TLegend *leg6 = new TLegend(0.9,0.9,0.48,0.6) ;
  h_MC_Delta->GetYaxis()->SetTitleSize(20);
  h_MC_Delta->GetYaxis()->SetTitleFont(43);
  h_MC_Delta->GetYaxis()->SetLabelFont(43);
  h_MC_Delta->GetYaxis()->SetLabelSize(15);
  h_MC_Delta->GetXaxis()->SetTitleSize(20);
  h_MC_Delta->GetXaxis()->SetTitleFont(43);
  h_MC_Delta->GetXaxis()->SetLabelFont(43);
  h_MC_Delta->GetXaxis()->SetLabelSize(15);
  h_MC_Delta-> SetStats(0);
  h_MC_Delta->GetYaxis()->SetTitle("Event Rate [a.u]");
  h_MC_Delta->SetFillColor( 39 );
  h_MC_Delta->SetFillStyle(3004);
  h_MC_Delta->SetTitle("");
  h_MC_Delta->GetXaxis()->SetTitle("Invariant Mass/GeV");
  h_MC_Delta->Scale(100/CountMC[3]);
  bw->SetLineStyle( 9 );
  h_MC_Delta->Draw();
  h_TReco_Delta->GetYaxis()->SetTitle("Event Rate [a.u]");
  h_TReco_Delta->GetXaxis()->SetTitle("Invariant mass/GeV");
  h_TReco_Delta->SetFillColor( 8 );
  h_TReco_Delta->SetFillStyle(3144);
  h_TReco_Delta-> SetStats(0);
  h_TReco_Delta->GetYaxis()->SetTitleSize(20);
  h_TReco_Delta->GetYaxis()->SetTitleFont(43);
  h_TReco_Delta->GetYaxis()->SetLabelFont(43);
  h_TReco_Delta->GetYaxis()->SetLabelSize(15);
  h_TReco_Delta->GetXaxis()->SetTitleSize(20);
  h_TReco_Delta->GetXaxis()->SetTitleFont(43);
  h_TReco_Delta->GetXaxis()->SetLabelFont(43);
  h_TReco_Delta->GetXaxis()->SetLabelSize(15);
  h_TReco_Delta->GetYaxis()->SetRangeUser(0.,150.);
  h_TReco_Delta->Scale(100/CountMC[3]);
  h_TReco_Delta ->Draw("same");
  leg6 -> AddEntry( h_MC_Delta  , " True " );
  leg6 -> AddEntry( h_TReco_Delta, " Signal " );
  leg6 -> AddEntry( bw          , " BW " );
  bw->Draw("same");
  leg6->Draw();
  c->Print( (filepath+"h_All_Delta.pdf").c_str() );
  c->SaveAs( (filepath+"h_All_Delta.root").c_str() );
  c->Clear();

  // ----------- Stack Histograms : Background study ( mu ) ------------
  TLegend *leg1 = new TLegend(0.9,0.7,0.48,0.9);
  h_mu_Reco_cc1pi_L->SetFillColor(kGreen);
  h_mu_Reco_cc0pi_L->SetFillColor(kRed);
  h_mu_Reco_ccpi0_L->SetFillColor(kBlue);
  h_mu_Reco_nc_L->SetFillColor(kOrange);
  h_mu_Reco_cc0pi_L->SetFillStyle(3004);
  h_mu_Reco_ccpi0_L->SetFillStyle(3004);
  h_mu_Reco_nc_L->SetFillStyle(3004);
  
  hS_mu_Reco_L->Add(h_mu_Reco_cc1pi_L);
  hS_mu_Reco_L->Add(h_mu_Reco_cc0pi_L);
  hS_mu_Reco_L->Add(h_mu_Reco_ccpi0_L);
  hS_mu_Reco_L->Add(h_mu_Reco_nc_L);

  leg1->AddEntry(h_mu_Reco_cc1pi_L, "Signal");
  leg1->AddEntry(h_mu_Reco_cc0pi_L, "BG: CC0#pi");
  leg1->AddEntry(h_mu_Reco_ccpi0_L, "BG: CC1#pi0");
  leg1->AddEntry(h_mu_Reco_nc_L,    "BG: other");

  TLegend *leg22 = new TLegend(0.9,0.7,0.48,0.9);
  h_mu_Reco_cc1pi_T->SetFillColor(kGreen);
  h_mu_Reco_cc0pi_T->SetFillColor(kRed);
  h_mu_Reco_ccpi0_T->SetFillColor(kBlue);
  h_mu_Reco_nc_T->SetFillColor(kOrange);

  h_mu_Reco_cc0pi_T->SetFillStyle(3004);
  h_mu_Reco_ccpi0_T->SetFillStyle(3004);
  h_mu_Reco_nc_T->SetFillStyle(3004);
  
  hS_mu_Reco_T->Add(h_mu_Reco_cc1pi_T);
  hS_mu_Reco_T->Add(h_mu_Reco_cc0pi_T);
  hS_mu_Reco_T->Add(h_mu_Reco_ccpi0_T);
  hS_mu_Reco_T->Add(h_mu_Reco_nc_T);
 

  leg22->AddEntry(h_mu_Reco_cc1pi_T, "Signal");
  leg22->AddEntry(h_mu_Reco_cc0pi_T, "BG: CC0#pi");
  leg22->AddEntry(h_mu_Reco_ccpi0_T, "BG: CC#pi0");
  leg22->AddEntry(h_mu_Reco_nc_T,    "BG: other");

  TLegend *leg33 = new TLegend(0.9,0.7,0.48,0.9);
  h_mu_Reco_cc1pi_E->SetFillColor(kGreen);
  h_mu_Reco_cc0pi_E->SetFillColor(kRed);
  h_mu_Reco_ccpi0_E->SetFillColor(kBlue);
  h_mu_Reco_nc_E->SetFillColor(kOrange);
  

  h_mu_Reco_cc0pi_E->SetFillStyle(3004);
  h_mu_Reco_ccpi0_E->SetFillStyle(3004);
  h_mu_Reco_nc_E->SetFillStyle(3004);
  hS_mu_Reco_E->Add(h_mu_Reco_cc1pi_E);
  hS_mu_Reco_E->Add(h_mu_Reco_cc0pi_E);
  hS_mu_Reco_E->Add(h_mu_Reco_ccpi0_E);
  hS_mu_Reco_E->Add(h_mu_Reco_nc_E);
 
  leg33->AddEntry(h_mu_Reco_cc1pi_E, "Signal");
  leg33->AddEntry(h_mu_Reco_cc0pi_E, "BG: CC0#pi");
  leg33->AddEntry(h_mu_Reco_ccpi0_E, "BG: CC#pi0");
  leg33->AddEntry(h_mu_Reco_nc_E,    "BG: other");

  hS_mu_Reco_L->SetTitle("");
  hS_mu_Reco_T->SetTitle("");
  hS_mu_Reco_E->SetTitle("");

  
  hS_mu_Reco_L->Draw();
  leg1->Draw();
  c->Print( (filepath+"h_mu_Reco_L.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_Reco_L.root").c_str() );
  c->Clear();
  hS_mu_Reco_T->Draw();
  leg22->Draw();
  c->Print( (filepath+"h_mu_Reco_T.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_Reco_T.root").c_str() );
  c->Clear();
  hS_mu_Reco_E->Draw();
  leg33->Draw();
  c->Print( (filepath+"h_mu_Reco_E.pdf").c_str() );
  c->SaveAs( (filepath+"h_mu_Reco_E.root").c_str() );
  c->Clear();

   // ----------- Stack Histograms : Background study ( pi ) ------------
  TLegend *leg111 = new TLegend(0.9,0.7,0.48,0.9);
  h_pi_Reco_cc1pi_L->SetFillColor(kGreen);
  h_pi_Reco_cc0pi_L->SetFillColor(kRed);
  h_pi_Reco_ccpi0_L->SetFillColor(kBlue);
  h_pi_Reco_nc_L->SetFillColor(kOrange);
  h_pi_Reco_cc0pi_L->SetFillStyle(3004);
  h_pi_Reco_ccpi0_L->SetFillStyle(3004);
  h_pi_Reco_nc_L->SetFillStyle(3004);
  
  hS_pi_Reco_L->Add(h_pi_Reco_cc1pi_L);
  hS_pi_Reco_L->Add(h_pi_Reco_cc0pi_L);
  hS_pi_Reco_L->Add(h_pi_Reco_ccpi0_L);
  hS_pi_Reco_L->Add(h_pi_Reco_nc_L);

  leg111->AddEntry(h_pi_Reco_cc1pi_L, "Signal");
  leg111->AddEntry(h_pi_Reco_cc0pi_L, "BG: CC0#pi");
  leg111->AddEntry(h_pi_Reco_ccpi0_L, "BG: CC1#pi0");
  leg111->AddEntry(h_pi_Reco_nc_L,    "BG: other");

  TLegend *leg222 = new TLegend(0.9,0.7,0.48,0.9);
  h_pi_Reco_cc1pi_T->SetFillColor(kGreen);
  h_pi_Reco_cc0pi_T->SetFillColor(kRed);
  h_pi_Reco_ccpi0_T->SetFillColor(kBlue);
  h_pi_Reco_nc_T->SetFillColor(kOrange);

  h_pi_Reco_cc0pi_T->SetFillStyle(3004);
  h_pi_Reco_ccpi0_T->SetFillStyle(3004);
  h_pi_Reco_nc_T->SetFillStyle(3004);
  
  hS_pi_Reco_T->Add(h_pi_Reco_cc1pi_T);
  hS_pi_Reco_T->Add(h_pi_Reco_cc0pi_T);
  hS_pi_Reco_T->Add(h_pi_Reco_ccpi0_T);
  hS_pi_Reco_T->Add(h_pi_Reco_nc_T);
 

  leg222->AddEntry(h_pi_Reco_cc1pi_T, "Signal");
  leg222->AddEntry(h_pi_Reco_cc0pi_T, "BG: CC0#pi");
  leg222->AddEntry(h_pi_Reco_ccpi0_T, "BG: CC#pi0");
  leg222->AddEntry(h_pi_Reco_nc_T,    "BG: other");

  TLegend *leg333 = new TLegend(0.9,0.7,0.48,0.9);
  h_pi_Reco_cc1pi_E->SetFillColor(kGreen);
  h_pi_Reco_cc0pi_E->SetFillColor(kRed);
  h_pi_Reco_ccpi0_E->SetFillColor(kBlue);
  h_pi_Reco_nc_E->SetFillColor(kOrange);
  

  h_pi_Reco_cc0pi_E->SetFillStyle(3004);
  h_pi_Reco_ccpi0_E->SetFillStyle(3004);
  h_pi_Reco_nc_E->SetFillStyle(3004);
  hS_pi_Reco_E->Add(h_pi_Reco_cc1pi_E);
  hS_pi_Reco_E->Add(h_pi_Reco_cc0pi_E);
  hS_pi_Reco_E->Add(h_pi_Reco_ccpi0_E);
  hS_pi_Reco_E->Add(h_pi_Reco_nc_E);
 
  leg333->AddEntry(h_pi_Reco_cc1pi_E, "Signal");
  leg333->AddEntry(h_pi_Reco_cc0pi_E, "BG: CC0#pi");
  leg333->AddEntry(h_pi_Reco_ccpi0_E, "BG: CC#pi0");
  leg333->AddEntry(h_pi_Reco_nc_E,    "BG: other");

  hS_pi_Reco_L->SetTitle("");
  hS_pi_Reco_T->SetTitle("");
  hS_pi_Reco_E->SetTitle("");
  hS_pi_Reco_L->Draw();
  leg111->Draw();
  c->Print( (filepath+"h_pi_Reco_L.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_Reco_L.root").c_str() );
  c->Clear();
  hS_pi_Reco_T->Draw();
  leg222->Draw();
  c->Print( (filepath+"h_pi_Reco_T.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_Reco_T.root").c_str() );
  c->Clear();
  hS_pi_Reco_E->Draw();
  leg33->Draw();
  c->Print( (filepath+"h_pi_Reco_E.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi_Reco_E.root").c_str() );
  c->Clear();

   // ----------- Stack Histograms : Background study ( pi0 ) ------------
  TLegend *leg1111 = new TLegend(0.9,0.7,0.48,0.9);
  h_pi0_Reco_cc1pi_L->SetFillColor(kGreen);
  h_pi0_Reco_cc0pi_L->SetFillColor(kRed);
  h_pi0_Reco_ccpi0_L->SetFillColor(kBlue);
  h_pi0_Reco_nc_L->SetFillColor(kOrange);
  h_pi0_Reco_cc0pi_L->SetFillStyle(3004);
  h_pi0_Reco_ccpi0_L->SetFillStyle(3004);
  h_pi0_Reco_nc_L->SetFillStyle(3004);
  
  hS_pi0_Reco_L->Add(h_pi0_Reco_cc1pi_L);
  hS_pi0_Reco_L->Add(h_pi0_Reco_cc0pi_L);
  hS_pi0_Reco_L->Add(h_pi0_Reco_ccpi0_L);
  hS_pi0_Reco_L->Add(h_pi0_Reco_nc_L);

  leg1111->AddEntry(h_pi0_Reco_cc1pi_L, "Signal");
  leg1111->AddEntry(h_pi0_Reco_cc0pi_L, "BG: CC0#pi");
  leg1111->AddEntry(h_pi0_Reco_ccpi0_L, "BG: CC1#pi0");
  leg1111->AddEntry(h_pi0_Reco_nc_L,    "BG: other");

  TLegend *leg2222 = new TLegend(0.9,0.7,0.48,0.9);
  h_pi0_Reco_cc1pi_T->SetFillColor(kGreen);
  h_pi0_Reco_cc0pi_T->SetFillColor(kRed);
  h_pi0_Reco_ccpi0_T->SetFillColor(kBlue);
  h_pi0_Reco_nc_T->SetFillColor(kOrange);

  h_pi0_Reco_cc0pi_T->SetFillStyle(3004);
  h_pi0_Reco_ccpi0_T->SetFillStyle(3004);
  h_pi0_Reco_nc_T->SetFillStyle(3004);
  
  hS_pi0_Reco_T->Add(h_pi0_Reco_cc1pi_T);
  hS_pi0_Reco_T->Add(h_pi0_Reco_cc0pi_T);
  hS_pi0_Reco_T->Add(h_pi0_Reco_ccpi0_T);
  hS_pi0_Reco_T->Add(h_pi0_Reco_nc_T);
 

  leg2222->AddEntry(h_pi0_Reco_cc1pi_T, "Signal");
  leg2222->AddEntry(h_pi0_Reco_cc0pi_T, "BG: CC0#pi");
  leg2222->AddEntry(h_pi0_Reco_ccpi0_T, "BG: CC#pi0");
  leg2222->AddEntry(h_pi0_Reco_nc_T,    "BG: other");

  hS_pi0_Reco_L->SetTitle("");
  hS_pi0_Reco_T->SetTitle("");
  hS_pi0_Reco_L->Draw();
  leg111->Draw();
  c->Print( (filepath+"h_pi0_Reco_L.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi0_Reco_L.root").c_str() );
  c->Clear();
  hS_pi0_Reco_T->Draw();
  leg222->Draw();
  c->Print( (filepath+"h_pi0_Reco_T.pdf").c_str() );
  c->SaveAs( (filepath+"h_pi0_Reco_T.root").c_str() );
  c->Clear();
 
  //--------------Neutrino Energy CC0pi
  h_Reco_energy_nu_cc0pi->GetYaxis()->SetTitle("Events [#]");
  h_Reco_energy_nu_cc0pi->GetXaxis()->SetTitle("E#nu[GeV]");
  h_Reco_energy_nu_cc0pi->SetFillColor(kRed-3);
  h_Reco_energy_nu_cc0pi->SetFillStyle(3004);
  h_Reco_energy_nu_cc0pi->SetStats(0);
  h_Reco_energy_nu_cc0pi->GetYaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc0pi->GetYaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc0pi->GetYaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc0pi->GetYaxis()->SetLabelSize(15);
  h_Reco_energy_nu_cc0pi->GetXaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc0pi->GetXaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc0pi->GetXaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc0pi->GetXaxis()->SetLabelSize(15);
  h_Reco_energy_nu_cc0pi->SetTitle("") ;
  h_Reco_energy_nu_cc0pi-> Draw();

  h_MC_energy_nu_cc0pi->GetYaxis()->SetTitle("Events [#]");
  h_MC_energy_nu_cc0pi->GetXaxis()->SetTitle("E#nu[GeV]");
  h_MC_energy_nu_cc0pi->SetLineColor(kBlue);
  h_MC_energy_nu_cc0pi->SetFillStyle(3004);
  h_MC_energy_nu_cc0pi-> SetStats(0);
  h_MC_energy_nu_cc0pi->SetLineWidth(1);
  h_MC_energy_nu_cc0pi->SetTitle("") ;
  leg -> AddEntry( h_MC_energy_nu_cc0pi , " True" );
  h_MC_energy_nu_cc0pi -> Draw("same");
  leg -> AddEntry( h_Reco_energy_nu_cc0pi , " Selected " ); 
  
  h_TReco_energy_nu_cc0pi->SetFillColor( 8 );
  h_TReco_energy_nu_cc0pi->SetFillStyle(3144);
  h_TReco_energy_nu_cc0pi-> SetStats(0);
  leg -> AddEntry( h_TReco_energy_nu_cc0pi , " Signal " );
  h_TReco_energy_nu_cc0pi->SetTitle("") ;
  h_TReco_energy_nu_cc0pi-> Draw("same");
  leg -> Draw();
  c->Print( (filepath+"h_energy_nu_cc0pi.pdf").c_str() );
  c->SaveAs( (filepath+"h_energy_nu_cc0pi.root").c_str() );
  
  c->Clear();
  leg->Clear();

 //--------------Neutrino Energy CC1pi
  h_Reco_energy_nu_cc1pi->GetYaxis()->SetTitle("Events [#]");
  h_Reco_energy_nu_cc1pi->GetXaxis()->SetTitle("E#nu[GeV]");
  h_Reco_energy_nu_cc1pi->SetFillColor(kRed-3);
  h_Reco_energy_nu_cc1pi->SetFillStyle(3004);
  h_Reco_energy_nu_cc1pi->SetStats(0);
  h_Reco_energy_nu_cc1pi->GetYaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc1pi->GetYaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc1pi->GetYaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc1pi->GetYaxis()->SetLabelSize(15);
  h_Reco_energy_nu_cc1pi->GetXaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc1pi->GetXaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc1pi->GetXaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc1pi->GetXaxis()->SetLabelSize(15);
  h_Reco_energy_nu_cc1pi->SetTitle("") ;
  h_Reco_energy_nu_cc1pi-> Draw();

  h_MC_energy_nu_cc1pi->GetYaxis()->SetTitle("Events [#]");
  h_MC_energy_nu_cc1pi->GetXaxis()->SetTitle("E#nu[GeV]");
  h_MC_energy_nu_cc1pi->SetFillColor(kCyan-6);
  h_MC_energy_nu_cc1pi->SetFillStyle(3004);
  h_MC_energy_nu_cc1pi-> SetStats(0);
  h_MC_energy_nu_cc1pi->SetLineWidth(1);
  h_MC_energy_nu_cc1pi->SetTitle("") ;
  leg -> AddEntry( h_MC_energy_nu_cc1pi , " True" );
  h_MC_energy_nu_cc1pi -> Draw("same");
  leg -> AddEntry( h_Reco_energy_nu_cc1pi , " Selected " ); 
  
  h_TReco_energy_nu_cc1pi->SetFillColor( 8 );
  h_TReco_energy_nu_cc1pi->SetFillStyle(3144);
  h_TReco_energy_nu_cc1pi-> SetStats(0);
  leg -> AddEntry( h_TReco_energy_nu_cc1pi , " Signal " );
  h_TReco_energy_nu_cc1pi->SetTitle("") ;
  h_TReco_energy_nu_cc1pi-> Draw("same");
  leg -> Draw();
  c->Print( (filepath+"h_energy_nu_cc1pi.pdf").c_str() );
  c->SaveAs( (filepath+"h_energy_nu_cc1pi.root").c_str() );
  
  c->Clear();
  leg->Clear();

 //--------------Neutrino Energy CC1pi Method 2
  h_Reco_energy_nu_cc1pi_M2->GetYaxis()->SetTitle("Events [#]");
  h_Reco_energy_nu_cc1pi_M2->GetXaxis()->SetTitle("E#nu[GeV]");
  h_Reco_energy_nu_cc1pi_M2->SetFillColor(kRed-3);
  h_Reco_energy_nu_cc1pi_M2->SetFillStyle(3004);
  h_Reco_energy_nu_cc1pi_M2->SetStats(0);
  h_Reco_energy_nu_cc1pi_M2->GetYaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc1pi_M2->GetYaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc1pi_M2->GetYaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc1pi_M2->GetYaxis()->SetLabelSize(15);
  h_Reco_energy_nu_cc1pi_M2->GetXaxis()->SetTitleSize(20);
  h_Reco_energy_nu_cc1pi_M2->GetXaxis()->SetTitleFont(43);
  h_Reco_energy_nu_cc1pi_M2->GetXaxis()->SetLabelFont(43);
  h_Reco_energy_nu_cc1pi_M2->GetXaxis()->SetLabelSize(15);

  h_Reco_energy_nu_cc1pi_M2->SetTitle("") ;
  h_Reco_energy_nu_cc1pi_M2-> Draw();

  h_MC_energy_nu_cc1pi_M2->GetYaxis()->SetTitle("Events [#]");
  h_MC_energy_nu_cc1pi_M2->GetXaxis()->SetTitle("E#nu[GeV]");
  h_MC_energy_nu_cc1pi_M2->SetFillColor(kCyan-6);
  h_MC_energy_nu_cc1pi_M2->SetFillStyle(3004);
  h_MC_energy_nu_cc1pi_M2-> SetStats(0);
  h_MC_energy_nu_cc1pi_M2->SetLineWidth(1);
  h_MC_energy_nu_cc1pi_M2->SetTitle("") ;
  leg -> AddEntry( h_MC_energy_nu_cc1pi_M2 , " True" );
  h_MC_energy_nu_cc1pi_M2 -> Draw("same");
  leg -> AddEntry( h_Reco_energy_nu_cc1pi_M2 , " Selected " ); 
  
  h_TReco_energy_nu_cc1pi_M2->SetFillColor( 8 );
  h_TReco_energy_nu_cc1pi_M2->SetFillStyle(3144);
  h_TReco_energy_nu_cc1pi_M2-> SetStats(0);
  leg -> AddEntry( h_TReco_energy_nu_cc1pi_M2 , " Signal " );

  h_TReco_energy_nu_cc1pi_M2->SetTitle("") ;
  h_TReco_energy_nu_cc1pi_M2-> Draw("same");
  leg -> Draw();
  c->Print( (filepath+"h_energy_nu_cc1pi_M2.pdf").c_str() );
  c->SaveAs( (filepath+"h_energy_nu_cc1pi_M2.root").c_str() );
  
  c->Clear();
  leg->Clear();

  return 0;
} // MainTest()

