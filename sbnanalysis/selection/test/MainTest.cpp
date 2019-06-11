#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include "TH1.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
using namespace selection;

int MainTest(){

  //std::string filename = "/hepstore/rjones/Samples/LArSoft_reconstruction_temp/single_file_tree_test.root";
  //std::string filename = "/hepstore/rjones/Samples/LArSoft_reconstruction_temp/30_tree_test.root" ;
   std::string filename = "/hepstore/rjones/Samples/FNAL/analysis_trees/analysis_trees_2.root" ;
  //  std::string filename = "/hepstore/rjones/Samples/FNAL/analysis_trees/analysis_trees_10.root" ;
  //std::string filename = "/hepstore/rjones/Samples/FNAL/analysis_trees/analysis_trees_100.root" ;
  //std::string filename = "/hepstore/rjones/Samples/FNAL/analysis_trees/analysis_trees_200.root" ;
  //std::string filename = "/hepstore/rjones/Samples/FNAL/analysis_trees/analysis_trees_398.root" ;

  ofstream rfile, efile, lfile, afile ;
  rfile.open( "results.txt" ) ;
  efile.open( "event_information.txt" ) ;
  lfile.open( "length_information_per_event.txt" ) ;
  afile.open( "angular_information_per_event.txt" );
  EventSelectionTool::EventList events;
  EventSelectionTool::LoadEventList(filename, events);
  //-----------------------------------
  // Parameters for Efficiency calculation ( MC, TReco and Reco ) for a given topology : 0-> NC , 1 -> CCinclusive, 2-> CC0pi, 3-> CC1pi+/-, 4-> CC1pi0
  std::vector< double > CountMC   ( 5 ) ;
  std::vector< double > CountTReco ( 5 ) ;
  std::vector< double > CountReco ( 5 ) ;

  //-----------------------------------
  // Count per ID
  // ID = mu, p, pi, pi0 [ 0, 1, 2, 3 ]
  std::vector< vector< double > > Count_MC_ID (4,vector<double>(4));   // Number of events with 0, 1, 2, >2 mu
  std::vector< vector< double > > Count_TReco_ID (4,vector<double>(4)); // Number of events with 0, 1, 2, >2 mu
  std::vector< vector< double > > Count_Reco_ID (4,vector<double>(4)); // Number of events with 0, 1, 2, >2 mu
  // Count ExChanged particles
  std::vector< vector< double > > Count_ExChange_MC (4,vector<double>(4)); // Count events in change
  std::vector< vector< double > > Count_ExChange_TReco (4,vector<double>(4)); // Change of pdg(i) -> pdg(j) (MC == Reco)
  std::vector< vector< double > > Count_ExChange_Reco (4,vector<double>(4)); // Change of pdg(i) -> pdg(j) (MC != Reco)
  //Topology matrix
  std::vector< vector< double> > Count_MC_Topology ( 5, vector< double > (5));
  std::vector< vector< double> > Count_TReco_Topology ( 5, vector< double > (5));
  std::vector< vector< double> > Count_Reco_Topology ( 5, vector< double > (5));
  //-----------------------------------
  // Length calculation variables
  // ---------------------------------
  std::vector< vector< double > > Count_L (3,vector<double>(3)); // [i][j]: i=MC,TReco,Reco, j=mu,pi,p
  TH1 * h_mu_MC_L       = new TH1D ( "h_mu_MC_L"      , " MC: mu track length "     , 30, 0., 1000. );
  TH1 * h_mu_TReco_L    = new TH1D ( "h_mu_TReco_L"   , "TReco: mu track length "   , 30, 0., 300.  );
  TH1 * h_mu_Reco_L     = new TH1D ( "h_mu_Reco_L"    , "Reco: mu track length "    , 30, 0., 300.  );

  TH1 * h_pi_MC_L       = new TH1D ( "h_pi_MC_L"      , " MC: pi track length "     , 20, 0., 200.  );
  TH1 * h_pi_TReco_L    = new TH1D ( "h_pi_TReco_L"   , "TReco: pi track length "   , 20, 0., 200.  );
  TH1 * h_pi_Reco_L     = new TH1D ( "h_pi_Reco_L"    , "Reco: pi track length "    , 20, 0., 200.  );

  TH1 * h_p_MC_L        = new TH1D ( "h_pi_MC_L"      , " MC: p track length "      , 20, 0., 100.  );
  TH1 * h_p_TReco_L     = new TH1D ( "h_pi_TReco_L"   , "TReco: p track length "    , 20, 0., 100.  );
  TH1 * h_p_Reco_L      = new TH1D ( "h_pi_Reco_L"    , "Reco: p track length "     , 20, 0., 40.   );

  TH1 * h_mu_pi_MC_L    = new TH1D ( "h_mu_pi_MC_L"   , " MC: mu/pi track length "  , 50, 0., 10.   );
  TH1 * h_mu_pi_TReco_L = new TH1D ( "h_mu_pi_TReco_L", "TReco: mu/pi track length ", 50, 0., 10.   );
  TH1 * h_mu_pi_Reco_L  = new TH1D ( "h_mu_pi_Reco_L" , "Reco: mu/pi track length " , 50, 0., 10.   );

  //-----------------------------------
  // Angle variables
  // ---------------------------------
  TH1 * h_mu_MC_T       = new TH1D ( "h_mu_MC_T"      , " MC: mu track angle "     , 30, -1., 1. );
  TH1 * h_mu_TReco_T    = new TH1D ( "h_mu_TReco_T"   , "TReco: mu track angle "   , 30, -1., 1. );
  TH1 * h_mu_Reco_T     = new TH1D ( "h_mu_Reco_T"    , "Reco: mu track angle "    , 30, -1., 1. );

  TH1 * h_pi_MC_T       = new TH1D ( "h_pi_MC_T"      , " MC: pi track angle "     , 30, -1., 1. );
  TH1 * h_pi_TReco_T    = new TH1D ( "h_pi_TReco_T"   , "TReco: pi track angle "   , 30, -1., 1. );
  TH1 * h_pi_Reco_T     = new TH1D ( "h_pi_Reco_T"    , "Reco: pi track angle "    , 30, -1., 1. );

  TH1 * h_p_MC_T        = new TH1D ( "h_pi_MC_T"      , " MC: p track angle "      , 30, -1., 1. );
  TH1 * h_p_TReco_T     = new TH1D ( "h_pi_TReco_T"   , "TReco: p track angle "    , 30, -1., 1. );
  TH1 * h_p_Reco_T      = new TH1D ( "h_pi_Reco_T"    , "Reco: p track angle "     , 30, -1., 1. );

  TH1 * h_mu_pi_MC_T    = new TH1D ( "h_mu_pi_MC_T"   , " MC: mu/pi track angle "  , 30, -1., 1. );
  TH1 * h_mu_pi_TReco_T = new TH1D ( "h_mu_pi_TReco_T", "TReco: mu/pi track angle ", 30, -1., 1. );
  TH1 * h_mu_pi_Reco_T  = new TH1D ( "h_mu_pi_Reco_T" , "Reco: mu/pi track angle " , 30, -1., 1. );


  for(unsigned int i = 0; i < events.size(); ++i){

    // Do analysis
    Event &e(events[i]);
    e.SetTopologies();

    efile << "-----------------------------------------------------------" << std::endl;
    efile << " RECO : " << std::endl;
    efile << "   muons   : " << e.CountRecoParticlesWithPdg(13) << "\n";
    efile << "   pi+/-   : " << e.CountRecoParticlesWithPdg(211) + e.CountRecoParticlesWithPdg(-211) << "\n";
    efile << "   pi0     : " << e.CountRecoParticlesWithPdg(111) << "\n";
    efile << "   protons : " << e.CountRecoParticlesWithPdg(2212) << "\n";
    efile << "   Topology: " << "\n";
    efile << "   NC      : " << e.CheckRecoTopology(e.signal_map_NC) << "\n";
    efile << "   ccincl. : " << e.CheckRecoTopology(e.signal_map_cc_inclusive) << "\n";
    efile << "   cc0pi   : " << e.CheckRecoTopology(e.signal_map_cc_0pi) << "\n";
    efile << "   cc1pi   : " << e.CheckRecoTopology(e.signal_map_cc_1pi) << "\n";
    efile << "   ccpi0   : " << e.CheckRecoTopology(e.signal_map_cc_pi0) << "\n";

    efile << "   MC      : " << "\n";
    efile << "   muons   : " << e.CountMCParticlesWithPdg(13) << "\n";
    efile << "   pi+/-   : " << e.CountMCParticlesWithPdg(211) + e.CountMCParticlesWithPdg(-211) << "\n";
    efile << "   pi0     : " << e.CountMCParticlesWithPdg(111) << "\n";
    efile << "   protons : " << e.CountMCParticlesWithPdg(2212) << "\n";
    efile << "   Topology: " << "\n";
    efile << "   NC      : " << e.CheckMCTopology(e.signal_map_NC) << "\n";
    efile << "   ccincl. : " << e.CheckMCTopology(e.signal_map_cc_inclusive) << "\n";
    efile << "   cc0pi   : " << e.CheckMCTopology(e.signal_map_cc_0pi) << "\n";
    efile << "   cc1pi   : " << e.CheckMCTopology(e.signal_map_cc_1pi) << "\n";
    efile << "   ccpi0   : " << e.CheckMCTopology(e.signal_map_cc_pi0) << "\n";


    // ----------------------------Eficiency calculation values---------------------------
    e.Count_per_Topology( e.signal_map_NC           , CountMC[0], CountTReco[0], CountReco[0] );
    e.Count_per_Topology( e.signal_map_cc_inclusive , CountMC[1], CountTReco[1], CountReco[1] );
    e.Count_per_Topology( e.signal_map_cc_0pi       , CountMC[2], CountTReco[2], CountReco[2] );
    e.Count_per_Topology( e.signal_map_cc_1pi       , CountMC[3], CountTReco[3], CountReco[3] );
    e.Count_per_Topology( e.signal_map_cc_pi0       , CountMC[4], CountTReco[4], CountReco[4] );


    //--------------------------- COUNT PARTICLES RECONSTRUCTED --------------------------
    /*
     e.ParticleReconstruction( Count_MC_ID, Count_TReco_ID, Count_Reco_ID, 13 );
     e.ParticleReconstruction( Count_MC_ID, Count_TReco_ID, Count_Reco_ID, 2212 );
     e.ParticleReconstruction( Count_MC_ID, Count_TReco_ID, Count_Reco_ID, 211 );
     e.ParticleReconstruction( Count_MC_ID, Count_TReco_ID, Count_Reco_ID, 111 );
     e.ParticleExChange( Count_ExChange_MC, Count_ExChange_TReco,  Count_ExChange_Reco ) ;*/
    //-------------------------- TOPOLOGY MATRIX ----------------------------------------
     e.TopologyMatrix(Count_MC_Topology, Count_TReco_Topology,Count_Reco_Topology ) ;

     //------------------------  LENGTH CC1PI TRACKS ------------------------------------
     e.CountLength_topology( e.signal_map_cc_1pi , Count_L );

     if( e.CheckMCTopology( e.signal_map_cc_1pi ) == 1 ) {
       lfile << "Event #" << i << "\n" ;
       lfile << "MC events :"  << "\n" ;
       lfile << "__________________________________________________________" << "\n";
       lfile << "Muon length:   " << e.GetMCLengthWithPdg( 13 )              << "\n" ;
       lfile << "Pion length:   " << e.GetMCLengthWithPdg( 211 )             << "\n" ;
       lfile << "Proton length: " << e.GetMCLengthWithPdg( 2212 )            << "\n" ;

       h_mu_MC_L -> Fill( e.GetMCLengthWithPdg( 13 ) );
       h_pi_MC_L -> Fill( e.GetMCLengthWithPdg( 211 ) );
       h_p_MC_L -> Fill( e.GetMCLengthWithPdg( 2212 ) );
       h_mu_pi_MC_L -> Fill( e.GetMCLengthWithPdg( 13 )/e.GetMCLengthWithPdg( 211 ) );

       afile << "__________________________________________________________" << "\n";
       afile << "Angle : Event #" << i << "\n" ;
       afile << "MC events :"  << "\n" ;
       afile << "__________________________________________________________" << "\n";
       afile << "Muon angle:   " << e.GetMCCosThetaWithPdg( 13 )              << "\n" ;
       afile << "Pion angle:   " << e.GetMCCosThetaWithPdg( 211 )             << "\n" ;
       afile << "Proton angle: " << e.GetMCCosThetaWithPdg( 2212 )            << "\n" ;
       h_mu_MC_T -> Fill( e.GetMCCosThetaWithPdg( 13 ) );
       h_pi_MC_T -> Fill( e.GetMCCosThetaWithPdg( 211 ) );
       h_p_MC_T  -> Fill( e.GetMCCosThetaWithPdg( 2212 ) );


       if( e.CheckRecoTopology( e.signal_map_cc_1pi ) == 1) {

	 lfile << "__________________________________________________________" << "\n";
	 lfile << "\n";
	 lfile << "True Reco events :"  << "\n" ;
	 lfile << "__________________________________________________________" << "\n";
	 lfile << "Muon length:   " << e.GetRecoLengthWithPdg( 13 )            << "\n" ;
	 lfile << "Pion length:   " << e.GetRecoLengthWithPdg( 211 )           << "\n" ;
	 lfile << "Proton length: " << e.GetRecoLengthWithPdg( 2212 )          << "\n" ;

	 h_mu_TReco_L -> Fill( e.GetRecoLengthWithPdg( 13 ) );
	 h_pi_TReco_L -> Fill( e.GetRecoLengthWithPdg( 211 ) );
	 h_p_TReco_L -> Fill( e.GetRecoLengthWithPdg( 2212 ) );
	 h_mu_pi_TReco_L -> Fill( e.GetRecoLengthWithPdg( 13 )/e.GetRecoLengthWithPdg( 211 ) );

	 afile << "__________________________________________________________" << "\n";
	 afile << "Angle : Event #" << i << "\n" ;
	 afile << "TReco events :"  << "\n" ;
	 afile << "__________________________________________________________" << "\n";
	 afile << "Muon angle:   " << e.GetRecoCosThetaWithPdg( 13 )              << "\n" ;
	 afile << "Pion angle:   " << e.GetRecoCosThetaWithPdg( 211 )             << "\n" ;
	 afile << "Proton angle: " << e.GetRecoCosThetaWithPdg( 2212 )            << "\n" ;
	 h_mu_TReco_T -> Fill( e.GetRecoCosThetaWithPdg( 13 ) );
	 h_pi_TReco_T -> Fill( e.GetRecoCosThetaWithPdg( 211 ) );
	 h_p_TReco_T  -> Fill( e.GetRecoCosThetaWithPdg( 2212 ) );


       }

     }
     if( e.CheckRecoTopology(e.signal_map_cc_1pi) == 1) {

       lfile << "__________________________________________________________" << "\n";
       lfile << "\n";
       lfile << "Reco events :"  << "\n" ;
       lfile << "__________________________________________________________" << "\n";
       lfile << "Muon length:   " << e.GetRecoLengthWithPdg( 13 )            << "\n" ;
       lfile << "Pion length:   " << e.GetRecoLengthWithPdg( 211 )           << "\n" ;
       lfile << "Proton length: " << e.GetRecoLengthWithPdg( 2212 )          << "\n" ;

       h_mu_Reco_L -> Fill( e.GetRecoLengthWithPdg( 13 ) );
       h_pi_Reco_L -> Fill( e.GetRecoLengthWithPdg( 211 ) );
       h_p_Reco_L -> Fill( e.GetRecoLengthWithPdg( 2212 ) );
       h_mu_pi_Reco_L -> Fill( e.GetRecoLengthWithPdg( 13 )/e.GetMCLengthWithPdg( 211 ));


       afile << "__________________________________________________________" << "\n";
       afile << "Angle : Event #" << i << "\n" ;
       afile << "Reco events :"  << "\n" ;
       afile << "__________________________________________________________" << "\n";
       afile << "Muon angle:   " << e.GetRecoCosThetaWithPdg( 13 )              << "\n" ;
       afile << "Pion angle:   " << e.GetRecoCosThetaWithPdg( 211 )             << "\n" ;
       afile << "Proton angle: " << e.GetRecoCosThetaWithPdg( 2212 )            << "\n" ;
       h_mu_Reco_T -> Fill( e.GetRecoCosThetaWithPdg( 13 ) );
       h_pi_Reco_T -> Fill( e.GetRecoCosThetaWithPdg( 211 ) );
       h_p_Reco_T  -> Fill( e.GetRecoCosThetaWithPdg( 2212 ) );
     }


  } //endfor


  for( int i = 0; i<5; ++i ){
    rfile << "__________________________________________________________" << "\n";
    rfile << "\n" ;
    rfile << "                 TOPOLOGY NUMBER " << i << "\n";
    rfile << "__________________________________________________________" << "\n";
    rfile << "CountMC = " << CountMC[i] << "\n";
    rfile << "CountReco = " << CountReco[i] << "\n";
    rfile << "SameCount = " << CountTReco[i] << "\n";
    rfile << "Background = " << CountReco[i] - CountTReco[i] << "\n";
    rfile << "Correct Reconstructed Events[%]=" << ( CountTReco[i] / CountMC[i] ) * 100 << "\n";
    rfile << "MisIdentified[%]=" << (( CountReco[i] - CountTReco[i] ) / CountMC[i] ) * 100 << "\n";
    rfile << "__________________________________________________________" << "\n";
 }


  rfile << "\n";
  rfile << "____________________________________________________________"<< "\n";
  rfile << "\n";
  rfile << "        TOPOLOGY MATRIX MC ( #_MC / Total_MC )  : "<< "\n";
  rfile << "____________________________________________________________"<< "\n";
  for( unsigned int i = 0 ; i < 5; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 5 ; ++k ) {
      if( i == k && Count_MC_Topology[i][k]!=0 ){
        rfile << ( Count_MC_Topology[i][k] / events.size() ) * 100 << "   ,   ";
      } else rfile << "   --   ";
    }
    rfile <<")"<< "\n";
  }

rfile << "\n";
  rfile << "____________________________________________________________"<< "\n";
  rfile << "\n";
  rfile << "   TOPOLOGY MATRIX - TRUE RECO  (#_TReco / Total_MC) : "<< "\n";
  rfile << "____________________________________________________________"<< "\n";
  for( unsigned int i = 0 ; i < 5; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 5 ; ++k ) {
      if( Count_TReco_Topology[i][k]!=0 ){
        rfile << ( Count_TReco_Topology[i][k] / events.size() ) * 100 << "   ,   ";
      } else rfile << "   --   ";
    }
    rfile <<")"<< "\n";
  }

  rfile << "\n";
  rfile << "____________________________________________________________"<< "\n";
  rfile << "\n";
  rfile << "TOPOLOGY MATRIX - BACKGROUND ((#_Reco-#_TReco) / Total_MC ) : "<< "\n";
  rfile << "____________________________________________________________"<< "\n";
  for( unsigned int i = 0 ; i < 5; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 5 ; ++k ) {
        rfile << (( Count_Reco_Topology[i][k] - Count_TReco_Topology[i][k] ) / events.size() ) * 100 << "   ,   ";
    }
    rfile <<")"<< "\n";
  }


  rfile << "\n";
  rfile << "____________________________________________________________"<< "\n";
  rfile << "\n";
  rfile << "  TOPOLOGY MATRIX - TRUE RECO  (#_TReco / #_Total_MC) : "<< "\n";
  rfile << "____________________________________________________________"<< "\n";
  for( unsigned int i = 0 ; i < 5; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 5 ; ++k ) {
      if( Count_TReco_Topology[i][k]!=0 ){
        rfile << ( Count_TReco_Topology[i][k] / Count_MC_Topology[i][k] ) * 100 << "   ,   ";
      } else rfile << "   --   ";
    }
    rfile <<")"<< "\n";
  }


  lfile << "\n";
  lfile << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << "\n";
  lfile << "                    LENGTH CHECK MC                         "         << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << " # L(mu) > L(pi) & L(p) =   " << ( Count_L[0][0] / CountMC[3] ) * 100 << "\n";
  lfile << " # L(p)  > L(pi)        =   " << ( Count_L[0][2] / CountMC[3] ) * 100 << "\n";
  lfile << " # L(pi) > L(p)         =   " << ( Count_L[0][1] / CountMC[3] ) * 100 << "\n";
  lfile << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << "\n";
  lfile << "                    LENGTH CHECK TReco                       "         << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << " # L(p)  > L(pi)        =   " << ( Count_L[1][2] / CountTReco[3] ) * 100 << "\n";
  lfile << " # L(pi) > L(p)         =   " << ( Count_L[1][1] / CountTReco[3] ) * 100 << "\n";
  lfile << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << "\n";
  lfile << "                    LENGTH CHECK Reco                       "         << "\n";
  lfile << "____________________________________________________________"         << "\n";
  lfile << " # L(mu) > L(pi) & L(p) =   " << ( Count_L[2][0] / CountReco[3] ) * 100 << "\n";
  lfile << " # L(p) > L(pi)         =   " << ( Count_L[2][2] / CountReco[3] ) * 100 << "\n";
  lfile << " # L(pi) > L(p)         =   " << ( Count_L[2][1] / CountReco[3] ) * 100 << "\n";


  // Plot histograms
  TCanvas * c = new TCanvas( "c_mu_MC_L", "c_mu_MC_L", 600, 600 );
  h_mu_MC_L->GetYaxis()->SetTitle("Events");
  h_mu_MC_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_MC_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_MC_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_MC_L.root");
  c->Clear();

  h_mu_TReco_L->GetYaxis()->SetTitle("Events");
  h_mu_TReco_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_TReco_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_TReco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_TReco_L.root");
  c->Clear();

  h_mu_Reco_L->GetYaxis()->SetTitle("Events");
  h_mu_Reco_L->GetXaxis()->SetTitle("Length [cm]");
  h_mu_Reco_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_Reco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_Reco_L.root");
  c->Clear();

  h_pi_MC_L->GetYaxis()->SetTitle("Events");
  h_pi_MC_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_MC_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_MC_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_MC_L.root");
  c->Clear();

  h_pi_TReco_L->GetYaxis()->SetTitle("Events");
  h_pi_TReco_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_TReco_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_TReco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_TReco_L.root");
  c->Clear();

  h_pi_Reco_L->GetYaxis()->SetTitle("Events");
  h_pi_Reco_L->GetXaxis()->SetTitle("Length [cm]");
  h_pi_Reco_L->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_Reco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_Reco_L.root");
  c->Clear();

  h_p_MC_L->GetYaxis()->SetTitle("Events");
  h_p_MC_L->GetXaxis()->SetTitle("Length [cm]");
  h_p_MC_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_MC_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_MC_L.root");
  c->Clear();

  h_p_TReco_L->GetYaxis()->SetTitle("Events");
  h_p_TReco_L->GetXaxis()->SetTitle("Length [cm]");
  h_p_TReco_L-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/h_p_TReco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/h_p_TReco_L.root");
  c->Clear();

  h_p_Reco_L->GetYaxis()->SetTitle("Events");
  h_p_Reco_L->GetXaxis()->SetTitle("Length [cm]");
  h_p_Reco_L->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_Reco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_Reco_L.root");
  c->Clear();

  h_mu_pi_MC_L->GetYaxis()->SetTitle("Events");
  h_mu_pi_MC_L->GetXaxis()->SetTitle("Ratio");
  h_mu_pi_MC_L->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_MC_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_MC_L.root");
  c->Clear();

  h_mu_pi_TReco_L->GetYaxis()->SetTitle("Events");
  h_mu_pi_TReco_L->GetXaxis()->SetTitle("Ratio");
  h_mu_pi_TReco_L->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_TReco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_TReco_L.root");
  c->Clear();

  h_mu_pi_Reco_L->GetYaxis()->SetTitle("Events");
  h_mu_pi_Reco_L->GetXaxis()->SetTitle("Ratio");
  h_mu_pi_Reco_L->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_Reco_L.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_pi_Reco_L.root");
  c->Clear();


  //ANGLE
  h_mu_MC_T->GetYaxis()->SetTitle("Events");
  h_mu_MC_T->GetXaxis()->SetTitle("Cos(theta)");
  h_mu_MC_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_MC_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_MC_T.root");
  c->Clear();

  h_mu_TReco_T->GetYaxis()->SetTitle("Events");
  h_mu_TReco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_mu_TReco_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_TReco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_TReco_T.root");
  c->Clear();

  h_mu_Reco_T->GetYaxis()->SetTitle("Events");
  h_mu_Reco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_mu_Reco_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_Reco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_mu_Reco_T.root");
  c->Clear();

  h_pi_MC_T->GetYaxis()->SetTitle("Events");
  h_pi_MC_T->GetXaxis()->SetTitle("Cos(theta)");
  h_pi_MC_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_MC_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_MC_T.root");
  c->Clear();

  h_pi_TReco_T->GetYaxis()->SetTitle("Events");
  h_pi_TReco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_pi_TReco_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_TReco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_TReco_T.root");
  c->Clear();

  h_pi_Reco_T->GetYaxis()->SetTitle("Events");
  h_pi_Reco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_pi_Reco_T->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_Reco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_pi_Reco_T.root");
  c->Clear();

  h_p_MC_T->GetYaxis()->SetTitle("Events");
  h_p_MC_T->GetXaxis()->SetTitle("Cos(theta)");
  h_p_MC_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_MC_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_MC_T.root");
  c->Clear();

  h_p_TReco_T->GetYaxis()->SetTitle("Events");
  h_p_TReco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_p_TReco_T-> Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_TReco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_TReco_T.root");
  c->Clear();

  h_p_Reco_T->GetYaxis()->SetTitle("Events");
  h_p_Reco_T->GetXaxis()->SetTitle("Cos(theta)");
  h_p_Reco_T->Draw();
  c->Print( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_Reco_T.jpg");
  c->SaveAs( "~/Desktop/reconstruction/Event_Selection_Tool/hist_output/test/h_p_Reco_T.root");
  c->Clear();




  return 0;

} // MainTest()




  /*  rfile <<"__________________________________________________________" << "\n";
  rfile << "\n";
  rfile <<"               PARTICLE IDENTIFICATION:   "<< "\n";
  rfile <<"__________________________________________________________" << "\n";
  rfile << "---------------------------------------"<< "\n";
  rfile << "-----------  Fraction W-R   -----------"<< "\n";
  rfile << "---------------------------------------"<< "\n";
  for( unsigned int i = 0 ; i < 3; ++i ){
    rfile << "(    ";
    for( unsigned int k = 0 ; k < 4 ; ++k ) {
      if(CountMC_ID[i][k] != 0) {
      	rfile << (CountSame_ID[i][k]/CountMC_ID[i][k])*100 << "    ,    ";
      }else  rfile << "  --  " ;
    }
    rfile << "    )" << "\n";
  }

  rfile << "---------------------------------------"<< "\n";
  rfile << "---------  Fraction Miss-R   ----------"<< "\n";
  rfile << "---------------------------------------"<< "\n";
  for( unsigned int i = 0 ; i < 3; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 4 ; ++k ) {
      if(CountReco_ID[i][k] != 0) {
	rfile << ((CountReco_ID[i][k]-CountSame_ID[i][k])/CountReco_ID[i][k])*100 << "   ,   ";
      } else   rfile << "  --  " ;
    }
    rfile << ")" << "\n";
  }
  rfile << "\n";
  rfile << "\n";

  rfile << "__________________________________________________________"<< "\n";
  rfile << "\n";
  rfile << "                   PARTICLE EXCHANGE: "<< "\n";
  rfile << "__________________________________________________________"<< "\n";
  rfile << "-----------------------SAME------------------------------"<< "\n";
  rfile << "------MC( 1p_i - 0p_j ) == Reco( 1p_i - 0p_j )-----------"<< "\n";
  rfile << "---------------------------------------------------------"<< "\n";
  for( unsigned int i = 0 ; i < 4; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 4 ; ++k ) {
      if( i != k  && Count_ExChange_MC[i][k]!=0  ){
	rfile << Count_ExChange_Same[i][k]/Count_ExChange_MC[i][k]*100 << "   ,   ";
      }else rfile << " -- "  ;
    }
    rfile << " )"<< "\n";
  }

  rfile << "-------------------------ExCHANGE------------------------" << "\n";
  rfile << "------MC( 1p_i - 0p_j ) == Reco( 0p_i - 1p_j )-----------" << "\n";
  rfile << "---------------------------------------------------------" << "\n";
  for( unsigned int i = 0 ; i < 4; ++i ){
    rfile << "(";
    for( unsigned int k = 0 ; k < 4 ; ++k ) {
      if( i != k && Count_ExChange_MC[i][k]!=0 ){
	rfile << Count_ExChange[i][k]/Count_ExChange_MC[i][k]*100 << "   ,   ";
      }else rfile << " -- " ;
    }
    rfile <<")"<< "\n";
  }
  */
