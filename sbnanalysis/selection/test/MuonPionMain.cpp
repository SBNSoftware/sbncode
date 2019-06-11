#include "../include/CC0piAnalysisHelper.h"
#include "../include/GeneralAnalysisHelper.h"
#include "../include/EventSelectionTool.h"
#include "../include/Event.h"
#include <iostream>
#include <sstream>
#include <numeric>
#include <time.h>
#include <algorithm>
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
  std::cout << " Start local time and date:  " << asctime(timeinfo)          << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;
 
  // Output file location
  std::string plots = "../Output_Selection_Tool/plots/muon_pion/";

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
   * Proton loss and mis-identification energy and length study
   *
   * Loop over events
   *   Loop over list of reconstructed particles
   *      If true == 2212 && reco == 2212: Matched
   *        Get kinetic energy of true
   *        Get length of true
   *      If reco != 2212 && true == 2212: Mis-identified true proton
   *        Get kinetic energy of true
   *        Get length of true
   *   Loop over list of MC particles
   *      If true == 2212
   *        Loop over reconstructed particles and see if any have corresponding MC ID
   *          If not: missed proton
   *            Get kinetic energy of true
   *            Get length of true
   *  
   */

  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();
  TopologyMap nc_signal_map    = GeneralAnalysisHelper::GetNCTopologyMap();
  TopologyMap cc0pi_signal_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  TopologyMap cc1pi_signal_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
  TopologyMap ccpi0_signal_map = GeneralAnalysisHelper::GetCCPi0TopologyMap();
  
  std::vector<float> signal_energy_mu, signal_length_mu, background_energy_mu, background_length_mu, missed_energy_mu, missed_length_mu;
  std::vector<float> signal_energy_pi, signal_length_pi, background_energy_pi, background_length_pi, missed_energy_pi, missed_length_pi;
  std::vector<int>   signal_hits_mu, background_hits_mu, missed_hits_mu;
  std::vector<int>   signal_hits_pi, background_hits_pi, missed_hits_pi;

  for(const Event &e : events){

    ParticleList reco_particles = e.GetRecoParticleList(); 
    ParticleList true_particles = e.GetMCParticleList();

    if(!e.IsSBNDTrueFiducial() || GeneralAnalysisHelper::NumberEscapingTracks(e) != 1) continue;
    
    for(Particle &p_reco : reco_particles){
      // If the particle is a reconstructed track (not pi0)
      if(p_reco.GetFromRecoTrack() && GeneralAnalysisHelper::ParticleHasAMatch(e, p_reco) >= 0){
        Particle p_match = GeneralAnalysisHelper::GetBestMCParticle(e,p_reco);
        if(p_match.GetPdgCode() == 13){
          if(p_reco.GetPdgCode() == 13){
            // True and reconstructed muon: signal
            signal_energy_mu.push_back(p_match.GetKineticEnergy());
            signal_length_mu.push_back(p_match.GetLength());
            signal_hits_mu.push_back(p_match.GetNumberOfHits());
          }
          else if(p_reco.GetPdgCode() != 13){
            // True and mis-indentified muon
            background_energy_mu.push_back(p_match.GetKineticEnergy());
            background_length_mu.push_back(p_match.GetLength());
            background_hits_mu.push_back(p_match.GetNumberOfHits());
          }
        }
        if(p_match.GetPdgCode() == 211){
          if(p_reco.GetPdgCode() == 211){
            // True and reconstructed pion: signal
            signal_energy_pi.push_back(p_match.GetKineticEnergy());
            signal_length_pi.push_back(p_match.GetLength());
            signal_hits_pi.push_back(p_match.GetNumberOfHits());
          }
          else if(p_reco.GetPdgCode() != 211){
            // True and mis-indentified pion
            background_energy_pi.push_back(p_match.GetKineticEnergy());
            background_length_pi.push_back(p_match.GetLength());
            background_hits_pi.push_back(p_match.GetNumberOfHits());
          }
        }
      }
    }
    for(Particle &p_true : true_particles){
      if(p_true.GetPdgCode() == 13 && !GeneralAnalysisHelper::HasBeenReconstructed(e, p_true)){
        // True muon, not reconstructed
        missed_energy_mu.push_back(p_true.GetKineticEnergy());
        missed_length_mu.push_back(p_true.GetLength());
        missed_hits_mu.push_back(p_true.GetNumberOfHits());
      }
      if(p_true.GetPdgCode() == 211 && !GeneralAnalysisHelper::HasBeenReconstructed(e, p_true)){
        // True pion, not reconstructed
        missed_energy_pi.push_back(p_true.GetKineticEnergy());
        missed_length_pi.push_back(p_true.GetLength());
        missed_hits_pi.push_back(p_true.GetNumberOfHits());
      }
    }
  }

  /*
   *
   * Fill
   *
   */
  TH1D *h_signal_energy_mu = new TH1D("h_signal_energy_mu","Correctly reconstructed muon kinetic energies",40,0,1.5);
  TH1D *h_signal_length_mu = new TH1D("h_signal_length_mu","Correctly reconstructed muon lengths",40,0,500);
  TH1D *h_signal_hits_mu   = new TH1D("h_signal_hits_mu",  "Correctly reconstructed muon hits",40,0,500);

  TH1D *h_missed_energy_mu = new TH1D("h_missed_energy_mu","Missed muon kinetic energies",40,0,1.5);
  TH1D *h_missed_length_mu = new TH1D("h_missed_length_mu","Missed muon lengths",40,0,500);
  TH1D *h_missed_hits_mu   = new TH1D("h_missed_hits_mu",  "Missed muon hits",40,0,500);

  TH1D *h_background_energy_mu = new TH1D("h_background_energy_mu","Mis-identified muon kinetic energies",40,0,1.5);
  TH1D *h_background_length_mu = new TH1D("h_background_length_mu","Mis-identified muon lengths",40,0,500);
  TH1D *h_background_hits_mu   = new TH1D("h_background_hits_mu",  "Mis-identified muon hits",40,0,500);

  TH1D *h_signal_energy_pi = new TH1D("h_signal_energy_pi","Correctly reconstructed pion kinetic energies",30,0,0.8);
  TH1D *h_signal_length_pi = new TH1D("h_signal_length_pi","Correctly reconstructed pion lengths",40,0,150);
  TH1D *h_signal_hits_pi   = new TH1D("h_signal_hits_pi",  "Correctly reconstructed pion hits",30,0,200);

  TH1D *h_missed_energy_pi = new TH1D("h_missed_energy_pi","Missed pion kinetic energies",30,0,0.8);
  TH1D *h_missed_length_pi = new TH1D("h_missed_length_pi","Missed pion lengths",40,0,150);
  TH1D *h_missed_hits_pi   = new TH1D("h_missed_hits_pi",  "Missed pion hits",30,0,200);

  TH1D *h_background_energy_pi = new TH1D("h_background_energy_pi","Mis-identified pion kinetic energies",30,0,0.8);
  TH1D *h_background_length_pi = new TH1D("h_background_length_pi","Mis-identified pion lengths",40,0,150);
  TH1D *h_background_hits_pi   = new TH1D("h_background_hits_pi",  "Mis-identified pion hits",30,0,200);

  for(unsigned int i = 0; i < signal_energy_mu.size(); ++i){
    h_signal_energy_mu->Fill(signal_energy_mu[i]);
    h_signal_length_mu->Fill(signal_length_mu[i]);
    h_signal_hits_mu->Fill(signal_hits_mu[i]);
  }
  for(unsigned int i = 0; i < background_energy_mu.size(); ++i){
    h_background_energy_mu->Fill(background_energy_mu[i]);
    h_background_length_mu->Fill(background_length_mu[i]);
    h_background_hits_mu->Fill(background_hits_mu[i]);
  }
  for(unsigned int i = 0; i < missed_energy_mu.size(); ++i){
    h_missed_energy_mu->Fill(missed_energy_mu[i]);
    h_missed_length_mu->Fill(missed_length_mu[i]);
    h_missed_hits_mu->Fill(missed_hits_mu[i]);
  }
  for(unsigned int i = 0; i < signal_energy_pi.size(); ++i){
    h_signal_energy_pi->Fill(signal_energy_pi[i]);
    h_signal_length_pi->Fill(signal_length_pi[i]);
    h_signal_hits_pi->Fill(signal_hits_pi[i]);
  }
  for(unsigned int i = 0; i < background_energy_pi.size(); ++i){
    h_background_energy_pi->Fill(background_energy_pi[i]);
    h_background_length_pi->Fill(background_length_pi[i]);
    h_background_hits_pi->Fill(background_hits_pi[i]);
  }
  for(unsigned int i = 0; i < missed_energy_pi.size(); ++i){
    h_missed_energy_pi->Fill(missed_energy_pi[i]);
    h_missed_length_pi->Fill(missed_length_pi[i]);
    h_missed_hits_pi->Fill(missed_hits_pi[i]);
  }

  /*
   *
   * Draw
   *
   */
  gStyle->SetPalette(55);
  gStyle->SetNumberContours(250);

  TCanvas *c = new TCanvas();
  TLegend *l = new TLegend( 0.48, 0.58, 0.88, 0.88 );
  l->SetBorderSize(0);

  l->AddEntry( h_signal_energy_mu,     " True, reconstructed muon", "l" );
  l->AddEntry( h_background_energy_mu, " True, misidentified muon", "l" );
  l->AddEntry( h_missed_energy_mu,     " True, missed muon", "l" );
  
  h_signal_energy_mu->SetLineColor(2);
  h_signal_energy_mu->SetStats(kFALSE);
  h_signal_energy_mu->Scale(1/double(signal_energy_mu.size()));
  h_background_energy_mu->SetLineColor(4);
  h_background_energy_mu->SetStats(kFALSE);
  h_background_energy_mu->GetXaxis()->SetTitle("Muon Kinetic Energy [GeV]");
  h_background_energy_mu->Scale(1/double(background_energy_mu.size()));
  h_missed_energy_mu->SetLineColor(8);
  h_missed_energy_mu->SetStats(kFALSE);
  h_missed_energy_mu->Scale(1/double(missed_energy_mu.size()));

  float max_y_e_mu = std::max(h_background_energy_mu->GetMaximum(), std::max(h_signal_energy_mu->GetMaximum(), h_missed_energy_mu->GetMaximum()));
  
  h_background_energy_mu->GetYaxis()->SetRangeUser(0,max_y_e_mu*1.1);
  h_background_energy_mu->Draw("hist");
  h_missed_energy_mu->Draw("hist same");
  h_signal_energy_mu->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"muon_kinetic_energy.root").c_str());
  c->SaveAs((plots+"muon_kinetic_energy.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_length_mu,     " True, reconstructed muon", "l" );
  l->AddEntry( h_background_length_mu, " True, misidentified muon", "l" );
  l->AddEntry( h_missed_length_mu,     " True, missed muon", "l" );
    
  h_signal_length_mu->SetLineColor(2);
  h_signal_length_mu->SetStats(kFALSE);
  h_signal_length_mu->Scale(1/double(signal_length_mu.size()));
  h_background_length_mu->SetLineColor(4);
  h_background_length_mu->SetStats(kFALSE);
  h_background_length_mu->GetXaxis()->SetTitle("Muon Length [cm]");
  h_background_length_mu->Scale(1/double(background_length_mu.size()));
  h_missed_length_mu->SetLineColor(8);
  h_missed_length_mu->SetStats(kFALSE);
  h_missed_length_mu->Scale(1/double(missed_length_mu.size()));

  float max_y_l_mu = std::max(h_background_length_mu->GetMaximum(), std::max(h_signal_length_mu->GetMaximum(),h_missed_length_mu->GetMaximum()));
  h_background_length_mu->GetYaxis()->SetRangeUser(0,max_y_l_mu*1.1);
  h_background_length_mu->Draw("hist");
  h_missed_length_mu->Draw("hist same");
  h_signal_length_mu->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"muon_length.root").c_str());
  c->SaveAs((plots+"muon_length.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_hits_mu,     " True, reconstructed muon", "l" );
  l->AddEntry( h_background_hits_mu, " True, misidentified muon", "l" );
  l->AddEntry( h_missed_hits_mu,     " True, missed muon", "l" );
  
  h_signal_hits_mu->SetLineColor(2);
  h_signal_hits_mu->SetStats(kFALSE);
  h_signal_hits_mu->Scale(1/double(signal_hits_mu.size()));
  h_background_hits_mu->SetLineColor(4);
  h_background_hits_mu->SetStats(kFALSE);
  h_background_hits_mu->GetXaxis()->SetTitle("Muon Hits");
  h_background_hits_mu->Scale(1/double(background_hits_mu.size()));
  h_missed_hits_mu->SetLineColor(8);
  h_missed_hits_mu->SetStats(kFALSE);
  h_missed_hits_mu->Scale(1/double(missed_hits_mu.size()));

  float max_y_h_mu = std::max(h_background_hits_mu->GetMaximum(), std::max(h_signal_hits_mu->GetMaximum(),h_missed_hits_mu->GetMaximum()));
  
  h_background_hits_mu->GetYaxis()->SetRangeUser(0,max_y_h_mu*1.1);
  h_background_hits_mu->Draw("hist");
  h_missed_hits_mu->Draw("hist same");
  h_signal_hits_mu->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"muon_hits.root").c_str());
  c->SaveAs((plots+"muon_hits.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_energy_pi,     " True, reconstructed pion", "l" );
  l->AddEntry( h_background_energy_pi, " True, misidentified pion", "l" );
  l->AddEntry( h_missed_energy_pi,     " True, missed pion", "l" );
  
  h_signal_energy_pi->SetLineColor(2);
  h_signal_energy_pi->SetStats(kFALSE);
  h_signal_energy_pi->Scale(1/double(signal_energy_pi.size()));
  h_background_energy_pi->SetLineColor(4);
  h_background_energy_pi->SetStats(kFALSE);
  h_background_energy_pi->GetXaxis()->SetTitle("Pion Kinetic Energy [GeV]");
  h_background_energy_pi->Scale(1/double(background_energy_pi.size()));
  h_missed_energy_pi->SetLineColor(8);
  h_missed_energy_pi->SetStats(kFALSE);
  h_missed_energy_pi->Scale(1/double(missed_energy_pi.size()));

  float max_y_e_pi = std::max(h_background_energy_pi->GetMaximum(), std::max(h_signal_energy_pi->GetMaximum(), h_missed_energy_pi->GetMaximum()));
  
  h_background_energy_pi->GetYaxis()->SetRangeUser(0,max_y_e_pi*1.1);
  h_background_energy_pi->Draw("hist");
  h_missed_energy_pi->Draw("hist same");
  h_signal_energy_pi->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"pion_kinetic_energy.root").c_str());
  c->SaveAs((plots+"pion_kinetic_energy.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_length_pi,     " True, reconstructed pion", "l" );
  l->AddEntry( h_background_length_pi, " True, misidentified pion", "l" );
  l->AddEntry( h_missed_length_pi,     " True, missed pion", "l" );
    
  h_signal_length_pi->SetLineColor(2);
  h_signal_length_pi->SetStats(kFALSE);
  h_signal_length_pi->Scale(1/double(signal_length_pi.size()));
  h_background_length_pi->SetLineColor(4);
  h_background_length_pi->SetStats(kFALSE);
  h_background_length_pi->GetXaxis()->SetTitle("Pion Length [cm]");
  h_background_length_pi->Scale(1/double(background_length_pi.size()));
  h_missed_length_pi->SetLineColor(8);
  h_missed_length_pi->SetStats(kFALSE);
  h_missed_length_pi->Scale(1/double(missed_length_pi.size()));

  float max_y_l_pi = std::max(h_background_length_pi->GetMaximum(), std::max(h_signal_length_pi->GetMaximum(),h_missed_length_pi->GetMaximum()));
  h_background_length_pi->GetYaxis()->SetRangeUser(0,max_y_l_pi*1.1);
  h_background_length_pi->Draw("hist");
  h_missed_length_pi->Draw("hist same");
  h_signal_length_pi->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"pion_length.root").c_str());
  c->SaveAs((plots+"pion_length.png").c_str());
  c->Clear();

  l->Clear();
  l->AddEntry( h_signal_hits_pi,     " True, reconstructed pion", "l" );
  l->AddEntry( h_background_hits_pi, " True, misidentified pion", "l" );
  l->AddEntry( h_missed_hits_pi,     " True, missed pion", "l" );
  
  h_signal_hits_pi->SetLineColor(2);
  h_signal_hits_pi->SetStats(kFALSE);
  h_signal_hits_pi->Scale(1/double(signal_hits_pi.size()));
  h_background_hits_pi->SetLineColor(4);
  h_background_hits_pi->SetStats(kFALSE);
  h_background_hits_pi->GetXaxis()->SetTitle("Pion Hits");
  h_background_hits_pi->Scale(1/double(background_hits_pi.size()));
  h_missed_hits_pi->SetLineColor(8);
  h_missed_hits_pi->SetStats(kFALSE);
  h_missed_hits_pi->Scale(1/double(missed_hits_pi.size()));

  float max_y_h_pi = std::max(h_background_hits_pi->GetMaximum(), std::max(h_signal_hits_pi->GetMaximum(),h_missed_hits_pi->GetMaximum()));
  
  h_background_hits_pi->GetYaxis()->SetRangeUser(0,max_y_h_pi*1.1);
  h_background_hits_pi->Draw("hist");
  h_missed_hits_pi->Draw("hist same");
  h_signal_hits_pi->Draw("hist same");
  l->Draw();

  c->SaveAs((plots+"pion_hits.root").c_str());
  c->SaveAs((plots+"pion_hits.png").c_str());
  c->Clear();

  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end)       << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
