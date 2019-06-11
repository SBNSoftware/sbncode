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
#include <cmath>

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
  std::string plots_location = "../Output_Selection_Tool/plots/proton_resolution/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 500;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){

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

  TopologyMap cc_signal_map    = GeneralAnalysisHelper::GetCCIncTopologyMap();

  // Counters for true protons
  unsigned int detectable_protons               = 0;
  unsigned int reconstructed_detectable_protons = 0;
  unsigned int missed_detectable_protons        = 0;

  // Histograms 
  TH1D *h_detectable_protons     = new TH1D("h_detectable_protons",     "", 10, 0, 10);
  TH1D *h_non_detectable_protons = new TH1D("h_non_detectable_protons", "", 10, 0, 10);
  TH1D *h_ke_reconstructed       = new TH1D("h_ke_reconstructed", "", 80, 0, 0.5);
  TH1D *h_ke_total               = new TH1D("h_ke_total", "", 80, 0, 0.5);
  TH1D *h_length_reconstructed   = new TH1D("h_length_reconstructed", "", 80, 0, 30);
  TH1D *h_length_total           = new TH1D("h_length_total", "", 80, 0, 30);
  TH1D *h_hits_reconstructed     = new TH1D("h_hits_reconstructed", "", 80, 0, 150);
  TH1D *h_hits_total             = new TH1D("h_hits_total", "", 80, 0, 150);
  TH1D *h_dist_reconstructed     = new TH1D("h_dist_reconstructed", "", 80, 0, 100);
  TH1D *h_dist_missed            = new TH1D("h_dist_missed", "", 80, 0, 100);
  TH1D *h_dist_total             = new TH1D("h_dist_total", "", 80, 0, 100);
  TH1D *h_angle_reconstructed    = new TH1D("h_angle_reconstructed", "", 80, -1, 1);
  TH1D *h_angle_missed           = new TH1D("h_angle_missed", "", 80, -1, 1);
  TH1D *h_angle_total            = new TH1D("h_angle_total", "", 80, -1, 1);

  // First, ensure a maximum of 1 track escapes and the topology is numu CC
  for(const Event &e : events){
    if(!e.IsSBNDTrueFiducial() || !GeneralAnalysisHelper::MaxOneEscapingTrack(e) || !e.CheckMCTopology(cc_signal_map)) continue;

    // Event-by-event counters for proton multiplicity
    unsigned int detectable_multiplicity     = 0;
    unsigned int non_detectable_multiplicity = 0;

    // Muon variables
    TVector3 muon_traj(0.,0.,0.);
    float    muon_length    = 0.;
    
    // Get true particle list and loop over them
    ParticleList true_particles = e.GetMCParticleList();
    for(const Particle &p : true_particles){
      if(p.GetPdgCode() != 13) continue;
      
      muon_traj      = p.GetEnd() - p.GetVertex();
      muon_length    = (muon_traj).Mag();
    }
    for(const Particle &p : true_particles){
      // Get the true protons
      if(p.GetPdgCode() != 2212) continue;

      // Find out if they are detectable
      if(p.GetKineticEnergy() >= 0.02) {
    
        // The proton/proton-muon variables
        TVector3 proton_traj(0.,0.,0.);
        float proton_length = 0.;
        float theta_mu_p    = 0.;
        float distance_mu_p = 0.;

        detectable_protons++;
        detectable_multiplicity++;

        // Proton variables
        proton_traj      = p.GetEnd() - p.GetVertex();
        proton_length    = (proton_traj).Mag();

        // Muon-proton variables
        theta_mu_p    = acos((1./(proton_length*muon_length)) * (proton_traj.Dot(muon_traj)));
        distance_mu_p = proton_length * sin(theta_mu_p);

        // Efficiency denominators
        h_ke_total->Fill(p.GetKineticEnergy());
        h_length_total->Fill(p.GetLength());
        h_hits_total->Fill(p.GetNumberOfHits());
        h_dist_total->Fill(distance_mu_p);
        h_angle_total->Fill(cos(theta_mu_p));
        // Find out if they were reconstructed
        if(GeneralAnalysisHelper::HasBeenReconstructed(e,p)){
          reconstructed_detectable_protons++;
          
          // Efficiency plots
          h_ke_reconstructed->Fill(p.GetKineticEnergy());
          h_length_reconstructed->Fill(p.GetLength());
          h_hits_reconstructed->Fill(p.GetNumberOfHits());
          h_dist_reconstructed->Fill(distance_mu_p);
          h_angle_reconstructed->Fill(cos(theta_mu_p));
        }
        else {
          missed_detectable_protons++;
          h_dist_missed->Fill(distance_mu_p);
          h_angle_missed->Fill(cos(theta_mu_p));

        }
      }
      else{
        non_detectable_multiplicity++;
      }
    }
    h_detectable_protons->Fill(detectable_multiplicity);
    h_non_detectable_protons->Fill(non_detectable_multiplicity);
  }

  // Statistics
  std::cout << std::endl;
  std::cout << " Number of detectable protons     : " << detectable_protons << std::endl;
  std::cout << " Reconstructed detectable protons : " << reconstructed_detectable_protons << std::endl;
  std::cout << " Missed detectable protons        : " << missed_detectable_protons << std::endl;
  std::cout << std::endl;

  // Plots
  float scale_detectable     = h_detectable_protons->Integral();
  float scale_non_detectable = h_non_detectable_protons->Integral();

  h_detectable_protons->Scale(1./scale_detectable);
  h_non_detectable_protons->Scale(1./scale_non_detectable);
  
  float max = 1.1*std::max(h_detectable_protons->GetMaximum(), h_non_detectable_protons->GetMaximum());
  h_detectable_protons->GetYaxis()->SetRangeUser(0,max);
  h_non_detectable_protons->GetYaxis()->SetRangeUser(0,max);
  
  h_detectable_protons->SetLineColor(4);
  h_non_detectable_protons->SetLineColor(2);
  h_detectable_protons->GetXaxis()->SetTitle("Proton multiplicity");
  h_detectable_protons->GetYaxis()->SetTitle("Area Normalised");

  TCanvas *c = new TCanvas("c", "", 800, 600);
  TLegend *l = new TLegend(0.48, 0.68, 0.88, 0.88);

  l->AddEntry(h_detectable_protons,     "Proton Ek #geq 20 MeV", "l");
  l->AddEntry(h_non_detectable_protons, "Proton Ek < 20 MeV", "l");
  l->SetLineWidth(0);

  h_detectable_protons->SetStats(0);
  h_detectable_protons->Draw("hist");
  h_non_detectable_protons->Draw("hist same");
  l->Draw("same");

  c->SaveAs((plots_location+"proton_multiplicity.root").c_str());
  c->SaveAs((plots_location+"proton_multiplicity.png").c_str());

  c->Clear();
  l->Clear();

  // KE efficiency
  h_ke_reconstructed->Divide(h_ke_total);
  h_ke_reconstructed->GetXaxis()->SetTitle("Kinetic energy [GeV]");
  h_ke_reconstructed->GetYaxis()->SetTitle("Reconstruction efficiency");

  h_ke_reconstructed->SetStats(0);
  h_ke_reconstructed->Draw("hist");
  c->SaveAs((plots_location+"kinetic_energy_efficiency.root").c_str());
  c->SaveAs((plots_location+"kinetic_energy_efficiency.png").c_str());

  c->Clear();

  // Length efficiency
  h_length_reconstructed->Divide(h_length_total);
  h_length_reconstructed->GetXaxis()->SetTitle("Length [cm]");
  h_length_reconstructed->GetYaxis()->SetTitle("Reconstruction efficiency");

  h_length_reconstructed->SetStats(0);
  h_length_reconstructed->Draw("hist");
  c->SaveAs((plots_location+"length_efficiency.root").c_str());
  c->SaveAs((plots_location+"length_efficiency.png").c_str());

  c->Clear();

  // Hits efficiency
  h_hits_reconstructed->Divide(h_hits_total);
  h_hits_reconstructed->GetXaxis()->SetTitle("Hits (all planes)");
  h_hits_reconstructed->GetYaxis()->SetTitle("Reconstruction efficiency");

  h_hits_reconstructed->SetStats(0);
  h_hits_reconstructed->Draw("hist");
  c->SaveAs((plots_location+"hits_efficiency.root").c_str());
  c->SaveAs((plots_location+"hits_efficiency.png").c_str());

  c->Clear();

  // Angle missed/reconstructed
  h_angle_reconstructed->SetLineColor(4);
  h_angle_missed->SetLineColor(2);

  l->AddEntry(h_angle_reconstructed, "Reconstructed protons", "l");
  l->AddEntry(h_angle_missed, "Missed protons", "l");
  
  h_angle_reconstructed->GetXaxis()->SetTitle("cos(#theta_{#mu p})");
  h_angle_reconstructed->GetYaxis()->SetTitle("Events");

  float max_a = 1.1*std::max(h_angle_reconstructed->GetMaximum(), h_angle_missed->GetMaximum());
  h_angle_reconstructed->GetYaxis()->SetRangeUser(0,max_a);
  h_angle_missed->GetYaxis()->SetRangeUser(0,max_a);
  
  h_angle_reconstructed->SetStats(0);
  h_angle_reconstructed->Draw("hist");
  h_angle_missed->Draw("hist same");
  l->SetLineWidth(0);
  l->Draw("same");
  c->SaveAs((plots_location+"angle_distribution.root").c_str());
  c->SaveAs((plots_location+"angle_distribution.png").c_str());

  c->Clear();
  l->Clear();

  // Distance missed/reconstructed
  h_dist_reconstructed->SetLineColor(4);
  h_dist_missed->SetLineColor(2);

  l->AddEntry(h_dist_reconstructed, "Reconstructed protons", "l");
  l->AddEntry(h_dist_missed, "Missed protons", "l");
  
  h_dist_reconstructed->GetXaxis()->SetTitle("Distance 'd' [cm]");
  h_dist_reconstructed->GetYaxis()->SetTitle("Events");

  float max_d = 1.1*std::max(h_dist_reconstructed->GetMaximum(), h_dist_missed->GetMaximum());
  h_dist_reconstructed->GetYaxis()->SetRangeUser(0,max_d);
  h_dist_missed->GetYaxis()->SetRangeUser(0,max_d);
  
  h_dist_reconstructed->SetStats(0);
  h_dist_reconstructed->Draw("hist");
  h_dist_missed->Draw("hist same");
  l->SetLineWidth(0);
  l->Draw("same");
  c->SaveAs((plots_location+"dist_distribution.root").c_str());
  c->SaveAs((plots_location+"dist_distribution.png").c_str());
  
  c->Clear();
  l->Clear();

  // Distance efficiency
  h_dist_reconstructed->Divide(h_dist_total);
  h_dist_reconstructed->GetXaxis()->SetTitle("Muon-proton distance [cm]");
  h_dist_reconstructed->GetYaxis()->SetTitle("Reconstruction efficiency");

  h_dist_reconstructed->SetStats(0);
  h_dist_reconstructed->Draw("hist");
  c->SaveAs((plots_location+"dist_efficiency.root").c_str());
  c->SaveAs((plots_location+"dist_efficiency.png").c_str());

  c->Clear();

  // Angle efficiency
  h_angle_reconstructed->Divide(h_angle_total);
  h_angle_reconstructed->GetXaxis()->SetTitle("cos(#theta_{#mu p})");
  h_angle_reconstructed->GetYaxis()->SetTitle("Reconstruction efficiency");

  h_angle_reconstructed->SetStats(0);
  h_angle_reconstructed->Draw("hist");
  c->SaveAs((plots_location+"angle_efficiency.root").c_str());
  c->SaveAs((plots_location+"angle_efficiency.png").c_str());

  c->Clear();
  
  delete h_detectable_protons;
  delete h_non_detectable_protons;
  delete h_ke_reconstructed;
  delete h_ke_total;
  delete h_length_reconstructed;
  delete h_length_total;
  delete h_hits_reconstructed;
  delete h_hits_total;
  delete h_angle_reconstructed;
  delete h_angle_missed;
  delete h_angle_total;
  delete h_dist_reconstructed;
  delete h_dist_missed;
  delete h_dist_total;
  delete c;
  delete l;

  // Run timing info
  time_t rawtime_end;
  struct tm * timeinfo_end;
  time (&rawtime_end);
  timeinfo_end = localtime (&rawtime_end);
  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << " End: Local time and date:  " << asctime(timeinfo_end) << std::endl;
  std::cout << "-----------------------------------------------------------" << std::endl;

  return 0;

} // MainTest()
