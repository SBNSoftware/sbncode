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
#include "TROOT.h"
#include "TAxis.h"

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
  std::string plots_location  = "../Output_Selection_Tool/plots/vertices/";

  //------------------------------------------------------------------------------------------
  //                                       Load events
  //------------------------------------------------------------------------------------------
  
  // Initialise event list and the topology maps
  EventSelectionTool::EventList events;
  
  // Maps
  TopologyMap ccinc = GeneralAnalysisHelper::GetCCIncTopologyMap();
  //TopologyMap cc0pi_map = GeneralAnalysisHelper::GetCC0PiTopologyMap();
  //TopologyMap cc1pi_map = GeneralAnalysisHelper::GetCC1PiTopologyMap();
 
  // Initiate histograms
  /*
  TH1D *h_truedirection_mu = new TH1D("h_truedirection_mu","", 50, -1, 1);
  TH1D *h_truedirection_pi = new TH1D("h_truedirection_pi","", 50, -1, 1);

  TH1D *h_truemomentum_mu = new TH1D("h_truemomentum_mu","", 50, 0, 2);
  TH1D *h_truemomentum_pi = new TH1D("h_truemomentum_pi","", 50, 0, 2);
  */
  TH1D *h_vertexdiff_x_sig = new TH1D("h_vertexdiff_x_sig","", 100, -100, 100);
  TH1D *h_vertexdiff_y_sig = new TH1D("h_vertexdiff_y_sig","", 100, -100, 100);
  TH1D *h_vertexdiff_z_sig = new TH1D("h_vertexdiff_z_sig","", 100, -100, 100);

  TH1D *h_vertexdiff_x_bac = new TH1D("h_vertexdiff_x_bac","", 100, -100, 100);
  TH1D *h_vertexdiff_y_bac = new TH1D("h_vertexdiff_y_bac","", 100, -100, 100);
  TH1D *h_vertexdiff_z_bac = new TH1D("h_vertexdiff_z_bac","", 100, -100, 100);

  // Counters
  unsigned int signal           = 0;
  unsigned int delta_x_10_sig   = 0;
  unsigned int delta_y_10_sig   = 0;
  unsigned int delta_z_10_sig   = 0;
  unsigned int delta_x_20_sig   = 0;
  unsigned int delta_y_20_sig   = 0;
  unsigned int delta_z_20_sig   = 0;
  unsigned int delta_xyz_10_sig = 0;
  unsigned int delta_xyz_20_sig = 0;
  
  unsigned int background       = 0;
  unsigned int delta_x_10_bac   = 0;
  unsigned int delta_y_10_bac   = 0;
  unsigned int delta_z_10_bac   = 0;
  unsigned int delta_x_20_bac   = 0;
  unsigned int delta_y_20_bac   = 0;
  unsigned int delta_z_20_bac   = 0;
  unsigned int delta_xyz_10_bac = 0;
  unsigned int delta_xyz_20_bac = 0;

  int start = static_cast<int>(time(NULL));
  unsigned int total_files = 20;

  // Load the events into the event list
  for( unsigned int i = 0; i < total_files; ++i ){
    // Get the filenames
    std::string name;
    name.clear();
    char file_name[1024];
    name = "/home/rhiannon/Samples/LocalSamples/analysis/200219_neutrino_only/selection/"+std::to_string(i)+"/output_file.root";
//    name = "/home/rhiannon/Samples/LiverpoolSamples/120918_analysis_sample/11509725_"+std::to_string(i)+"/output_file.root";
    strcpy( file_name, name.c_str() );

    EventSelectionTool::LoadEventList(file_name, events, i);
    EventSelectionTool::GetTimeLeft(start,total_files,i);
  }
  std::cout << std::endl;
  
  // Loop over events and perform vertexing study
  for(const Event &e : events){

    // Particles 
//    ParticleList mc = e.GetMCParticleList();

    // TPC track criteria
    if(!GeneralAnalysisHelper::MaxOneEscapingTrack(e)) continue;

    // Neutrino vertex is within the ficudial border
    if(e.IsSBNDTrueFiducial()){
      if(e.CheckRecoTopology(ccinc)){
        double xdiff = e.GetMCNuVertex()[0] - e.GetRecoNuVertex()[0];
        double ydiff = e.GetMCNuVertex()[1] - e.GetRecoNuVertex()[1];
        double zdiff = e.GetMCNuVertex()[2] - e.GetRecoNuVertex()[2];
        if(e.CheckMCTopology(ccinc)){
          signal++;
          h_vertexdiff_x_sig->Fill(xdiff);
          h_vertexdiff_y_sig->Fill(ydiff);
          h_vertexdiff_z_sig->Fill(zdiff);
          if(std::fabs(xdiff) < 10) delta_x_10_sig++;
          if(std::fabs(ydiff) < 10) delta_y_10_sig++;
          if(std::fabs(zdiff) < 10) delta_z_10_sig++;
          if(std::fabs(xdiff) < 20) delta_x_20_sig++;
          if(std::fabs(ydiff) < 20) delta_y_20_sig++;
          if(std::fabs(zdiff) < 20) delta_z_20_sig++;
          if(std::fabs(xdiff) < 10 && ydiff < 10 && zdiff < 10) delta_xyz_10_sig++;
          if(std::fabs(xdiff) < 20 && ydiff < 20 && zdiff < 20) delta_xyz_20_sig++;
        }
        else{
          background++;
          h_vertexdiff_x_bac->Fill(xdiff);
          h_vertexdiff_y_bac->Fill(ydiff);
          h_vertexdiff_z_bac->Fill(zdiff);
          if(std::fabs(xdiff) < 10) delta_x_10_bac++;
          if(std::fabs(ydiff) < 10) delta_y_10_bac++;
          if(std::fabs(zdiff) < 10) delta_z_10_bac++;
          if(std::fabs(xdiff) < 20) delta_x_20_bac++;
          if(std::fabs(ydiff) < 20) delta_y_20_bac++;
          if(std::fabs(zdiff) < 20) delta_z_20_bac++;
          if(std::fabs(xdiff) < 10 && ydiff < 10 && zdiff < 10) delta_xyz_10_bac++;
          if(std::fabs(xdiff) < 20 && ydiff < 20 && zdiff < 20) delta_xyz_20_bac++;
        }
      }

      /*
      // True muon-pion comparison of z-angle
      for(const Particle &p : mc){
        if(p.GetPdgCode() == 13){
          h_truedirection_mu->Fill(p.GetCosTheta());
          h_truemomentum_mu->Fill(p.GetModulusMomentum());
        }
        if(p.GetPdgCode() == 211 || p.GetPdgCode() == -211){
          h_truedirection_pi->Fill(p.GetCosTheta());
          h_truemomentum_pi->Fill(p.GetModulusMomentum());
        }
      }
      */

      // Get the reconstructed vertex position in the numu CCInc. final state and compare it to the true vertex position in x,y,z
    }
  }
  // Length vs dist 2D plot
  TCanvas *c = new TCanvas("c","",800,600);

  h_vertexdiff_x_sig->SetLineColor(kSpring-3);
  h_vertexdiff_x_bac->SetLineColor(kOrange+7);
  h_vertexdiff_y_sig->SetLineColor(kSpring-3);
  h_vertexdiff_y_bac->SetLineColor(kOrange+7);
  h_vertexdiff_z_sig->SetLineColor(kSpring-3);
  h_vertexdiff_z_bac->SetLineColor(kOrange+7);
  
  // Legend
  TLegend *l      = new TLegend(0.58,0.68,0.88,0.88);
  l->AddEntry(h_vertexdiff_x_sig,  "#nu_{#mu} CC Inc. Signal",  "l");
  l->AddEntry(h_vertexdiff_x_bac,  "#nu_{#mu} CC Inc. Background",  "l");
  l->SetLineWidth(0);
  l->SetTextAlign(12);
  l->SetTextFont(133);

  // Max y
  double maxx = 1.1 * std::max(h_vertexdiff_x_sig->GetMaximum(), h_vertexdiff_x_bac->GetMaximum());
  std::cout << " Max x : " << maxx << ", sig max : " << h_vertexdiff_x_sig->GetMaximum() << ", bac max : " << h_vertexdiff_x_bac->GetMaximum() << std::endl;

  h_vertexdiff_x_sig->SetStats(0);
  h_vertexdiff_x_sig->GetXaxis()->SetTitle("Vertex x-position difference: MC - Reco");
  h_vertexdiff_x_sig->GetYaxis()->SetTitle("Event rate");
  h_vertexdiff_x_sig->GetXaxis()->SetTitleFont(132);
  h_vertexdiff_x_sig->GetXaxis()->SetLabelFont(132);
  h_vertexdiff_x_sig->GetYaxis()->SetTitleFont(132);
  h_vertexdiff_x_sig->GetYaxis()->SetLabelFont(132);
  h_vertexdiff_x_sig->GetYaxis()->SetTitleOffset(1.2);
  h_vertexdiff_x_sig->GetYaxis()->SetMaxDigits(3);
  h_vertexdiff_x_sig->GetYaxis()->SetRangeUser(0, maxx);
  h_vertexdiff_x_sig->Draw();
  h_vertexdiff_x_bac->Draw("same");
  l->Draw("same");

  c->SaveAs((plots_location+"vertex_x_diff_ccinc.root").c_str());
  c->Clear();
  l->Clear();

  l->AddEntry(h_vertexdiff_y_sig,  "#nu_{#mu} CC Inc. Signal",  "l");
  l->AddEntry(h_vertexdiff_y_bac,  "#nu_{#mu} CC Inc. Background",  "l");
  l->SetLineWidth(0);
  l->SetTextAlign(12);
  l->SetTextFont(133);

  // Max y
  double maxy = 1.1 * std::max(h_vertexdiff_y_sig->GetMaximum(), h_vertexdiff_y_bac->GetMaximum());

  h_vertexdiff_y_sig->SetStats(0);
  h_vertexdiff_y_sig->GetXaxis()->SetTitle("Vertex y-position difference: MC - Reco");
  h_vertexdiff_y_sig->GetYaxis()->SetTitle("Event rate");
  h_vertexdiff_y_sig->GetXaxis()->SetTitleFont(132);
  h_vertexdiff_y_sig->GetXaxis()->SetLabelFont(132);
  h_vertexdiff_y_sig->GetYaxis()->SetTitleFont(132);
  h_vertexdiff_y_sig->GetYaxis()->SetLabelFont(132);
  h_vertexdiff_y_sig->GetYaxis()->SetTitleOffset(1.2);
  h_vertexdiff_y_sig->GetYaxis()->SetMaxDigits(3);
  h_vertexdiff_y_sig->GetYaxis()->SetRangeUser(0, maxy);
  h_vertexdiff_y_sig->Draw();
  h_vertexdiff_y_bac->Draw("same");
  l->Draw("same");

  c->SaveAs((plots_location+"vertex_y_diff_ccinc.root").c_str());
  c->Clear();
  l->Clear();

  l->AddEntry(h_vertexdiff_z_sig,  "#nu_{#mu} CC Inc. Signal",  "l");
  l->AddEntry(h_vertexdiff_z_bac,  "#nu_{#mu} CC Inc. Background",  "l");
  l->SetLineWidth(0);
  l->SetTextAlign(12);
  l->SetTextFont(133);

  // Max y
  double maxz = 1.1 * std::max(h_vertexdiff_z_sig->GetMaximum(), h_vertexdiff_z_bac->GetMaximum());

  h_vertexdiff_z_sig->SetStats(0);
  h_vertexdiff_z_sig->GetXaxis()->SetTitle("Vertex z-position difference: MC - Reco");
  h_vertexdiff_z_sig->GetYaxis()->SetTitle("Event rate");
  h_vertexdiff_z_sig->GetXaxis()->SetTitleFont(132);
  h_vertexdiff_z_sig->GetXaxis()->SetLabelFont(132);
  h_vertexdiff_z_sig->GetYaxis()->SetTitleFont(132);
  h_vertexdiff_z_sig->GetYaxis()->SetLabelFont(132);
  h_vertexdiff_z_sig->GetYaxis()->SetTitleOffset(1.2);
  h_vertexdiff_z_sig->GetYaxis()->SetMaxDigits(3);
  h_vertexdiff_z_sig->GetYaxis()->SetRangeUser(0, maxz);
  h_vertexdiff_z_sig->Draw();
  h_vertexdiff_z_bac->Draw("same");
  l->Draw("same");

  c->SaveAs((plots_location+"vertex_z_diff_ccinc.root").c_str());
  c->Clear();
  l->Clear();

  std::cout << " Signal : " << signal << std::endl;
  std::cout << "   |MC - Reco| < 10 cm : " << std::endl;
  std::cout << "                     x : " << delta_x_10_sig << std::endl;
  std::cout << "                     y : " << delta_y_10_sig << std::endl;
  std::cout << "                     z : " << delta_z_10_sig << std::endl;
  std::cout << "             x & y & z : " << delta_xyz_10_sig << std::endl;
  std::cout << "   |MC - Reco| < 20 cm : " << std::endl;
  std::cout << "                     x : " << delta_x_20_sig << std::endl;
  std::cout << "                     y : " << delta_y_20_sig << std::endl;
  std::cout << "                     z : " << delta_z_20_sig << std::endl;
  std::cout << "             x & y & z : " << delta_xyz_20_sig << std::endl;
  std::cout << " Background : " << background << std::endl;
  std::cout << "   |MC - Reco| < 10 cm : " << std::endl;
  std::cout << "                     x : " << delta_x_10_bac << std::endl;
  std::cout << "                     y : " << delta_y_10_bac << std::endl;
  std::cout << "                     z : " << delta_z_10_bac << std::endl;
  std::cout << "             x & y & z : " << delta_xyz_10_bac << std::endl;
  std::cout << "   |MC - Reco| < 20 cm : " << std::endl;
  std::cout << "                     x : " << delta_x_20_bac << std::endl;
  std::cout << "                     y : " << delta_y_20_bac << std::endl;
  std::cout << "                     z : " << delta_z_20_bac << std::endl;
  std::cout << "             x & y & z : " << delta_xyz_20_bac << std::endl;


  /*
  h_truedirection_mu->Scale(1./double(h_truedirection_mu->Integral()));
  h_truedirection_pi->Scale(1./double(h_truedirection_pi->Integral()));
  h_truedirection_mu->SetFillColor(kSpring-3);
  h_truedirection_pi->SetFillColor(kOrange+7);
  h_truedirection_mu->SetFillStyle(3001);
  h_truedirection_pi->SetFillStyle(3001);
  h_truedirection_mu->SetLineColor(kSpring-3);
  h_truedirection_pi->SetLineColor(kOrange+7);

  h_truemomentum_mu->Scale(1./double(h_truemomentum_mu->Integral()));
  h_truemomentum_pi->Scale(1./double(h_truemomentum_pi->Integral()));
  h_truemomentum_mu->SetFillColor(kSpring-3);
  h_truemomentum_pi->SetFillColor(kOrange+7);
  h_truemomentum_mu->SetFillStyle(3001);
  h_truemomentum_pi->SetFillStyle(3001);
  h_truemomentum_mu->SetLineColor(kSpring-3);
  h_truemomentum_pi->SetLineColor(kOrange+7);

  TCanvas *c      = new TCanvas("c", "", 800, 600);
  TLegend *l      = new TLegend(0.22,0.68,0.52,0.88);
  
  // Legend
  l->AddEntry(h_truedirection_mu,  "#mu direction",  "f");
  l->AddEntry(h_truedirection_pi,  "#pi direction",  "f");
  l->SetLineWidth(0);
  l->SetTextAlign(12);
  l->SetTextFont(133);

  // Max y
  double max_direction = 1.1 * std::max(h_truedirection_mu->GetMaximum(), h_truedirection_pi->GetMaximum());

  // Draw plots
  h_truedirection_mu->SetStats(0);
  h_truedirection_mu->GetXaxis()->SetTitle("Track angle, cos(#theta)");
  h_truedirection_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_truedirection_mu->GetXaxis()->SetTitleFont(132);
  h_truedirection_mu->GetXaxis()->SetLabelFont(132);
  h_truedirection_mu->GetYaxis()->SetTitleFont(132);
  h_truedirection_mu->GetYaxis()->SetLabelFont(132);
  h_truedirection_mu->GetYaxis()->SetTitleOffset(1.2);
  h_truedirection_mu->GetYaxis()->SetMaxDigits(3);
  h_truedirection_mu->GetYaxis()->SetRangeUser(0, max_direction);
  h_truedirection_mu->Draw("hist");
  h_truedirection_pi->Draw("same hist");
  l->Draw("same");

  c->SaveAs((plots_location+"track_direction.root").c_str());
  c->Clear();
  l->Clear();

  // Legend
  TLegend *l_mom = new TLegend(0.58,0.68,0.88,0.88);
  l_mom->AddEntry(h_truemomentum_mu,  "#mu momentum",  "f");
  l_mom->AddEntry(h_truemomentum_pi,  "#pi momentum",  "f");
  l_mom->SetLineWidth(0);
  l_mom->SetTextAlign(12);
  l_mom->SetTextFont(133);

  // Max y
  double max_momentum = 1.1 * std::max(h_truemomentum_mu->GetMaximum(), h_truemomentum_pi->GetMaximum());

  // Draw plots
  h_truemomentum_mu->SetStats(0);
  h_truemomentum_mu->GetXaxis()->SetTitle("Track momentum, [GeV]");
  h_truemomentum_mu->GetYaxis()->SetTitle("Normalised event rate");
  h_truemomentum_mu->GetXaxis()->SetTitleFont(132);
  h_truemomentum_mu->GetXaxis()->SetLabelFont(132);
  h_truemomentum_mu->GetYaxis()->SetTitleFont(132);
  h_truemomentum_mu->GetYaxis()->SetLabelFont(132);
  h_truemomentum_mu->GetYaxis()->SetTitleOffset(1.2);
  h_truemomentum_mu->GetYaxis()->SetMaxDigits(3);
  h_truemomentum_mu->GetYaxis()->SetRangeUser(0, max_momentum);
  h_truemomentum_mu->Draw("hist");
  h_truemomentum_pi->Draw("same hist");
  l_mom->Draw("same");

  c->SaveAs((plots_location+"track_momentum.root").c_str());
  c->Clear();
  l_mom->Clear();

  */
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

} // MainTest()
