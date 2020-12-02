#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Systs/SBNWeightSysts.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TGraph.h"

#include <fstream>
#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void allowed(const std::string anatype = numuStr)
{
  const char* name_in;
  const char* name_out;
  const char* graph_out;
  std::string dir;
  if (anatype == numuStr) {
    name_in = "surfaces_numu.root";
    name_out = "allowed_numu.pdf";
    graph_out = "CAFAna_numu_disapp_allowed.root";
    dir = "allowed/";
  }
  else if (anatype == nueStr) {
    name_in = "surfaces_nue.root";
    name_out = "allowed_nue.pdf";
    graph_out = "CAFAna_nue_app_allowed.root";
    dir = "LSND/";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  std::string suffix = "prop_systs_mec";

  TFile fin(name_in);

  Surface& surf_syst_nd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nd_"+suffix).c_str())).release(); 
  Surface& surf_syst_fd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"fd_"+suffix).c_str())).release(); 
  Surface& surf_syst_ub = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"ub_"+suffix).c_str())).release(); 
  Surface& surf_syst_nd_fd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nd_fd_"+suffix).c_str())).release();
  Surface& surf_syst_all = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"allexpt_"+suffix).c_str())).release(); 
  
  Surface& surf_nom_nd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nom_nd").c_str())).release(); 
  Surface& surf_nom_ub = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nom_ub").c_str())).release(); 
  Surface& surf_nom_fd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nom_fd").c_str())).release(); 
  Surface& surf_nom_nd_fd = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nom_nd_fd").c_str())).release(); 
  Surface& surf_nom = *ana::LoadFrom<Surface>(fin.GetDirectory((dir+"nom").c_str())).release(); 

  TH2* crit5sig = Gaussian5Sigma2D(surf_nom);
  TH2* crit3sig = Gaussian3Sigma2D(surf_nom);
  TH2* crit90 = Gaussian90Percent2D(surf_nom);
  TH2* crit95 = Gaussian95Percent2D(surf_nom);
  TH2* crit99 = Gaussian99Percent2D(surf_nom);

  surf_nom.SetTitle("5#sigma Allowed Regions");
 
  surf_nom.DrawContour(crit5sig, 7, kBlack);
  surf_nom_nd_fd.DrawContour(crit5sig, 7, kMagenta);
  surf_nom_nd.DrawContour(crit5sig,7,kRed);
  surf_nom_ub.DrawContour(crit5sig,7,kGreen);
  surf_nom_fd.DrawContour(crit5sig, 7, kBlue);

  surf_syst_all.DrawContour(crit5sig, kSolid, kBlack);
  surf_syst_nd.DrawContour(crit5sig, kSolid, kRed);
  surf_syst_ub.DrawContour(crit5sig, kSolid, kGreen);
  surf_syst_fd.DrawContour(crit5sig, kSolid, kBlue);
  surf_syst_nd_fd.DrawContour(crit5sig, kSolid, kMagenta);
    
  TLegend * lgdis = new TLegend(0.11,0.11,0.3,0.3);
  lgdis->SetFillColor(0);
  lgdis->SetBorderSize(0);
  TH1* dummy = new TH1F("", "", 1, 0, 1);
  dummy->SetLineColor(kRed);
  lgdis->AddEntry(dummy->Clone(), "SBND");
  dummy->SetLineColor(kGreen);
  lgdis->AddEntry(dummy->Clone(), "MicroBoone");
  dummy->SetLineColor(kBlue);
  lgdis->AddEntry(dummy->Clone(), "Icarus");
  dummy->SetLineColor(kMagenta);
  lgdis->AddEntry(dummy->Clone(), "SBND + Icarus");
  dummy->SetLineColor(kBlack);
  lgdis->AddEntry(dummy->Clone(), "All");
  dummy->SetLineStyle(7);
  lgdis->AddEntry(dummy->Clone(), "Stats Only");

  lgdis->Draw("same");
  gPad->Print(name_out);

  TFile fout(graph_out, "RECREATE");
  
  std::vector<Surface> surfaces {surf_nom, surf_nom_nd, surf_nom_ub, surf_nom_fd, surf_nom_nd_fd, surf_syst_all, surf_syst_nd, surf_syst_ub, surf_syst_fd, surf_syst_nd_fd};

  std::vector<TH2*> con_lev {crit90, crit95, crit99, crit3sig, crit5sig};

  std::vector<std::string> s_name {"all_stat", "nd_stat", "ub_stat", "fd_stat", "nd_fd_stat", "all_syst", "nd_syst", "ub_syst", "fd_syst", "nd_fd_syst"};

  std::vector<std::string> c_name {"_90pct", "_95pct", "_99pct", "_3sig", "_5sig"};
  
  for (int i = 0; i < 10; ++i) {
    for (int j = 0; j < 5; ++j) {
      std::vector<TGraph*> graphs = surfaces[i].GetGraphs(con_lev[j]);
      int k = 1;
      for (auto g : graphs) {
        g->Write((s_name[i]+c_name[j]+"_"+std::to_string(k)).c_str());
        ++k;
      }
    }
  } 
  
}
