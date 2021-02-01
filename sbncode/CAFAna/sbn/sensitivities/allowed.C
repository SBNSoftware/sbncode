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

#include "OscLib/IOscCalc.h"

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

void allowed()
{
  TFile fin("surfaces.root");  

  Surface& surf_syst_nd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nd_prop_systs")).release(); 
  Surface& surf_syst_fd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/fd_prop_systs")).release(); 
  Surface& surf_syst_ub = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/ub_prop_systs")).release(); 
  Surface& surf_syst_nd_fd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nd_fd_prop_systs")).release();
  Surface& surf_syst_all = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/allexpt_prop_systs")).release(); 
  
  Surface& surf_nom_nd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nom_nd")).release(); 
  Surface& surf_nom_ub = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nom_ub")).release(); 
  Surface& surf_nom_fd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nom_fd")).release(); 
  Surface& surf_nom_nd_fd = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nom_nd_fd")).release(); 
  Surface& surf_nom = *ana::LoadFrom<Surface>(fin.GetDirectory("allowed/nom")).release(); 

  TH2* crit5sig = Gaussian5Sigma2D(surf_nom);
  TH2* crit3sig = Gaussian3Sigma2D(surf_nom);
  TH2* crit90 = Gaussian90Percent2D(surf_nom);
  TH2* crit95 = Gaussian95Percent2D(surf_nom);
  TH2* crit99 = Gaussian99Percent2D(surf_nom);

  surf_nom.SetTitle("3 Sigma Allowed Regions");
 
  surf_nom.DrawContour(crit3sig, 7, kBlack);
  surf_nom_nd_fd.DrawContour(crit3sig, 7, kMagenta);
  surf_nom_nd.DrawContour(crit3sig,7,kRed);
  surf_nom_ub.DrawContour(crit3sig,7,kGreen);
  surf_nom_fd.DrawContour(crit3sig, 7, kBlue);

  surf_syst_all.DrawContour(crit3sig, kSolid, kBlack);
  surf_syst_nd.DrawContour(crit3sig, kSolid, kRed);
  surf_syst_ub.DrawContour(crit3sig, kSolid, kGreen);
  surf_syst_fd.DrawContour(crit3sig, kSolid, kBlue);
  surf_syst_nd_fd.DrawContour(crit3sig, kSolid, kMagenta);
    
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
  gPad->Print("allowed_3sig.pdf");

//  TFile fout("allowed_graphs.root", "RECREATE");
//  
//  std::vector<Surface> surfaces {surf_nom, surf_nom_nd, surf_nom_ub, surf_nom_fd, surf_nom_nd_fd, surf_syst_all, surf_syst_nd, surf_syst_ub, surf_syst_fd, surf_syst_nd_fd};
//
//  std::vector<TH2*> con_lev {crit90, crit95, crit99, crit3sig, crit5sig};
//
//  std::vector<std::string> s_name {"all_stat", "nd_stat", "ub_stat", "fd_stat", "nd_fd_stat", "all_syst", "nd_syst", "ub_syst", "fd_syst", "nd_fd_syst"};
//
//  std::vector<std::string> c_name {"_90pct", "_95pct", "_99pct", "_3sig", "_5sig"};
//  
//  std::cout<<"Before Loop"<<std::endl;
//
//  ofstream txtfile;
//  txtfile.open("CAFAna_allowed.txt");
// 
//  for (int i = 0; i < 10; ++i) {
//    txtfile<<s_name[i]<<"\n"<<"=========="<<"\n";
//    for (int j = 0; j < 5; ++j) {
//      std::vector<TGraph*> graphs = surfaces[i].GetGraphs(con_lev[j]);
//      std::cout<<"Got Graph"<<std::endl;
//      txtfile<<c_name[j]<<"\n"<<"------"<<"\n";
//      for (auto g : graphs) {
//        g->Write((s_name[i]+c_name[j]).c_str());
//        int n = g->GetN();
//        double *x = g->GetX();
//        double *y = g->GetY();
//        for (int k = 0; k < n; ++k) {
//	  //x_v.push_back(x[k]);
//	  //y_v.push_back(y[k]);
//          txtfile<<x[k]<<" "<<y[k]<<"\n";
//	}
//      }
//      txtfile<<"\n";
//    }
//    txtfile<<"\n\n";
//  } 
//
//  txtfile.close();
  
}
