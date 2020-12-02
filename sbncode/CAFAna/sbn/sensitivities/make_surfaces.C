#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/SBNFluxSysts.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Systs/UniverseOracle.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"

#include <vector>
#include <fstream>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_surfaces(const std::string anatype = numuStr)
{
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  //auto systs_flux = GetSBNFluxWeightSysts();
  //auto systs_genie = GetSBNGenieWeightSysts();

  GetSBNFluxHadronSysts(30);
  GetMECSyst();
 
  std::vector<const ISyst*> systs_to_process;

  // Proposal Systs
  std::vector<std::string> syst_names{"expskin_FluxUnisim","horncurrent_FluxUnisim","kminus_PrimaryHadronNormalization","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piminus_PrimaryHadronSWCentralSplineVariation","pioninexsec_FluxUnisim","pionqexsec_FluxUnisim","piontotxsec_FluxUnisim","piplus_PrimaryHadronSWCentralSplineVariation","genie_ccresAxial_Genie","genie_ncresAxial_Genie","genie_qema_Genie","genie_NC_Genie","genie_NonResRvbarp1pi_Genie","genie_NonResRvbarp2pi_Genie","genie_NonResRvp1pi_Genie","genie_NonResRvp2pi_Genie","genie_NonResRvbarp1piAlt_Genie","genie_NonResRvbarp2piAlt_Genie","genie_NonResRvp1piAlt_Genie","genie_NonResRvp2piAlt_Genie"};

  // Select only proposal systs + MEC
  for (auto s : systs) {
    for (auto n : syst_names) if (n == s->ShortName()) systs_to_process.push_back(s);
  }

  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs_to_process.push_back(s);

  systs_to_process.push_back(&GetMECSyst());

  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "surfaces_numu.root";
  }
  else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "surfaces_nue.root";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);
  TFile fout(name_out,"RECREATE");

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();

  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  OscCalcSterileApproxAdjustable* seed = DefaultSterileApproxCalc();

  if (anatype == nueStr) {
    seed->calc.SetSinSq2ThetaMuE(0.001);
    seed->calc.SetDmsq(1.32);
    p_nd->SetOscSeed(seed);
    p_fd->SetOscSeed(seed);
    p_ub->SetOscSeed(seed);
  }
    
  //Define fit axes
  //smaller fit axes for no
  const IFitVar* SterileTh;
  double thlim = 1e-3;
  if (anatype == numuStr) {
    SterileTh = &kFitSinSq2ThetaMuMu;
  }
  else {
    SterileTh = &kFitSinSq2ThetaMuE;
    thlim = 5e-5;
  }
  const FitAxis kAxForTh(SterileTh, 40, thlim, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 1e-2, 1e2, true);

  // We'll call zero nominal
  const Spectrum data_nd = p_nd->Predict(calc).FakeData(sbndPOT);
  const Spectrum data_fd = p_fd->Predict(calc).FakeData(icarusPOT);
  const Spectrum data_ub = p_ub->Predict(calc).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  SingleSampleExperiment expt_ub(p_ub, data_ub);

  MultiExperimentSBN multiExpt({&expt_nd, &expt_fd, &expt_ub}, {kSBND, kICARUS, kMicroBoone});
  MultiExperimentSBN fd_nd({&expt_nd, &expt_fd}, {kSBND, kICARUS});

  Surface surf_nom(&multiExpt, calc, kAxForTh, kAxDmSq);
  Surface surf_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_nd(&expt_nd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_fd(&expt_fd, calc, kAxForTh, kAxDmSq);
  Surface surf_nom_ub(&expt_ub, calc, kAxForTh, kAxDmSq);

  fout.mkdir("exclusion");
  fout.cd("exclusion");    

  surf_nom.SaveTo(gDirectory->mkdir("nom"));
  surf_nom_nd.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd.SaveTo(gDirectory->mkdir("nom_fd"));
  surf_nom_ub.SaveTo(gDirectory->mkdir("nom_ub"));
  surf_nd_fd.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  std::map<std::string, std::vector<const ISyst*>> slists;
  slists["prop_systs_mec"] = systs_to_process;
  //slists[systs[0]->ShortName()] = {systs[0]};

  for(auto syst_pair: slists){
    std::string suffix = syst_pair.first;
    std::vector<const ISyst*> slist = syst_pair.second;

    Surface surf_syst(&multiExpt, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_nd_fd(&fd_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_nd(&expt_nd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_fd(&expt_fd, calc, kAxForTh, kAxDmSq, {}, slist);
    Surface surf_syst_ub(&expt_ub, calc, kAxForTh, kAxDmSq, {}, slist);

    surf_syst_nd.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
    surf_syst_ub.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
    surf_syst.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
    surf_syst_nd_fd.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));

  } // end for s


  // Allowed Region

  OscCalcSterileApproxAdjustable* calc2 = DefaultSterileApproxCalc(); // Global Best Fit
  OscCalcSterileApproxAdjustable* calc3 = DefaultSterileApproxCalc(); // LSND (nue only)
  if (anatype == numuStr) {
    calc2->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
    calc2->calc.SetDmsq(1.32);
  }
  else {
    calc2->calc.SetSinSq2ThetaMuE(0.001);
    calc2->calc.SetDmsq(1.32);
    calc3->calc.SetSinSq2ThetaMuE(0.003);
    calc3->calc.SetDmsq(1.2);
  }

  const Spectrum data_nd2 = p_nd->Predict(calc2).FakeData(sbndPOT);
  const Spectrum data_fd2 = p_fd->Predict(calc2).FakeData(icarusPOT);
  const Spectrum data_ub2 = p_ub->Predict(calc2).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd2(p_nd, data_nd2);
  SingleSampleExperiment expt_fd2(p_fd, data_fd2);
  SingleSampleExperiment expt_ub2(p_ub, data_ub2);

  MultiExperimentSBN multiExpt2({&expt_nd2, &expt_fd2, &expt_ub2}, {kSBND, kICARUS, kMicroBoone});
  MultiExperimentSBN fd_nd2({&expt_nd2, &expt_fd2}, {kSBND, kICARUS});

  double mlim_hi;
  double mlim_lo;
  if (anatype == numuStr) {
    thlim = 1e-2;
    mlim_hi = 1e1;
    mlim_lo = 1e-1;
  }
  else {
    thlim = 5e-5;
    mlim_hi = 1e2;
    mlim_lo = 1e-2;
  }
  const FitAxis kAxForTh2(SterileTh, 40, thlim, 1, true);
  const FitAxis kAxDmSq2(&kFitDmSqSterile, 40, mlim_lo, mlim_hi, true);

  Surface surf_nom2(&multiExpt2, calc2, kAxForTh2, kAxDmSq2);
  Surface surf_nd_fd2(&fd_nd2, calc2, kAxForTh2, kAxDmSq2);
  Surface surf_nom_nd2(&expt_nd2, calc2, kAxForTh2, kAxDmSq2);
  Surface surf_nom_fd2(&expt_fd2, calc2, kAxForTh2, kAxDmSq2);
  Surface surf_nom_ub2(&expt_ub2, calc2, kAxForTh2, kAxDmSq2);
    
  fout.cd("..");
  fout.mkdir("allowed");
  fout.cd("allowed");

  if(anatype == nueStr) {
    fout.mkdir("glb_fit");
    fout.cd("glb_fit");
  }

  surf_nom2.SaveTo(gDirectory->mkdir("nom"));
  surf_nom_nd2.SaveTo(gDirectory->mkdir("nom_nd"));
  surf_nom_fd2.SaveTo(gDirectory->mkdir("nom_fd"));
  surf_nom_ub2.SaveTo(gDirectory->mkdir("nom_ub"));
  surf_nd_fd2.SaveTo(gDirectory->mkdir("nom_nd_fd"));

  for(auto syst_pair: slists){
    std::string suffix = syst_pair.first;
    std::vector<const ISyst*> slist = syst_pair.second;

    Surface surf_syst2(&multiExpt2, calc2, kAxForTh2, kAxDmSq2, {}, slist);
    Surface surf_syst_nd_fd2(&fd_nd2, calc2, kAxForTh2, kAxDmSq2, {}, slist);
    Surface surf_syst_nd2(&expt_nd2, calc2, kAxForTh2, kAxDmSq2, {}, slist);
    Surface surf_syst_fd2(&expt_fd2, calc2, kAxForTh2, kAxDmSq2, {}, slist);
    Surface surf_syst_ub2(&expt_ub2, calc2, kAxForTh2, kAxDmSq2, {}, slist);

    surf_syst_nd2.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
    surf_syst_fd2.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
    surf_syst_ub2.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
    surf_syst2.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
    surf_syst_nd_fd2.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));

  } // end for s
  

  if (anatype == nueStr) {  
    const Spectrum data_nd3 = p_nd->Predict(calc3).FakeData(sbndPOT);
    const Spectrum data_fd3 = p_fd->Predict(calc3).FakeData(icarusPOT);
    const Spectrum data_ub3 = p_ub->Predict(calc3).FakeData(uboonePOT);
  
    SingleSampleExperiment expt_nd3(p_nd, data_nd3);
    SingleSampleExperiment expt_fd3(p_fd, data_fd3);
    SingleSampleExperiment expt_ub3(p_ub, data_ub3);
  
    MultiExperimentSBN multiExpt3({&expt_nd3, &expt_fd3, &expt_ub3}, {kSBND, kICARUS, kMicroBoone});
    MultiExperimentSBN fd_nd3({&expt_nd3, &expt_fd3}, {kSBND, kICARUS});
  
    Surface surf_nom3(&multiExpt3, calc3, kAxForTh2, kAxDmSq2);
    Surface surf_nd_fd3(&fd_nd3, calc3, kAxForTh2, kAxDmSq2);
    Surface surf_nom_nd3(&expt_nd3, calc3, kAxForTh2, kAxDmSq2);
    Surface surf_nom_fd3(&expt_fd3, calc3, kAxForTh2, kAxDmSq2);
    Surface surf_nom_ub3(&expt_ub3, calc3, kAxForTh2, kAxDmSq2);
      
    fout.cd("..");
    fout.mkdir("LSND");
    fout.cd("LSND");
  
    surf_nom3.SaveTo(gDirectory->mkdir("nom"));
    surf_nom_nd3.SaveTo(gDirectory->mkdir("nom_nd"));
    surf_nom_fd3.SaveTo(gDirectory->mkdir("nom_fd"));
    surf_nom_ub3.SaveTo(gDirectory->mkdir("nom_ub"));
    surf_nd_fd3.SaveTo(gDirectory->mkdir("nom_nd_fd"));
  
    for(auto syst_pair: slists){
      std::string suffix = syst_pair.first;
      std::vector<const ISyst*> slist = syst_pair.second;

      Surface surf_syst3(&multiExpt3, calc3, kAxForTh2, kAxDmSq2, {}, slist);
      Surface surf_syst_nd_fd3(&fd_nd3, calc3, kAxForTh2, kAxDmSq2, {}, slist);
      Surface surf_syst_nd3(&expt_nd3, calc3, kAxForTh2, kAxDmSq2, {}, slist);
      Surface surf_syst_fd3(&expt_fd3, calc3, kAxForTh2, kAxDmSq2, {}, slist);
      Surface surf_syst_ub3(&expt_ub3, calc3, kAxForTh2, kAxDmSq2, {}, slist);
  
      surf_syst_nd3.SaveTo(gDirectory->mkdir(("nd_"+suffix).c_str()));
      surf_syst_fd3.SaveTo(gDirectory->mkdir(("fd_"+suffix).c_str()));
      surf_syst_ub3.SaveTo(gDirectory->mkdir(("ub_"+suffix).c_str()));
      surf_syst3.SaveTo(gDirectory->mkdir(("allexpt_"+suffix).c_str()));
      surf_syst_nd_fd3.SaveTo(gDirectory->mkdir(("nd_fd_"+suffix).c_str()));
  
    } // end for s
  }
}
