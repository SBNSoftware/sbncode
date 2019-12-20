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
#include "mySysts.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"

#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void syst_fit(const std::string anatype = numuStr)
{
  //  GetSBNWeightSysts(); // initialize
  const std::vector<const ISyst*>& wtsysts = GetSBNWeightSysts();

  std::vector<const ISyst*> systs_to_process;

  std::vector<std::string> syst_names{"expskin_FluxUnisim","horncurrent_FluxUnisim","kminus_PrimaryHadronNormalization","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piminus_PrimaryHadronSWCentralSplineVariation","pioninexsec_FluxUnisim","pionqexsec_FluxUnisim","piontotxsec_FluxUnisim","piplus_PrimaryHadronSWCentralSplineVariation","genie_ccresAxial_Genie","genie_ncresAxial_Genie","genie_qema_Genie","genie_NC_Genie","genie_NonResRvbarp1pi_Genie","genie_NonResRvbarp2pi_Genie","genie_NonResRvp1pi_Genie","genie_NonResRvp2pi_Genie","genie_NonResRvbarp1piAlt_Genie","genie_NonResRvbarp2piAlt_Genie","genie_NonResRvp1piAlt_Genie","genie_NonResRvp2piAlt_Genie"};

  for (auto s : wtsysts) {
    for (auto n : syst_names) if (n == s->ShortName()) systs_to_process.push_back(s);
  }
  
  std::vector<const ISyst*> systs;
  //for (const ISyst* s : det_systs) systs.push_back(s);
  for (const ISyst* s : systs_to_process) systs.push_back(s);

  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "output_syst_fit_numu.root";
  }
  else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "output_syst_fit_nue.root";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }
  
  TFile fin(name_in);
  
  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();
  
  //std::vector<const ISyst*> bigsysts;
  
  //osc::NoOscillations unosc;
  //for(const ISyst* s: systs){
  //if(fabs(p_fd->PredictSyst(&unosc, SystShifts(s, +1)).Integral(1e20)/p_fd->Predict(&unosc).Integral(1e20)-1) > .01) bigsysts.push_back(s);
  //}
  
  //std::cout << bigsysts.size() << " big systs out of " << systs.size() << std::endl;
  //for(const ISyst* s: bigsysts) std::cout << s->ShortName() << " ";
  //std::cout << std::endl;
  
  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  
  calc->SetL(kBaselineSBND);
  const Spectrum data_nd = p_nd->Predict(calc).FakeData(sbndPOT);
  calc->SetL(kBaselineIcarus);
  const Spectrum data_fd = p_fd->Predict(calc).FakeData(icarusPOT);
  calc->SetL(kBaselineMicroBoone);
  const Spectrum data_ub = p_ub->Predict(calc).FakeData(uboonePOT);
  
  SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  SingleSampleExperiment expt_ub(p_ub, data_ub);
   
  MultiExperimentSBN multiExpt({&expt_nd, &expt_fd, &expt_ub}, {kSBND, kICARUS, kMicroBoone});
  
   std::vector<std::vector<const ISyst*>> slists;
   //slists.push_back(bigsysts);
   slists.push_back(systs);
   //for(const ISyst* s: systs) slists.emplace_back(1, s); // and then each

   std::vector<const IFitVar*> oscVars;
   if (anatype == numuStr) {   
     oscVars = {&kFitDmSqSterile, &kFitSinSq2ThetaMuMu};
   }
   else {
     oscVars = {&kFitDmSqSterile, &kFitSinSq2ThetaMuE};
   }

   for(const std::vector<const ISyst*> slist: slists){

     std::cout << "Fitting all 3 expts together..." << std::endl;
     Fitter fit_syst(&multiExpt, oscVars, slist, Fitter::kCareful);
     OscCalcSterileApproxAdjustable* calc_multi = DefaultSterileApproxCalc();
     SystShifts bestSysts_multi = kNoShift;
     double chi_multi = fit_syst.Fit(calc_multi, bestSysts_multi);

     std::vector<std::string> pnames = fit_syst.GetParamNames();
     std::vector<double> prefit = fit_syst.GetPreFitValues();
     std::vector<double> prefit_err = fit_syst.GetPreFitErrors();
     std::vector<double> postfit = fit_syst.GetPostFitValues();
     std::vector<double> postfit_err = fit_syst.GetPostFitErrors();
  
     // std::cout << "Parameters considered: " << std::endl;
     // for (const string st: pnames) {
     //   std::cout << st << std::endl;
     // }
     // std::cout << "Prefit params: " << std::endl;
     // for(const double d: prefit) {
     //   std::cout << d << std::endl;
     // }
     // std::cout << "Postfit params: " << std::endl;
     // for(const double d: postfit) {
     //   std::cout << d << std::endl;
     // }
     // std::cout << "Prefit param errs: " << std::endl;
     // for(const double d: prefit_err) {
     //   std::cout << d << std::endl;
     // }
     // std::cout << "Postfit param errs: " << std::endl;
     // for(const double d: postfit_err) {
     //   std::cout << d << std::endl;
     // }

     int nsyst = systs.size();
     int nvar = oscVars.size();
     TH1D *hsyst_pre = new TH1D("hsyst_pre","hsyst_pre",nsyst, 0, nsyst);
     TH1D *hsyst_post = new TH1D("hsyst_post","hsyst_post",nsyst, 0, nsyst);
     for (int i=0;i<nsyst+nvar;i++) {
       if (i>nvar-1) {
	 hsyst_pre->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	 hsyst_pre->Fill(i-nvar,prefit_err[i]);
	 hsyst_post->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	 hsyst_post->Fill(i-nvar,postfit_err[i]);
       }
     }

     std::cout << "Fitting ND only..." << std::endl;
     Fitter fit_nd(&expt_nd, oscVars, slist, Fitter::kCareful);
     OscCalcSterileApproxAdjustable* calc_nd = DefaultSterileApproxCalc();
     calc_nd->SetL(kBaselineSBND);
     SystShifts bestSysts_nd;
     double chi_nd = fit_nd.Fit(calc_nd, bestSysts_nd);
     std::vector<double> prefitnd = fit_nd.GetPreFitValues();
     std::vector<double> prefitnd_err = fit_nd.GetPreFitErrors();
     std::vector<double> postfitnd = fit_nd.GetPostFitValues();
     std::vector<double> postfitnd_err = fit_nd.GetPostFitErrors();
     TH1D *hsystnd_pre = new TH1D("hsystnd_pre","hsystnd_pre",nsyst-nvar, 0, nsyst-nvar);
     TH1D *hsystnd_post = new TH1D("hsystnd_post","hsystnd_post",nsyst-nvar, 0, nsyst-nvar);
     for (int i=0;i<nsyst;i++) {
      if (i>=nvar) {
	hsystnd_pre->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	hsystnd_pre->Fill(i-nvar,prefitnd_err[i]);
	hsystnd_post->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	hsystnd_post->Fill(i-nvar,postfitnd_err[i]);
      }
     }

     std::cout << "Fitting FD only..." << std::endl;
     Fitter fit_fd(&expt_fd, oscVars, slist, Fitter::kCareful);
     OscCalcSterileApproxAdjustable* calc_fd = DefaultSterileApproxCalc();
     calc_fd->SetL(kBaselineIcarus);
     SystShifts bestSysts_fd;
     double chi_fd = fit_fd.Fit(calc_fd, bestSysts_fd);
     std::vector<double> prefitfd = fit_fd.GetPreFitValues();
     std::vector<double> prefitfd_err = fit_fd.GetPreFitErrors();
     std::vector<double> postfitfd = fit_fd.GetPostFitValues();
     std::vector<double> postfitfd_err = fit_fd.GetPostFitErrors();
     TH1D *hsystfd_pre = new TH1D("hsystfd_pre","hsystfd_pre",nsyst-nvar, 0, nsyst-nvar);
     TH1D *hsystfd_post = new TH1D("hsystfd_post","hsystfd_post",nsyst-nvar, 0, nsyst-nvar);
     for (int i=0;i<nsyst;i++) {
       if (i>=nvar) {
	 hsystfd_pre->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	 hsystfd_pre->Fill(i-nvar,prefitfd_err[i]);
	 hsystfd_post->GetXaxis()->SetBinLabel(i+1-nvar,pnames[i].c_str());
	 hsystfd_post->Fill(i-nvar,postfitfd_err[i]);
       }
     }


     TFile *outfile = new TFile(name_out,"RECREATE");
     hsyst_pre->Write();
     hsyst_post->Write();
     hsystnd_pre->Write();
     hsystnd_post->Write();
     hsystfd_pre->Write();
     hsystfd_post->Write();

     outfile->Close();

   } // end for s
}
