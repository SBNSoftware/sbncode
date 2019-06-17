// Get Numbers for various interaction types and flavours
// cafe event_numbers.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

// #include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Core/MultiVar.h"
#include "OscLib/OscCalculatorSterile.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"

#include "toysysts.h"

using namespace ana;

//const char* basicFname = "cafe_state_smear_numu.root";

//void spec(const char* stateFname = basicFname)
void spec()
{

  const char* stateFname = "cafe_state_smear_numu.root";
  const char* stateFname2 = "cafe_state_smear_nue.root";
  if (TFile(stateFname).IsZombie()){
    std:: cout << "Run make_state.C first!" << std::endl;
    return;
  }
  if (TFile(stateFname2).IsZombie()){
    std:: cout << "Run make_state.C(nue) first!" << std::endl;
    return;
  }

  std::cout << "Loading state from " << stateFname << std::endl; 
  std::cout << "Loading state from " << stateFname2 << std::endl; 
  TFile fin(stateFname);
  TFile fin2(stateFname2);
  PredictionInterp& pred_nd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd_numu")).release();
  PredictionInterp& pred_fd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu")).release();
  PredictionInterp& pred_nd_nue = *ana::LoadFrom<PredictionInterp>(fin2.GetDirectory("pred_nd_nue")).release();
  PredictionInterp& pred_fd_nue = *ana::LoadFrom<PredictionInterp>(fin2.GetDirectory("pred_fd_nue")).release();
  
  std::cout << "read in done" << std::endl;


  const double pot = kPOTnominal;

  OscCalcSterileApproxAdjustable* noosc_nd = DefaultSterileApproxCalc();
  noosc_nd->SetL(kBaselineSBND);

  //SBND
  OscCalcSterileApproxAdjustable* osc_nd_opt1 = DefaultSterileApproxCalc();
  osc_nd_opt1->SetL(kBaselineSBND);
  osc_nd_opt1->calc.SetSinSq2ThetaMuMu(0.1);
  osc_nd_opt1->calc.SetDmsq(0.44);

  OscCalcSterileApproxAdjustable* osc_nd_opt2 = DefaultSterileApproxCalc();
  osc_nd_opt2->SetL(kBaselineSBND);
  osc_nd_opt2->calc.SetSinSq2ThetaMuMu(0.1);
  osc_nd_opt2->calc.SetDmsq(1.1);

  OscCalcSterileApproxAdjustable* osc_nd_nue1 = DefaultSterileApproxCalc();
  osc_nd_nue1->SetL(kBaselineSBND);
  osc_nd_nue1->calc.SetSinSq2ThetaMuE(0.013);
  osc_nd_nue1->calc.SetDmsq(0.43);

  //Numu
  TH1* hnumu_nd_signal_unosc = pred_nd_numu.PredictComponent(noosc_nd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_signal_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_signal_osc2 = pred_nd_numu.PredictComponent(osc_nd_opt2, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_ncbg_unosc = pred_nd_numu.PredictComponent(noosc_nd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_ncbg_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_ncbg_osc2 = pred_nd_numu.PredictComponent(osc_nd_opt2, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);

  //Nue
  TH1* hnue_nd_signal_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_signal_osc1 = pred_nd_nue.PredictComponent(osc_nd_nue1, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_ncbg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_numubg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_nuebg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);


  //ICARUS

  OscCalcSterileApproxAdjustable* noosc_fd = DefaultSterileApproxCalc();
  noosc_fd->SetL(kBaselineIcarus);

  OscCalcSterileApproxAdjustable* osc_fd_opt1 = DefaultSterileApproxCalc();
  osc_fd_opt1->SetL(kBaselineIcarus);
  osc_fd_opt1->calc.SetSinSq2ThetaMuMu(0.1);
  osc_fd_opt1->calc.SetDmsq(0.44);

  OscCalcSterileApproxAdjustable* osc_fd_opt2 = DefaultSterileApproxCalc();
  osc_fd_opt2->SetL(kBaselineIcarus);
  osc_fd_opt2->calc.SetSinSq2ThetaMuMu(0.1);
  osc_fd_opt2->calc.SetDmsq(1.1);

  OscCalcSterileApproxAdjustable* osc_fd_nue1 = DefaultSterileApproxCalc();
  osc_fd_nue1->SetL(kBaselineIcarus);
  osc_fd_nue1->calc.SetSinSq2ThetaMuE(0.013);
  osc_fd_nue1->calc.SetDmsq(0.43);
  
  //Numu
  TH1* hnumu_fd_signal_unosc = pred_fd_numu.PredictComponent(noosc_fd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_signal_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_signal_osc2 = pred_fd_numu.PredictComponent(osc_fd_opt2, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_ncbg_unosc = pred_fd_numu.PredictComponent(noosc_fd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_ncbg_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_ncbg_osc2 = pred_fd_numu.PredictComponent(osc_fd_opt2, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);

  //Nue
  TH1* hnue_fd_signal_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_signal_osc1 = pred_fd_nue.PredictComponent(osc_fd_nue1, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_ncbg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_numubg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_nuebg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);

  TFile* fOutput = new TFile("output_spec.root","RECREATE");

  hnumu_nd_signal_unosc->Write("hnumu_nd_signal_unosc");
  hnumu_nd_signal_osc1->Write("hnumu_nd_signal_osc1");
  hnumu_nd_signal_osc2->Write("hnumu_nd_signal_osc2");
  hnumu_nd_ncbg_unosc->Write("hnumu_nd_ncbg_unosc");
  hnumu_nd_ncbg_osc1->Write("hnumu_nd_ncbg_osc1");
  hnumu_nd_ncbg_osc2->Write("hnumu_nd_ncbg_osc2");

  hnumu_fd_signal_unosc->Write("hnumu_fd_signal_unosc");
  hnumu_fd_signal_osc1->Write("hnumu_fd_signal_osc1");
  hnumu_fd_signal_osc2->Write("hnumu_fd_signal_osc2");
  hnumu_fd_ncbg_unosc->Write("hnumu_fd_ncbg_unosc");
  hnumu_fd_ncbg_osc1->Write("hnumu_fd_ncbg_osc1");
  hnumu_fd_ncbg_osc2->Write("hnumu_fd_ncbg_osc2");

  hnue_nd_signal_unosc->Write("hnue_nd_signal_unosc");
  hnue_nd_signal_osc1->Write("hnue_nd_signal_osc1");
  hnue_nd_ncbg_unosc->Write("hnue_nd_ncbg_unosc");
  hnue_nd_numubg_unosc->Write("hnue_nd_numubg_unosc");
  hnue_nd_nuebg_unosc->Write("hnue_nd_nuebg_unosc");

  hnue_fd_signal_unosc->Write("hnue_fd_signal_unosc");
  hnue_fd_signal_osc1->Write("hnue_fd_signal_osc1");
  hnue_fd_ncbg_unosc->Write("hnue_fd_ncbg_unosc");
  hnue_fd_numubg_unosc->Write("hnue_fd_numubg_unosc");
  hnue_fd_nuebg_unosc->Write("hnue_fd_nuebg_unosc");

  fOutput->Close();

}
