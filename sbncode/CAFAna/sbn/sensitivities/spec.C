// Get Numbers for various interaction types and flavours
// cafe event_numbers.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

// #include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionExtrap.h"
#include "CAFAna/Extrap/IExtrap.h"
//#include "StandardRecord/StandardRecord.h"
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
#include "TObjString.h"

#include "CAFAna/Systs/SBNWeightSysts.h"

using namespace ana;

//const char* basicFname = "cafe_state_smear_numu.root";

//void spec(const char* stateFname = basicFname)
void spec()
{

  const char* stateFname = "cafe_state_syst_numu.root";
  const char* stateFname2 = "cafe_state_syst_nue.root";
  if (TFile(stateFname).IsZombie()){
    std:: cout << "Run make_state_syst.C first!" << std::endl;
    return;
  }
  if (TFile(stateFname2).IsZombie()){
    std:: cout << "Run make_state_syst.C(nue) first!" << std::endl;
    return;
  }

  std::cout << "Loading state from " << stateFname << std::endl; 
  std::cout << "Loading state from " << stateFname2 << std::endl; 
  TFile fin(stateFname);
  TFile fin2(stateFname2);
  PredictionInterp& pred_nd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp& pred_fd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp& pred_nd_nue  = *ana::LoadFrom<PredictionInterp>(fin2.GetDirectory("pred_nd")).release();
  PredictionInterp& pred_fd_nue  = *ana::LoadFrom<PredictionInterp>(fin2.GetDirectory("pred_fd")).release();
  PredictionInterp& pred_ub_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();
  PredictionInterp& pred_ub_nue  = *ana::LoadFrom<PredictionInterp>(fin2.GetDirectory("pred_ub")).release();
  
  std::cout << "read in done" << std::endl;


  const double pot = kPOTnominal;
  const double pot_ub = 1.32e21;

  OscCalcSterileApproxAdjustable* noosc_nd = DefaultSterileApproxCalc();
  noosc_nd->SetL(kBaselineSBND);

  //SBND
  OscCalcSterileApproxAdjustable* osc_nd_opt1 = DefaultSterileApproxCalc();
  osc_nd_opt1->SetL(kBaselineSBND);
  osc_nd_opt1->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
  osc_nd_opt1->calc.SetDmsq(1.32);

  OscCalcSterileApproxAdjustable* osc_nd_nue1 = DefaultSterileApproxCalc();
  osc_nd_nue1->SetL(kBaselineSBND);
  osc_nd_nue1->calc.SetSinSq2ThetaMuE(0.003);
  osc_nd_nue1->calc.SetDmsq(1.2);

  OscCalcSterileApproxAdjustable* osc_nd_nue2 = DefaultSterileApproxCalc();
  osc_nd_nue2->SetL(kBaselineSBND);
  osc_nd_nue2->calc.SetSinSq2ThetaMuE(0.001);
  osc_nd_nue2->calc.SetDmsq(1.32);
  
  //Numu
  TH1* hnumu_nd_signal_unosc = pred_nd_numu.PredictComponent(noosc_nd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_signal_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_ncbg_unosc = pred_nd_numu.PredictComponent(noosc_nd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_nd_ncbg_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  
  TH1* hnumu_nd_osc_ratio = (TH1*)hnumu_nd_signal_osc1->Clone();
  hnumu_nd_osc_ratio->Add(hnumu_nd_ncbg_osc1);
  TH1* hnumu_nd_tot_unosc = (TH1*)hnumu_nd_signal_unosc->Clone();
  hnumu_nd_tot_unosc->Add(hnumu_nd_ncbg_unosc);
  hnumu_nd_osc_ratio->Divide(hnumu_nd_tot_unosc);
  

  //Nue
  TH1* hnue_nd_signal_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_signal_osc1 = pred_nd_nue.PredictComponent(osc_nd_nue1, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_signal_osc2 = pred_nd_nue.PredictComponent(osc_nd_nue2, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_ncbg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_numubg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_nd_nuebg_unosc = pred_nd_nue.PredictComponent(noosc_nd, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);

  TH1* hnue_nd_osc1_ratio = (TH1*)hnue_nd_signal_osc1->Clone();
  TH1* hnue_nd_osc2_ratio = (TH1*)hnue_nd_signal_osc2->Clone();
  TH1* hnue_nd_tot_unosc = (TH1*)hnue_nd_signal_unosc->Clone();
  hnue_nd_tot_unosc->Add(hnue_nd_ncbg_unosc);
  hnue_nd_tot_unosc->Add(hnue_nd_numubg_unosc);
  hnue_nd_tot_unosc->Add(hnue_nd_nuebg_unosc);
  hnue_nd_osc1_ratio->Add(hnue_nd_ncbg_unosc);
  hnue_nd_osc1_ratio->Add(hnue_nd_numubg_unosc);
  hnue_nd_osc1_ratio->Add(hnue_nd_nuebg_unosc);
  hnue_nd_osc2_ratio->Add(hnue_nd_ncbg_unosc);
  hnue_nd_osc2_ratio->Add(hnue_nd_numubg_unosc);
  hnue_nd_osc2_ratio->Add(hnue_nd_nuebg_unosc);
  hnue_nd_osc1_ratio->Divide(hnue_nd_tot_unosc);
  hnue_nd_osc2_ratio->Divide(hnue_nd_tot_unosc);

  //for (int i = 1; i < 20; ++i) {
  //  std::cout<<hnue_nd_osc1_ratio->GetBinContent(i)<<" divided by "<<hnue_nd_tot_unosc->GetBinContent(i)<<std::endl;
  //  hnue_nd_osc1_ratio->SetBinContent(i, hnue_nd_osc1_ratio->GetBinContent(i)/hnue_nd_tot_unosc->GetBinContent(i));
  //  std::cout<<hnue_nd_osc1_ratio->GetBinContent(i)<<std::endl;
  // }

//  auto fExtrap = dynamic_cast<PredictionExtrap*>(pred_nd_nue.fPredNom.release())->GetExtrap(); 
//  TH1* hnue_nd_nuebg_fromMu_unosc = fExtrap->NueSurvFromMuComponent().Oscillated(noosc_nd, +12, +12).ToTH1(pot);
//  TH1* hnue_nd_nuebg_fromKZ_unosc = fExtrap->NueSurvFromKZComponent().Oscillated(noosc_nd, +12, +12).ToTH1(pot);
//  TH1* hnue_nd_nuebg_fromKP_unosc = fExtrap->NueSurvFromKPComponent().Oscillated(noosc_nd, +12, +12).ToTH1(pot);


  //ICARUS

  OscCalcSterileApproxAdjustable* noosc_fd = DefaultSterileApproxCalc();
  noosc_fd->SetL(kBaselineIcarus);

  OscCalcSterileApproxAdjustable* osc_fd_opt1 = DefaultSterileApproxCalc();
  osc_fd_opt1->SetL(kBaselineIcarus);
  osc_fd_opt1->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
  osc_fd_opt1->calc.SetDmsq(1.32);

  OscCalcSterileApproxAdjustable* osc_fd_nue1 = DefaultSterileApproxCalc();
  osc_fd_nue1->SetL(kBaselineIcarus);
  osc_fd_nue1->calc.SetSinSq2ThetaMuE(0.003);
  osc_fd_nue1->calc.SetDmsq(1.2);
  
  OscCalcSterileApproxAdjustable* osc_fd_nue2 = DefaultSterileApproxCalc();
  osc_fd_nue2->SetL(kBaselineIcarus);
  osc_fd_nue2->calc.SetSinSq2ThetaMuE(0.001);
  osc_fd_nue2->calc.SetDmsq(1.32);

  //Numu
  TH1* hnumu_fd_signal_unosc = pred_fd_numu.PredictComponent(noosc_fd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_signal_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_ncbg_unosc = pred_fd_numu.PredictComponent(noosc_fd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnumu_fd_ncbg_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  
  TH1* hnumu_fd_osc_ratio = (TH1*)hnumu_fd_signal_osc1->Clone();
  hnumu_fd_osc_ratio->Add(hnumu_fd_ncbg_osc1);
  TH1* hnumu_fd_tot_unosc = (TH1*)hnumu_fd_signal_unosc->Clone();
  hnumu_fd_tot_unosc->Add(hnumu_fd_ncbg_unosc);
  hnumu_fd_osc_ratio->Divide(hnumu_fd_tot_unosc);

  //Nue
  TH1* hnue_fd_signal_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_signal_osc1 = pred_fd_nue.PredictComponent(osc_fd_nue1, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_signal_osc2 = pred_fd_nue.PredictComponent(osc_fd_nue2, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_ncbg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_numubg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* hnue_fd_nuebg_unosc = pred_fd_nue.PredictComponent(noosc_fd, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  
  TH1* hnue_fd_osc1_ratio = (TH1*)hnue_fd_signal_osc1->Clone();
  TH1* hnue_fd_osc2_ratio = (TH1*)hnue_fd_signal_osc2->Clone();
  TH1* hnue_fd_tot_unosc = (TH1*)hnue_fd_signal_unosc->Clone();
  hnue_fd_tot_unosc->Add(hnue_fd_ncbg_unosc);
  hnue_fd_tot_unosc->Add(hnue_fd_numubg_unosc);
  hnue_fd_tot_unosc->Add(hnue_fd_nuebg_unosc);
  hnue_fd_osc1_ratio->Add(hnue_fd_ncbg_unosc);
  hnue_fd_osc1_ratio->Add(hnue_fd_numubg_unosc);
  hnue_fd_osc1_ratio->Add(hnue_fd_nuebg_unosc);
  hnue_fd_osc2_ratio->Add(hnue_fd_ncbg_unosc);
  hnue_fd_osc2_ratio->Add(hnue_fd_numubg_unosc);
  hnue_fd_osc2_ratio->Add(hnue_fd_nuebg_unosc);
  hnue_fd_osc1_ratio->Divide(hnue_fd_tot_unosc);
  hnue_fd_osc2_ratio->Divide(hnue_fd_tot_unosc);
  

//  auto fExtrap2 = dynamic_cast<PredictionExtrap*>(pred_fd_nue.fPredNom.release())->GetExtrap(); 
//  TH1* hnue_fd_nuebg_fromMu_unosc = fExtrap2->NueSurvFromMuComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot);
//  TH1* hnue_fd_nuebg_fromKZ_unosc = fExtrap2->NueSurvFromKZComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot);
//  TH1* hnue_fd_nuebg_fromKP_unosc = fExtrap2->NueSurvFromKPComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot);
  
  //MicroBoone

  OscCalcSterileApproxAdjustable* noosc_ub = DefaultSterileApproxCalc();
  noosc_ub->SetL(kBaselineMicroBoone);

  OscCalcSterileApproxAdjustable* osc_ub_opt1 = DefaultSterileApproxCalc();
  osc_ub_opt1->SetL(kBaselineMicroBoone);
  osc_ub_opt1->calc.SetSinSq2ThetaMuMu(4*0.135*0.135*(1-0.135*0.135));
  osc_ub_opt1->calc.SetDmsq(1.32);

  OscCalcSterileApproxAdjustable* osc_ub_nue1 = DefaultSterileApproxCalc();
  osc_ub_nue1->SetL(kBaselineMicroBoone);
  osc_ub_nue1->calc.SetSinSq2ThetaMuE(0.003);
  osc_ub_nue1->calc.SetDmsq(1.2);
  
  OscCalcSterileApproxAdjustable* osc_ub_nue2 = DefaultSterileApproxCalc();
  osc_ub_nue2->SetL(kBaselineMicroBoone);
  osc_ub_nue2->calc.SetSinSq2ThetaMuE(0.001);
  osc_ub_nue2->calc.SetDmsq(1.32);

  //Numu
  TH1* hnumu_ub_signal_unosc = pred_ub_numu.PredictComponent(noosc_ub, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnumu_ub_signal_osc1 = pred_ub_numu.PredictComponent(osc_ub_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnumu_ub_ncbg_unosc = pred_ub_numu.PredictComponent(noosc_ub, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnumu_ub_ncbg_osc1 = pred_ub_numu.PredictComponent(osc_ub_opt1, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot_ub);

  TH1* hnumu_ub_osc_ratio = (TH1*)hnumu_ub_signal_osc1->Clone();
  hnumu_ub_osc_ratio->Add(hnumu_ub_ncbg_osc1);
  TH1* hnumu_ub_tot_unosc = (TH1*)hnumu_ub_signal_unosc->Clone();
  hnumu_ub_tot_unosc->Add(hnumu_ub_ncbg_unosc);
  hnumu_ub_osc_ratio->Divide(hnumu_ub_tot_unosc);
  
  //Nue
  TH1* hnue_ub_signal_unosc = pred_ub_nue.PredictComponent(noosc_ub, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnue_ub_signal_osc1 = pred_ub_nue.PredictComponent(osc_ub_nue1, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnue_ub_signal_osc2 = pred_ub_nue.PredictComponent(osc_ub_nue2, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnue_ub_ncbg_unosc = pred_ub_nue.PredictComponent(noosc_ub, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnue_ub_numubg_unosc = pred_ub_nue.PredictComponent(noosc_ub, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot_ub);
  TH1* hnue_ub_nuebg_unosc = pred_ub_nue.PredictComponent(noosc_ub, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(pot_ub);

  TH1* hnue_ub_osc1_ratio = (TH1*)hnue_ub_signal_osc1->Clone();
  TH1* hnue_ub_osc2_ratio = (TH1*)hnue_ub_signal_osc2->Clone();
  TH1* hnue_ub_tot_unosc = (TH1*)hnue_ub_signal_unosc->Clone();
  hnue_ub_tot_unosc->Add(hnue_ub_ncbg_unosc);
  hnue_ub_tot_unosc->Add(hnue_ub_numubg_unosc);
  hnue_ub_tot_unosc->Add(hnue_ub_nuebg_unosc);
  hnue_ub_osc1_ratio->Add(hnue_ub_ncbg_unosc);
  hnue_ub_osc1_ratio->Add(hnue_ub_numubg_unosc);
  hnue_ub_osc1_ratio->Add(hnue_ub_nuebg_unosc);
  hnue_ub_osc2_ratio->Add(hnue_ub_ncbg_unosc);
  hnue_ub_osc2_ratio->Add(hnue_ub_numubg_unosc);
  hnue_ub_osc2_ratio->Add(hnue_ub_nuebg_unosc);
  hnue_ub_osc1_ratio->Divide(hnue_ub_tot_unosc);
  hnue_ub_osc2_ratio->Divide(hnue_ub_tot_unosc);
  
//  auto fExtrap3 = dynamic_cast<PredictionExtrap*>(pred_ub_nue.fPredNom.release())->GetExtrap(); 
//  TH1* hnue_ub_nuebg_fromMu_unosc = fExtrap3->NueSurvFromMuComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot_ub);
//  TH1* hnue_ub_nuebg_fromKZ_unosc = fExtrap3->NueSurvFromKZComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot_ub);
//  TH1* hnue_ub_nuebg_fromKP_unosc = fExtrap3->NueSurvFromKPComponent().Oscillated(noosc_fd, +12, +12).ToTH1(pot_ub);
  
  // Far over Near Ratios
  TH1* hnumu_ub_nd_ratio = (TH1*)hnumu_ub_signal_osc1->Clone();
  hnumu_ub_nd_ratio->Add(hnumu_ub_ncbg_osc1);
  TH1* hnumu_nd_tot_osc = (TH1*)hnumu_nd_signal_osc1->Clone();
  hnumu_nd_tot_osc->Add(hnumu_nd_ncbg_osc1);
  hnumu_ub_nd_ratio->Divide(hnumu_nd_tot_osc);

  TH1* hnumu_fd_nd_ratio = (TH1*)hnumu_fd_signal_osc1->Clone();
  hnumu_fd_nd_ratio->Add(hnumu_fd_ncbg_osc1);
  hnumu_fd_nd_ratio->Divide(hnumu_nd_tot_osc);

  TH1* hnue_ub_nd_ratio1 = (TH1*)hnue_ub_signal_osc1->Clone();
  hnue_ub_nd_ratio1->Add(hnue_ub_numubg_unosc);
  hnue_ub_nd_ratio1->Add(hnue_ub_ncbg_unosc);
  hnue_ub_nd_ratio1->Add(hnue_ub_nuebg_unosc);
  TH1* hnue_nd_tot_nue1 = (TH1*)hnue_nd_signal_osc1->Clone();
  hnue_nd_tot_nue1->Add(hnue_nd_numubg_unosc);
  hnue_nd_tot_nue1->Add(hnue_nd_ncbg_unosc);
  hnue_nd_tot_nue1->Add(hnue_nd_nuebg_unosc);
  hnue_ub_nd_ratio1->Divide(hnue_nd_tot_nue1);
  TH1* hnue_ub_nd_ratio2 = (TH1*)hnue_ub_signal_osc2->Clone();
  hnue_ub_nd_ratio2->Add(hnue_ub_numubg_unosc);
  hnue_ub_nd_ratio2->Add(hnue_ub_ncbg_unosc);
  hnue_ub_nd_ratio2->Add(hnue_ub_nuebg_unosc);
  TH1* hnue_nd_tot_nue2 = (TH1*)hnue_nd_signal_osc2->Clone();
  hnue_nd_tot_nue2->Add(hnue_nd_numubg_unosc);
  hnue_nd_tot_nue2->Add(hnue_nd_ncbg_unosc);
  hnue_nd_tot_nue2->Add(hnue_nd_nuebg_unosc);
  hnue_ub_nd_ratio2->Divide(hnue_nd_tot_nue2);

  TH1* hnue_fd_nd_ratio1 = (TH1*)hnue_fd_signal_osc1->Clone();
  hnue_fd_nd_ratio1->Add(hnue_fd_numubg_unosc);
  hnue_fd_nd_ratio1->Add(hnue_fd_ncbg_unosc);
  hnue_fd_nd_ratio1->Add(hnue_fd_nuebg_unosc);
  hnue_fd_nd_ratio1->Divide(hnue_nd_tot_nue1);
  TH1* hnue_fd_nd_ratio2 = (TH1*)hnue_fd_signal_osc2->Clone();
  hnue_fd_nd_ratio2->Add(hnue_fd_numubg_unosc);
  hnue_fd_nd_ratio2->Add(hnue_fd_ncbg_unosc);
  hnue_fd_nd_ratio2->Add(hnue_fd_nuebg_unosc);
  hnue_fd_nd_ratio2->Divide(hnue_nd_tot_nue2);

  std::cout << "Writing File" << std::endl;

  TFile* fOutput = new TFile("output_spec.root","RECREATE");

  hnumu_nd_signal_unosc->Write("hnumu_nd_signal_unosc");
  hnumu_nd_signal_osc1->Write("hnumu_nd_signal_osc");
  hnumu_nd_ncbg_unosc->Write("hnumu_nd_ncbg_unosc");
  hnumu_nd_ncbg_osc1->Write("hnumu_nd_ncbg_osc");
  hnumu_nd_osc_ratio->Write("hnumu_nd_osc_ratio");  

  hnumu_fd_signal_unosc->Write("hnumu_fd_signal_unosc");
  hnumu_fd_signal_osc1->Write("hnumu_fd_signal_osc");
  hnumu_fd_ncbg_unosc->Write("hnumu_fd_ncbg_unosc");
  hnumu_fd_ncbg_osc1->Write("hnumu_fd_ncbg_osc");
  hnumu_fd_osc_ratio->Write("hnumu_fd_osc_ratio");
  
  hnumu_ub_signal_unosc->Write("hnumu_ub_signal_unosc");
  hnumu_ub_signal_osc1->Write("hnumu_ub_signal_osc");
  hnumu_ub_ncbg_unosc->Write("hnumu_ub_ncbg_unosc");
  hnumu_ub_ncbg_osc1->Write("hnumu_ub_ncbg_osc");
  hnumu_ub_osc_ratio->Write("hnumu_ub_osc_ratio");

  hnumu_ub_nd_ratio->Write("hnumu_ub_nd_ratio");
  hnumu_fd_nd_ratio->Write("hnumu_fd_nd_ratio");

  hnue_nd_signal_unosc->Write("hnue_nd_signal_unosc");
  hnue_nd_signal_osc1->Write("hnue_nd_signal_osc1");
  hnue_nd_signal_osc2->Write("hnue_nd_signal_osc2");
  hnue_nd_ncbg_unosc->Write("hnue_nd_ncbg_unosc");
  hnue_nd_numubg_unosc->Write("hnue_nd_numubg_unosc");
  hnue_nd_nuebg_unosc->Write("hnue_nd_nuebg_unosc");
//  hnue_nd_nuebg_fromMu_unosc->Write("hnue_nd_nuebg_fromMu_unosc");
//  hnue_nd_nuebg_fromKZ_unosc->Write("hnue_nd_nuebg_fromKZ_unosc");
//  hnue_nd_nuebg_fromKP_unosc->Write("hnue_nd_nuebg_fromKP_unosc");
  hnue_nd_osc1_ratio->Write("hnue_nd_osc1_ratio");
  hnue_nd_osc2_ratio->Write("hnue_nd_osc2_ratio");

  hnue_fd_signal_unosc->Write("hnue_fd_signal_unosc");
  hnue_fd_signal_osc1->Write("hnue_fd_signal_osc1");
  hnue_fd_signal_osc2->Write("hnue_fd_signal_osc2");
  hnue_fd_ncbg_unosc->Write("hnue_fd_ncbg_unosc");
  hnue_fd_numubg_unosc->Write("hnue_fd_numubg_unosc");
  hnue_fd_nuebg_unosc->Write("hnue_fd_nuebg_unosc");
//  hnue_fd_nuebg_fromMu_unosc->Write("hnue_fd_nuebg_fromMu_unosc");
//  hnue_fd_nuebg_fromKZ_unosc->Write("hnue_fd_nuebg_fromKZ_unosc");
//  hnue_fd_nuebg_fromKP_unosc->Write("hnue_fd_nuebg_fromKP_unosc");
  hnue_fd_osc1_ratio->Write("hnue_fd_osc1_ratio");
  hnue_fd_osc2_ratio->Write("hnue_fd_osc2_ratio");

  hnue_ub_signal_unosc->Write("hnue_ub_signal_unosc");
  hnue_ub_signal_osc1->Write("hnue_ub_signal_osc1");
  hnue_ub_signal_osc2->Write("hnue_ub_signal_osc2");
  hnue_ub_ncbg_unosc->Write("hnue_ub_ncbg_unosc");
  hnue_ub_numubg_unosc->Write("hnue_ub_numubg_unosc");
  hnue_ub_nuebg_unosc->Write("hnue_ub_nuebg_unosc");
//  hnue_ub_nuebg_fromMu_unosc->Write("hnue_ub_nuebg_fromMu_unosc");
//  hnue_ub_nuebg_fromKZ_unosc->Write("hnue_ub_nuebg_fromKZ_unosc");
//  hnue_ub_nuebg_fromKP_unosc->Write("hnue_ub_nuebg_fromKP_unosc");
  hnue_ub_osc1_ratio->Write("hnue_ub_osc1_ratio");
  hnue_ub_osc2_ratio->Write("hnue_ub_osc2_ratio");

  hnue_ub_nd_ratio1->Write("hnue_ub_nd_ratio1");
  hnue_ub_nd_ratio2->Write("hnue_ub_nd_ratio2");
  hnue_fd_nd_ratio1->Write("hnue_fd_nd_ratio1");
  hnue_fd_nd_ratio2->Write("hnue_fd_nd_ratio2");

  fOutput->Close();

}
