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
#include "OscLib/func/OscCalculatorSterile.h"
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

const char* basicFname = "cafe_state_smear_numu.root";

void spec(const char* stateFname = basicFname)
{
  if (TFile(stateFname).IsZombie()){
    std:: cout << "Run make_state.C first!" << std::endl;
    return;
  }
  else{

    std::cout << "Loading state from " << stateFname << std::endl; 
    TFile fin(stateFname);
    PredictionInterp& pred_nd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd_numu")).release();
    PredictionInterp& pred_fd_numu = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu")).release();
    //PredictionInterp& pred_fd_nue = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue")).release();

    std::cout << "read in done" << std::endl;


    const double pot = POTnominal;

    osc::NoOscillations* osc_none = new osc::NoOscillations;

    //SBND
    osc::OscCalculatorSterile* osc_nd_opt1 = DefaultSterileCalc(4);
    osc_nd_opt1->SetL(BaselineSBND);
    osc_nd_opt1->SetAngle(2, 4, 0.1);
    osc_nd_opt1->SetDm(4, 0.44);

    osc::OscCalculatorSterile* osc_nd_opt2 = DefaultSterileCalc(4);
    osc_nd_opt2->SetL(BaselineSBND);
    osc_nd_opt2->SetAngle(2, 4, 0.1);
    osc_nd_opt2->SetDm(4, 1.1);

    //Numu
    TH1* hnumu_nd_signal_unosc = pred_nd_numu.PredictComponent(osc_none, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_nd_signal_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_nd_signal_osc2 = pred_nd_numu.PredictComponent(osc_nd_opt2, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_nd_ncbg_unosc = pred_nd_numu.PredictComponent(osc_none, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
    //TH1* hnumu_nd_ncbg_osc1 = pred_nd_numu.PredictComponent(osc_nd_opt1, Flavors::kAllNuMu, Current::kNC, Sign::kBoth).ToTH1(pot);
    //TH1* hnumu_nd_ncbg_osc2 = pred_nd_numu.PredictComponent(osc_nd_opt2, Flavors::kAllNuMu, Current::kNC, Sign::kBoth).ToTH1(pot);


    //ICARUS

    osc::OscCalculatorSterile* osc_fd_opt1 = DefaultSterileCalc(4);
    osc_fd_opt1->SetL(BaselineIcarus);
    osc_fd_opt1->SetAngle(2, 4, 0.1);
    osc_fd_opt1->SetDm(4, 0.44);

    osc::OscCalculatorSterile* osc_fd_opt2 = DefaultSterileCalc(4);
    osc_fd_opt2->SetL(BaselineIcarus);
    osc_fd_opt2->SetAngle(2, 4, 0.1);
    osc_fd_opt2->SetDm(4, 1.1);

    //osc::OscCalculatorSterile* osc_fd_nue = DefaultSterileCalc(4);
    //osc_fd_nue->SetL(BaselineIcarus);
    //osc_fd_nue->SetDm(4, 0.43);
    //uhoh - make new calculator?
    //osc_fd_nue->SetAngle 

    //Numu
    TH1* hnumu_fd_signal_unosc = pred_fd_numu.PredictComponent(osc_none, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_fd_signal_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_fd_signal_osc2 = pred_fd_numu.PredictComponent(osc_fd_opt2, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
    TH1* hnumu_fd_ncbg_unosc = pred_fd_numu.PredictComponent(osc_none, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);
    //TH1* hnumu_fd_ncbg_osc1 = pred_fd_numu.PredictComponent(osc_fd_opt1, Flavors::kAllNuMu, Current::kNC, Sign::kBoth).ToTH1(pot);
    //TH1* hnumu_fd_ncbg_osc2 = pred_fd_numu.PredictComponent(osc_fd_opt2, Flavors::kAllNuMu, Current::kNC, Sign::kBoth).ToTH1(pot);

    //Nue
    //TH1* hnue_fd_signal_unosc = pred_fd_nue.PredictComponent(osc_fd_nue, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
    //TH1* hnue_fd_signal_osc = pred_fd_nue.PredictComponent(osc_fd_nue, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(pot);

    TFile* fOutput = new TFile("output_spec.root","RECREATE");

    hnumu_nd_signal_unosc->Write("hnumu_nd_signal_unosc");
    hnumu_nd_signal_osc1->Write("hnumu_nd_signal_osc1");
    hnumu_nd_signal_osc2->Write("hnumu_nd_signal_osc2");
    hnumu_nd_ncbg_unosc->Write("hnumu_nd_ncbg_unosc");
    //hnumu_nd_ncbg_osc1->Write("hnumu_nd_ncbg_osc1");
    //hnumu_nd_ncbg_osc2->Write("hnumu_nd_ncbg_osc2");

    hnumu_fd_signal_unosc->Write("hnumu_fd_signal_unosc");
    hnumu_fd_signal_osc1->Write("hnumu_fd_signal_osc1");
    hnumu_fd_signal_osc2->Write("hnumu_fd_signal_osc2");
    hnumu_fd_ncbg_unosc->Write("hnumu_fd_ncbg_unosc");
    //hnumu_fd_ncbg_osc1->Write("hnumu_fd_ncbg_osc1");
    //hnumu_fd_ncbg_osc2->Write("hnumu_fd_ncbg_osc2");

    fOutput->Close();
  }
}
