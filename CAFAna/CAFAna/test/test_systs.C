#include "CAFAna/Systs/GenieSysts.h"

#include "CAFAna/Analysis/TDRLoaders.h"
#include "CAFAna/Analysis/CalcsNuFit.h"

#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"

#include "CAFAna/Core/Ratio.h"

#include "CAFAna/Cuts/AnaCuts.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEFluxSysts.h"
#include "CAFAna/Systs/EnergySysts.h"

#include "OscLib/func/IOscCalculator.h"

using namespace ana;

#include "TCanvas.h"

const Binning binsFDEreco = Binning::Simple(80, 0, 10);
const Var kRecoE_numu = SIMPLEVAR(dune.Ev_reco_numu);
const HistAxis axRecoEnuFDnumu("Reco energy (GeV)", binsFDEreco, kRecoE_numu);

const Var kGENIEWeights = SIMPLEVAR(dune.total_cv_wgt);

std::vector<const ISyst*> GetListOfSysts()
{

  std::vector<const ISyst*> systlist;

  std::vector<const ISyst*> fluxlist = GetDUNEFluxSysts(10, true);
  systlist.insert(systlist.end(), fluxlist.begin(), fluxlist.end());

  std::vector<const ISyst*> detlist_dis  = {&keScaleMuLArSyst, &kChargedHadCorrSyst,
                                            &kChargedHadUncorrFDSyst, &kNUncorrFDSyst,
                                            &kEnergyScalePi0Syst, &kPi0UncorrFDSyst};

  systlist.insert(systlist.end(), detlist_dis.begin(), detlist_dis.end());

  std::vector<const ISyst*> xseclist = GetGenieSysts({}, true);
  systlist.insert(systlist.end(), xseclist.begin(), xseclist.end());

  return systlist;
};

void test_systs(bool reload = false)
{
  osc::IOscCalculatorAdjustable* calc = NuFitOscCalc(1);

  const std::vector<const ISyst*> systs = GetListOfSysts();

  if(reload){
    TDRLoaders loaders(Loaders::kFHC);

    NoExtrapPredictionGenerator gen(axRecoEnuFDnumu, kPassFD_CVN_NUMU && kIsTrueFV, kGENIEWeights);

    PredictionInterp pred(systs, calc, gen, loaders);

    loaders.Go();

    SaveToFile(pred, "test_systs.root", "pred");
  }

  IPrediction* pred = LoadFromFile<IPrediction>("test_systs.root", "pred").release();

  Spectrum nom = pred->Predict(calc);

  new TCanvas;
  gPad->Print("test_systs.pdf[");

  for(const ISyst* s: systs){
    Spectrum up = pred->PredictSyst(calc, SystShifts(s, +1));
    Spectrum dn = pred->PredictSyst(calc, SystShifts(s, -1));

    Ratio rup(up, nom);
    Ratio rdn(dn, nom);

    TH1* hup = rup.ToTH1(kRed);
    hup->Draw("hist");
    hup->SetTitle(s->LatexName().c_str());
    hup->GetYaxis()->SetRangeUser(.8, 1.2);
    rdn.ToTH1(kBlue)->Draw("hist same");

    gPad->Print("test_systs.pdf");
  }

  gPad->Print("test_systs.pdf]");
}
