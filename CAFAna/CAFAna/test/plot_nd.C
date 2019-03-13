#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Analysis/Fit.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"

void plot_nd()
{
  rootlogon(); // style

  SpectrumLoader loaderFHC("/dune/data/users/marshalc/NDTF_FGT_FHC.root");
  SpectrumLoader loaderRHC("/dune/data/users/marshalc/NDTF_FGT_RHC.root");

  // TODO - this whole dance isn't really necessary for the saner ND
  // files. Work on SpectrumLoader so that it can just count all the POT in
  // this case without the intermediaries.
  auto* loaderFHCPOT = loaderFHC.LoaderForRunPOT(1);
  auto* loaderRHCPOT = loaderRHC.LoaderForRunPOT(1);

  const Var kRecoE = SIMPLEVAR(dune.Ev_reco);
  const Var kPIDmu = SIMPLEVAR(dune.numu_pid);
  const Var kPIDe = SIMPLEVAR(dune.nue_pid);
  const Var kQ = SIMPLEVAR(dune.reco_q);

  const HistAxis axis("Reconstructed energy (GeV)",
                      Binning::Simple(40, 0, 10),
                      kRecoE);

  const HistAxis axisPIDmu("PIDmu",
                         Binning::Simple(30, -1.5, +1.5),
                         kPIDmu);
  const HistAxis axisPIDe("PIDe",
                         Binning::Simple(30, -1.5, +1.5),
                         kPIDe);

  Spectrum sMuPIDFHCNumu(*loaderFHCPOT, axisPIDmu, kIsNumuCC);
  Spectrum sMuPIDFHCNC(*loaderFHCPOT, axisPIDmu, kIsNC);
  Spectrum sMuPIDFHCNue(*loaderFHCPOT, axisPIDmu, kIsBeamNue);

  Spectrum sEPIDFHCNumu(*loaderFHCPOT, axisPIDe, kIsNumuCC);
  Spectrum sEPIDFHCNC(*loaderFHCPOT, axisPIDe, kIsNC);
  Spectrum sEPIDFHCNue(*loaderFHCPOT, axisPIDe, kIsBeamNue);
  
  Spectrum sRecoEFHCNumuSelTot(*loaderFHCPOT, axis, kPIDmu > 0.5 && kQ < 0.);
  Spectrum sRecoEFHCNumuSelNumu(*loaderFHCPOT, axis, kIsNumuCC && !kIsAntiNu && kPIDmu > 0.5 && kQ < 0.);
  Spectrum sRecoEFHCNumuSelNumubar(*loaderFHCPOT, axis, kIsNumuCC && kIsAntiNu && kPIDmu > 0.5 && kQ < 0.);
  Spectrum sRecoEFHCNumuSelNC(*loaderFHCPOT, axis, kIsNC && kPIDmu > 0.5 && kQ < 0.);
  Spectrum sRecoEFHCNumuSelNue(*loaderFHCPOT, axis, kIsBeamNue && kPIDmu > 0.5 && kQ < 0.);

  Spectrum sRecoEFHCNueSelTot(*loaderFHCPOT, axis, kPIDe > 0.5);
  Spectrum sRecoEFHCNueSelNumu(*loaderFHCPOT, axis, kIsNumuCC && kPIDe > 0.5);
  Spectrum sRecoEFHCNueSelNC(*loaderFHCPOT, axis, kIsNC && kPIDe > 0.5);
  Spectrum sRecoEFHCNueSelNue(*loaderFHCPOT, axis, kIsBeamNue && kPIDe > 0.5);

  DUNEXSecSyst nushift(nu_MEC_dummy);
  DUNEXSecSyst nubarshift(nubar_MEC_dummy);
  std::map<const ISyst*,double> shiftmap;
  shiftmap.emplace(&nushift, +1.);
  shiftmap.emplace(&nubarshift, +1.);
  SystShifts shifts(shiftmap);
  Spectrum sRecoEFHCNumuSelTotUp(*loaderFHCPOT, axis, kPIDmu > 0.5 && kQ < 0., shifts);
  Spectrum sRecoEFHCNumuSelNumuUp(*loaderFHCPOT, axis, kIsNumuCC && !kIsAntiNu && kPIDmu > 0.5 && kQ < 0., shifts);
  Spectrum sRecoEFHCNumuSelNumubarUp(*loaderFHCPOT, axis, kIsNumuCC && kIsAntiNu && kPIDmu > 0.5 && kQ < 0., shifts);
  Spectrum sRecoEFHCNumuSelNCUp(*loaderFHCPOT, axis, kIsNC && kPIDmu > 0.5 && kQ < 0., shifts);
  Spectrum sRecoEFHCNumuSelNueUp(*loaderFHCPOT, axis, kIsBeamNue && kPIDmu > 0.5 && kQ < 0., shifts);

  PredictionScaleComp pred(*loaderFHCPOT, 
                           axis,
                           kPIDmu > 0.5 && kQ < 0.,
                           GetDUNEXSecSysts());

  loaderFHC.Go();

  new TCanvas;
  sMuPIDFHCNumu.ToTH1(1.47e21, kRed)->Draw("hist");
  sMuPIDFHCNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sMuPIDFHCNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "numu_PID.png" );

  new TCanvas;
  sEPIDFHCNumu.ToTH1(1.47e21, kRed)->Draw("hist");
  sEPIDFHCNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sEPIDFHCNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "nue_PID.png" );

  new TCanvas;
  sRecoEFHCNumuSelTot.ToTH1(1.47e21)->Draw("hist");
  sRecoEFHCNumuSelNumu.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoEFHCNumuSelNumubar.ToTH1(1.47e21, kGreen+2)->Draw("hist same");
  sRecoEFHCNumuSelNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoEFHCNumuSelNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "fhc_numu_selection.png" );

  new TCanvas;
  sRecoEFHCNueSelTot.ToTH1(1.47e21)->Draw("hist");
  sRecoEFHCNueSelNumu.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoEFHCNueSelNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoEFHCNueSelNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "fhc_nue_selection.png" );

  new TCanvas;
  sRecoEFHCNumuSelTotUp.ToTH1(1.47e21)->Draw("hist");
  sRecoEFHCNumuSelNumuUp.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoEFHCNumuSelNumubarUp.ToTH1(1.47e21, kGreen+2)->Draw("hist same");
  sRecoEFHCNumuSelNCUp.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoEFHCNumuSelNueUp.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "fhc_numu_selection_shift.png" );

  new TCanvas;

  osc::NoOscillations noosc;

  Spectrum fake = pred.PredictSyst(&noosc, shiftmap).FakeData(1.47e21);

  fake.ToTH1(1.47e21, kRed)->Draw("hist");
  pred.Predict(&noosc).ToTH1(1.47e21)->Draw("hist same");

  SingleSampleExperiment expt(&pred, fake);

  Fitter fit(&expt, {}, GetDUNEXSecSysts());
  SystShifts seed = SystShifts::Nominal();
  fit.Fit(seed);

  Spectrum bf = pred.PredictSyst(&noosc, shiftmap);
  Spectrum bf2 = pred.PredictSyst(&noosc, seed);
  bf.ToTH1(1.47e21, kBlue)->Draw("hist same");
  bf2.ToTH1(1.47e21, kBlue, 7)->Draw("hist same");

  gPad->Print("fhc_fit.png");


  // RHC
  Spectrum sMuPIDRHCNumu(*loaderRHCPOT, axisPIDmu, kIsNumuCC);
  Spectrum sMuPIDRHCNC(*loaderRHCPOT, axisPIDmu, kIsNC);
  Spectrum sMuPIDRHCNue(*loaderRHCPOT, axisPIDmu, kIsBeamNue);

  Spectrum sEPIDRHCNumu(*loaderRHCPOT, axisPIDe, kIsNumuCC);
  Spectrum sEPIDRHCNC(*loaderRHCPOT, axisPIDe, kIsNC);
  Spectrum sEPIDRHCNue(*loaderRHCPOT, axisPIDe, kIsBeamNue);
  
  Spectrum sRecoERHCNumuSelTot(*loaderRHCPOT, axis, kPIDmu > 0.5 && kQ > 0.);
  Spectrum sRecoERHCNumuSelNumu(*loaderRHCPOT, axis, kIsNumuCC && kPIDmu > 0.5 && kQ > 0. && !kIsAntiNu);
  Spectrum sRecoERHCNumuSelNumubar(*loaderRHCPOT, axis, kIsNumuCC && kPIDmu > 0.5 && kQ > 0. && kIsAntiNu);
  Spectrum sRecoERHCNumuSelNC(*loaderRHCPOT, axis, kIsNC && kPIDmu > 0.5 && kQ > 0.);
  Spectrum sRecoERHCNumuSelNue(*loaderRHCPOT, axis, kIsBeamNue && kPIDmu > 0.5 && kQ > 0.);

  Spectrum sRecoERHCNueSelTot(*loaderRHCPOT, axis, kPIDe > 0.5);
  Spectrum sRecoERHCNueSelNumu(*loaderRHCPOT, axis, kIsNumuCC && kPIDe > 0.5);
  Spectrum sRecoERHCNueSelNC(*loaderRHCPOT, axis, kIsNC && kPIDe > 0.5);
  Spectrum sRecoERHCNueSelNue(*loaderRHCPOT, axis, kIsBeamNue && kPIDe > 0.5);

//  DUNEXSecSyst nubarshift(_nubar_MEC_dummy);
//  SystShifts nubarshifts( &nubarshift, +1);
  Spectrum sRecoERHCNumuSelTotUp(*loaderRHCPOT, axis, kPIDmu > 0.5 && kQ > 0., shifts);
  Spectrum sRecoERHCNumuSelNumuUp(*loaderRHCPOT, axis, kIsNumuCC && kPIDmu > 0.5 && kQ > 0. && !kIsAntiNu, shifts);
  Spectrum sRecoERHCNumuSelNumubarUp(*loaderRHCPOT, axis, kIsNumuCC && kPIDmu > 0.5 && kQ > 0. && kIsAntiNu, shifts);
  Spectrum sRecoERHCNumuSelNCUp(*loaderRHCPOT, axis, kIsNC && kPIDmu > 0.5 && kQ > 0., shifts);
  Spectrum sRecoERHCNumuSelNueUp(*loaderRHCPOT, axis, kIsBeamNue && kPIDmu > 0.5 && kQ > 0., shifts);

  PredictionScaleComp RHCpred(*loaderRHCPOT, 
                            axis,
                            kPIDmu > 0.5 && kQ > 0.,
                            GetDUNEXSecSysts());

  loaderRHC.Go();

  new TCanvas;
  sRecoERHCNumuSelTot.ToTH1(1.47e21)->Draw("hist");
  sRecoERHCNumuSelNumu.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoERHCNumuSelNumubar.ToTH1(1.47e21, kGreen+2)->Draw("hist same");
  sRecoERHCNumuSelNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoERHCNumuSelNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "rhc_numu_selection.png" );

  new TCanvas;
  sRecoERHCNueSelTot.ToTH1(1.47e21)->Draw("hist");
  sRecoERHCNueSelNumu.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoERHCNueSelNC.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoERHCNueSelNue.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "rhc_nue_selection.png" );

  new TCanvas;
  sRecoERHCNumuSelTotUp.ToTH1(1.47e21)->Draw("hist");
  sRecoERHCNumuSelNumuUp.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoERHCNumuSelNumubarUp.ToTH1(1.47e21, kGreen+2)->Draw("hist same");
  sRecoERHCNumuSelNCUp.ToTH1(1.47e21, kBlue)->Draw("hist same");
  sRecoERHCNumuSelNueUp.ToTH1(1.47e21, kMagenta)->Draw("hist same");
  gPad->Print( "rhc_numu_selection_shift.png" );

  new TCanvas;

  //osc::NoOscillations noosc;

  Spectrum fakeRHC = RHCpred.PredictSyst(&noosc, shiftmap).FakeData(1.47e21);

  fakeRHC.ToTH1(1.47e21, kRed)->Draw("hist");
  RHCpred.Predict(&noosc).ToTH1(1.47e21)->Draw("hist same");

  SingleSampleExperiment exptRHC(&RHCpred, fakeRHC);

  Fitter RHCfit(&exptRHC, {}, GetDUNEXSecSysts());
  SystShifts RHCseed = SystShifts::Nominal();
  RHCfit.Fit(RHCseed);

  Spectrum RHCbf = RHCpred.PredictSyst(&noosc, shiftmap);
  Spectrum RHCbf2 = RHCpred.PredictSyst(&noosc, RHCseed);
  RHCbf.ToTH1(1.47e21, kBlue)->Draw("hist same");
  RHCbf2.ToTH1(1.47e21, kBlue, 7)->Draw("hist same");

  gPad->Print("rhc_fit.png");

}
