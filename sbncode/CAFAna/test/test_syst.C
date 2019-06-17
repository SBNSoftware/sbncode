#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
using namespace ana;

#include "Utilities/rootlogon.C"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"

void test_syst()
{
  rootlogon(); // style

  SpectrumLoader loaderNumu("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.1/numutest.root");

  auto* loaderNumuBeam  = loaderNumu.LoaderForRunPOT(20000001);

  const Var kRecoE = SIMPLEVAR(dune.Ev_reco);
  const Var kPID = SIMPLEVAR(dune.mvaresult);

  const HistAxis axis("Reconstructed energy (GeV)",
                      Binning::Simple(40, 0, 10),
                      kRecoE);

  const HistAxis axisPID("mvaresult",
                         Binning::Simple(30, -1.5, +1.5),
                         kPID);

  SystShifts up(&kFirstFluxSyst, +1);
  SystShifts dn(&kFirstFluxSyst, -1);

  Spectrum sPIDFHCTot(*loaderNumuBeam, axisPID, kNoCut);
  Spectrum sPIDFHCTotUp(*loaderNumuBeam, axisPID, kNoCut, up);
  Spectrum sPIDFHCTotDn(*loaderNumuBeam, axisPID, kNoCut, dn);

  Spectrum sRecoEFHCNumuSelTot(*loaderNumuBeam, axis, kPID > 0.5);
  Spectrum sRecoEFHCNumuSelTotUp(*loaderNumuBeam, axis, kPID > 0.5, up);
  Spectrum sRecoEFHCNumuSelTotDn(*loaderNumuBeam, axis, kPID > 0.5, dn);

  loaderNumu.Go();

  new TCanvas;
  sPIDFHCTot.ToTH1(1.47e21)->Draw("hist");
  sPIDFHCTotUp.ToTH1(1.47e21, kRed)->Draw("hist same");
  sPIDFHCTotDn.ToTH1(1.47e21, kBlue)->Draw("hist same");

  new TCanvas;
  sRecoEFHCNumuSelTot.ToTH1(1.47e21)->Draw("hist");
  sRecoEFHCNumuSelTotUp.ToTH1(1.47e21, kRed)->Draw("hist same");
  sRecoEFHCNumuSelTotDn.ToTH1(1.47e21, kBlue)->Draw("hist same");
}
