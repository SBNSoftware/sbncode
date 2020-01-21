// Introduction to systematics
// cafe demo3.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Analysis/ExpInfo.h"

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"

void demo3()
{
  const std::string fname = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";
  SpectrumLoader loader(fname);
  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);
  const HistAxis axEnergy("True energy (GeV)", Binning::Simple(50, 0, 5), kRecoE);
  Spectrum sNominal(loader, axEnergy, kNoCut, kNoShift, kWeight);

  // Systematics are implemented as objects that either compute a weight for
  // each event...
  class TailScaleSyst: public ISyst
  {
  public:
    TailScaleSyst() : ISyst("tailScale", "High energy scale factor") {}

    // 50% 1sigma shift on normalization of events above 2GeV
    void Shift(double sigma, caf::SRProxy* sr, double& weight) const
    {
      const double E = sr->reco.reco_energy;
      /**/ if(E > 2) weight *= 1+.5*sigma;
      else if(E > 1) weight *= 1+.5*sigma*(E-1);
    }
  };
  const TailScaleSyst kTailScaleSyst;


  // ...or rewrite the event record itself
  class EnergyScaleSyst: public ISyst
  {
  public:
    EnergyScaleSyst() : ISyst("energyScale", "Energy scale") {}

    // 20% 1sigma shift on all event energies
    void Shift(double sigma, caf::SRProxy* sr, double& weight) const
    {
      sr->reco.reco_energy *= 1+.2*sigma;
    }
  };
  const EnergyScaleSyst kEnergyScaleSyst;

  // A SystShifts packages up a bunch of ISyst objects with the values they
  // should be shifted by.
  SystShifts tailUp(&kTailScaleSyst, +1);
  SystShifts tailDn(&kTailScaleSyst, -1);

  // Specify SystShifts in the Spectrum constructor and it will be filled with
  // the shifted variables.
  Spectrum sTailUp(loader, axEnergy, kNoCut, tailUp, kWeight);
  Spectrum sTailDn(loader, axEnergy, kNoCut, tailDn, kWeight);


  SystShifts energyUp(&kEnergyScaleSyst, +1);
  SystShifts energyDn(&kEnergyScaleSyst, -1);

  Spectrum sEnergyUp(loader, axEnergy, kNoCut, energyUp, kWeight);
  Spectrum sEnergyDn(loader, axEnergy, kNoCut, energyDn, kWeight);

  // Fill the spectrum and all the systematic variants
  loader.Go();

  sTailUp.ToTH1(kPOTnominal, kRed)->Draw("hist");
  sTailDn.ToTH1(kPOTnominal, kBlue)->Draw("hist same");
  sNominal.ToTH1(kPOTnominal)->Draw("hist same");

  new TCanvas;

  sEnergyDn.ToTH1(kPOTnominal, kBlue)->Draw("hist");
  sEnergyUp.ToTH1(kPOTnominal, kRed)->Draw("hist same");
  sNominal.ToTH1(kPOTnominal)->Draw("hist same");
}
