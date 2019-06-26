#include <iostream>
#include <vector>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"
#include "NueSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

NueSelection::NueSelection() : SelectionBase(), fEventCounter(0), fNuCount(0) {}


void NueSelection::Initialize(fhicl::ParameterSet* config) {
  // Load configuration parameters
  fTruthTag = { "generator" };

  fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("SBNOsc");

  if (config) {
    fTruthTag = { pconfig.get<std::string>("MCTruthTag", "generator") };
  }

  hello();
}


void NueSelection::Finalize() {}


bool NueSelection::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
  if (fEventCounter % 10 == 0) {
    std::cout << "NueSelection: Processing event " << fEventCounter << " "
              << "(" << fNuCount << " neutrinos selected)"
              << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    if (nu.CCNC() == simb::kCC && nu.Mode() == 0 && nu.Nu().PdgCode() == 12) {
      event::RecoInteraction interaction(i);
      reco.push_back(interaction);
    }
  }

  bool selected = !reco.empty();

  if (selected) {
    fNuCount++;
  }

  return selected;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NueSelection)

