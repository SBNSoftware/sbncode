#include <iostream>
#include <vector>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "MySelection.h"
#include "ExampleTools.h"
#include "core/Event.hh"
#include "util/Interaction.hh"
#include "core/ProviderManager.hh"

namespace ana {
  namespace ExampleAnalysis {

MySelection::MySelection()
    : SelectionBase(), fNuCount(0), fEventCounter(0) {}


void MySelection::Initialize(fhicl::ParameterSet* config) {
  // Make a histogram

  // Load configuration parameters
  fTruthTag = { "generator" };

}


void MySelection::Finalize() {
  // Output our histograms to the ROOT file
  fOutputFile->cd();
}


bool MySelection::ProcessEvent(
    const gallery::Event& ev,
    const std::vector<event::Interaction>& truth,
    std::vector<event::RecoInteraction>& reco) {

  if (fEventCounter % 10 == 0) {
    std::cout << "MySelection: Processing event "
              << fEventCounter << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);

  // Fill in the custom branches
  fNuCount = mctruths.size();  // Number of neutrinos in this event

  return true;
}
    
  }  // namespace ExampleAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::ExampleAnalysis::MySelection)

