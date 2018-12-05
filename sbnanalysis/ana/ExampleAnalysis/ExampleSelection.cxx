#include <iostream>
#include <vector>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ExampleSelection.h"
#include "ExampleTools.h"
#include "core/Event.hh"
#include "core/ServiceManager.hh"
#include "larcorealg/Geometry/GeometryCore.h"

namespace ana {
  namespace ExampleAnalysis {

ExampleSelection::ExampleSelection()
    : SelectionBase(), fNuCount(0), fEventCounter(0) {}


void ExampleSelection::Initialize(fhicl::ParameterSet* config) {
  // Make a histogram
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", "",
                             100, -1000, 1000, 100, -1000, 1000);

  // Load configuration parameters
  fMyParam = 0;
  fTruthTag = { "generator" };

  if (config) {
    fhicl::ParameterSet pconfig = \
      config->get<fhicl::ParameterSet>("ExampleAnalysis");

    fMyParam = pconfig.get<int>("parameter", 0);

    fTruthTag = \
      { pconfig.get<std::string>("MCTruthTag", "generator") };
  }

  // Add custom branches
  AddBranch("nucount", &fNuCount);
  AddBranch("myvar", &fMyVar);

  // Get access to services
  fServiceManager = new core::ServiceManager(core::kICARUS);
  std::cout << "Detector: " << fServiceManager->GetGeometryService()->DetectorName() << std::endl;

  core::ServiceManager* fServiceManager2 = new core::ServiceManager(core::kSBND);
  std::cout << "Detector: " << fServiceManager2->GetGeometryService()->DetectorName() << std::endl;

  core::ServiceManager* fServiceManager3 = new core::ServiceManager(core::kUBOONE);
  std::cout << "Detector: " << fServiceManager3->GetGeometryService()->DetectorName() << std::endl;

  // Use some library code
  hello();
}


void ExampleSelection::Finalize() {
  // Output our histograms to the ROOT file
  fOutputFile->cd();
  fNuVertexXZHist->Write();
}


bool ExampleSelection::ProcessEvent(
    const gallery::Event& ev,
    const std::vector<Event::Interaction>& truth,
    std::vector<Event::RecoInteraction>& reco) {
  if (fEventCounter % 10 == 0) {
    std::cout << "ExampleSelection: Processing event "
              << fEventCounter << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);

  // Fill in the custom branches
  fNuCount = mctruths.size();  // Number of neutrinos in this event
  fMyVar = fMyParam;

  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);

    // Fill neutrino vertex position histogram
    if (mctruth.NeutrinoSet()) {
      fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(),
                            mctruth.GetNeutrino().Nu().Vz());
    }
  }

  return true;
}

  }  // namespace ExampleAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::ExampleAnalysis::ExampleSelection)

