#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH2D.h>
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "ExampleSelection.h"
#include "ExampleTools.h"
#include "sbncode/sbnanalysis/util/Interaction.hh"
#include "sbncode/sbnanalysis/core/Processor/ProcessorBase.hh"
#include "sbncode/sbnanalysis/core/Processor/ProcessorCreator.hh"
#include "sbncode/sbnanalysis/core/DataTypes/Event.hh"

namespace sbnanalysis {

namespace ana {
  namespace ExampleAnalysis {

ExampleSelection::ExampleSelection()
    : ProcessorBase(), fEventCounter(0), fNuCount(0) {}


void ExampleSelection::Initialize(fhicl::ParameterSet* config) {
  // Make a histogram
  art::ServiceHandle<art::TFileService> tfs;
  fNuVertexXZHist = tfs->make<TH2D>("nu_vtx_XZ", "",
                                    100, -1000, 1000, 100, -1000, 1000);

  // Load configuration parameters
  fMyParam = 0;
  fTruthTag = { "generator" };

  if (config) {
    fMyParam = config->get<int>("parameter", 0);

    fTruthTag = \
      { config->get<std::string>("MCTruthTag", "generator") };
  }

  // Add custom branches
  AddBranch("nucount", &fNuCount);
  AddBranch("myvar", &fMyVar);

  // Use some library code
  hello();
}


void ExampleSelection::Finalize() {
  // Output our histograms to the ROOT file
  //fOutputFile->cd();
  fNuVertexXZHist->Write();
}


bool ExampleSelection::ProcessEvent(
    const art::Event& ev,
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

    // Add in the "reconstructed" interaction
    // Contruct truth information from the provided vector
    Event::RecoInteraction interaction(truth[i], i);
    interaction.reco_energy = \
      util::ECCQE(mctruth.GetNeutrino().Nu().Momentum().Vect(), mctruth.GetNeutrino().Lepton().E());
    reco.push_back(interaction);
  }

  return true;
}

REGISTER_PROCESSOR(ExampleSelection)

  }  // namespace ExampleAnalysis
}  // namespace ana

}

