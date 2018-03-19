#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "SBNProposalSelection.h"
#include "SBNProposalTools.h"

namespace ana {
  namespace SBNProposalAnalysis {

SBNProposalSelection::SBNProposalSelection() : SelectionBase(), fNuCount(0), fEventCounter(0) {}


void SBNProposalSelection::Initialize(Json::Value* config) {
  // Make a histogram
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", "",
                             100, -1000, 1000, 100, -1000, 1000);

  // Load configuration parameters
  fMyParam = 0;
  fTruthTag = { "generator" };

  if (config) {
    fMyParam = (*config)["SBNProposalAnalysis"].get("parameter", 0).asInt();
    fTruthTag = { (*config)["SBNProposalAnalysis"].get("MCTruthTag", "generator").asString() };
  }

  // Add custom branches
  AddBranch("nucount", &fNuCount);
  AddBranch("myvar", &fMyVar);

  // Use some library code
  hello();
}


void SBNProposalSelection::Finalize() {
  // Output our histograms to the ROOT file
  fOutputFile->cd();
  fNuVertexXZHist->Write();
}


bool SBNProposalSelection::ProcessEvent(gallery::Event& ev) {
  if (fEventCounter % 10 == 0) {
    std::cout << "SBNProposalSelection: Processing event " << fEventCounter << std::endl;
  }
  fEventCounter++;

  // Grab a data product from the event
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(fTruthTag);

  // Fill in the custom branches
  fNuCount = mctruths.size();  // Number of neutrinos in this event
  fMyVar = fMyParam;
  
  // Iterate through the neutrinos
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);

    // Fill neutrino vertex position histogram
    fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(),
                          mctruth.GetNeutrino().Nu().Vz());
  }

  return true;
}

  }  // namespace SBNProposalAnalysis
}  // namespace ana


// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::SBNProposalAnalysis::SBNProposalSelection)

