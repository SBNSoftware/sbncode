#include <iostream>
#include <vector>
#include <TH2D.h>
#include <json/json.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "ExampleSelection.h"
#include "ExampleTools.h"

namespace ana {
  namespace ExampleAnalysis {

ExampleSelection::ExampleSelection() : SelectionBase(), fNuCount(0) {}


void ExampleSelection::Initialize(Json::Value* config) {
  // Here you name the thing that "produced" the data product that you want to
  // look at. In our event dump we see two things:
  //
  // GenieGen.. | corsika..... | .... | std::vector<simb::MCTruth>.... | ....1
  // GenieGen.. | generator... | .... | std::vector<simb::MCTruth>.... | ....1
  //
  // This means that if you want to look at the simb::MCTruth data product
  // there are two "producers." We want to look at neutrinos so we choose
  // "generator" if you wanted to look at cosmics you could pick "coriska."
  mctruths_tag = { "generator" };

  // Make a histogram
  fNuVertexXZHist = new TH2D("nu_vtx_XZ", "",
                             100, -1000, 1000, 100, -1000, 1000);

  // Load configuration parameters
  fMyParam = (*config)["ExampleAnalysis"].get("parameter", 0).asInt();
  
  // Add custom branches
  AddBranch("nucount", &fNuCount);
  AddBranch("myvar", &fMyVar);

  // Use some library code
  hello();
}


void ExampleSelection::Finalize() {
  fOutputFile->cd();
  fNuVertexXZHist->Write();
}


void ExampleSelection::ProcessEvent(gallery::Event& ev) {
  // Grab the data product that you want from the event
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth>>(mctruths_tag);

  // Fill in the custom branch (e.g. number of neutrino interactions)
  fNuCount = mctruths.size();
  fMyVar = fMyParam;
  
  // Now we'll iterate through these 
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);

    // Now for each simb::MCTruth we look at two things      
    fNuVertexXZHist->Fill(mctruth.GetNeutrino().Nu().Vx(),
                          mctruth.GetNeutrino().Nu().Vz());
  }
}


  }  // namespace ExampleAnalysis
}  // namespace ana

DECLARE_SBN_PROCESSOR(ana::ExampleAnalysis::ExampleSelection)

