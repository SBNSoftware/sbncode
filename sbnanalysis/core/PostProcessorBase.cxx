#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include "fhiclcpp/ParameterSet.h"
#include "Event.hh"
#include "Loader.hh"
#include "PostProcessorBase.hh"

namespace core {

PostProcessorBase::PostProcessorBase(): fEvent(NULL) {}


PostProcessorBase::~PostProcessorBase() {}


void PostProcessorBase::Initialize(char* config) {
  fhicl::ParameterSet* cfg = LoadConfig(config);
  Initialize(cfg);
}


void PostProcessorBase::Run(std::vector<std::string> inputFiles) {
  for (auto const& fname: inputFiles) {
    // get ROOT file
    TFile f(fname.c_str());
    if (f.IsZombie()) {
      std::cerr << "Failed openning file: " << fname << ". "
                << "Cleaning up and exiting." << std::endl;
      break;
    }

    // set Event
    fEventTree = (TTree *) f.Get("sbnana");
    fEventTree->SetBranchAddress("events", &fEvent);
    FileSetup(fEventTree);

    // process all events
    for (int event_ind = 0; event_ind < fEventTree->GetEntries(); event_ind++) {
      fEventTree->GetEntry(event_ind);
      ProcessEvent(fEvent);
    }

    FileCleanup(fEventTree);
  } 

  // teardown
  Finalize();
}

}  // namespace core

