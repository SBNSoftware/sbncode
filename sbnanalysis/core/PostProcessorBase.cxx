#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include "fhiclcpp/ParameterSet.h"
#include "Event.hh"
#include "Loader.hh"
#include "PostProcessorBase.hh"
#include "Experiment.hh"

#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

namespace core {

PostProcessorBase::PostProcessorBase(): fEvent(NULL), fProviderManager(NULL) {}


PostProcessorBase::~PostProcessorBase() {}


void PostProcessorBase::Initialize(char* config, const std::string &output_fname) {
  fhicl::ParameterSet* cfg = LoadConfig(config);
  if (output_fname.size() != 0) cfg->put("OutputFile", output_fname);
  Initialize(cfg);
  fSubRun = 0;
  fEvent = 0;
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
    // get the Experiment ID
    f.GetObject("experiment", fExperimentID);
    if ((Experiment)fExperimentID->GetVal() != kExpOther) {
      fProviderManager = new ProviderManager((Experiment)fExperimentID->GetVal(), "", false); 
    }

    // set Event
    f.GetObject("sbnana", fEventTree);
    fEventTree->SetBranchAddress("events", &fEvent);

    FileSetup(fEventTree);
    // process all events
    for (int event_ind = 0; event_ind < fEventTree->GetEntries(); event_ind++) {
      fEventTree->GetEntry(event_ind);
      ProcessEvent(fEvent);
    }

    // process all subruns
    f.GetObject("sbnsubrun", fSubRunTree);
    if (fSubRunTree == NULL) {
      std::cerr << "Error: NULL subrun tree" << std::endl;
    }
    fSubRunTree->SetBranchAddress("subruns", &fSubRun);

    for (int subrun_ind = 0; subrun_ind < fSubRunTree->GetEntries(); subrun_ind++) {
      fSubRunTree->GetEntry(subrun_ind);
      ProcessSubRun(fSubRun);
    }

    FileCleanup(fEventTree);

    if (fProviderManager != NULL) {
      delete fProviderManager;
      fProviderManager = NULL;
    }

  } 

  // teardown
  Finalize();
}

}  // namespace core

