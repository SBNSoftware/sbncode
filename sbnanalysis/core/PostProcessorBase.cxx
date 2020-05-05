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
#include "larcorealg/Geometry/AuxDetGeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

namespace core {

PostProcessorBase::PostProcessorBase(): fProviderManager(NULL), fConfigExperimentID(-1), fNWorkers(1) {}


PostProcessorBase::~PostProcessorBase() {}


void PostProcessorBase::Initialize(char* config, const std::string &output_fname, unsigned n_workers) {
  fhicl::ParameterSet* cfg = LoadConfig(config);
  if (cfg == NULL) cfg = new fhicl::ParameterSet;
  if (output_fname.size() != 0) cfg->put("OutputFile", output_fname);
  fConfigExperimentID = cfg->get("ExperimentID", -1);

  if (fConfigExperimentID >= 0) {
    fProviderManager = new ProviderManager((Experiment)fConfigExperimentID, "", false); 
  }

  fNWorkers = n_workers;

  Initialize(cfg);
}

void PostProcessorBase::ProcessFile(const std::string &fname) {
  // get ROOT file
  TFile f(fname.c_str());
  if (f.IsZombie()) {
    std::cerr << "Failed openning file: " << fname << ". "
              << "Cleaning up and exiting." << std::endl;
    return;
  }

  TTree *event_tree = 0;
  event::Event *event = 0;

  TTree *subrun_tree = 0;
  SubRun *subrun = 0;
  
  TTree *filemeta_tree = 0;
  FileMeta *filemeta = 0;

  f.GetObject("sbnana", event_tree);
  event_tree->SetBranchAddress("events", &event);
  FileSetup(&f, event_tree);
  // process all events
  for (int event_ind = 0; event_ind < event_tree->GetEntries(); event_ind++) {
    event_tree->GetEntry(event_ind);
    ProcessEvent(event);
  }
  // process all subruns
  f.GetObject("sbnsubrun", subrun_tree);
  if (subrun_tree == NULL) {
    std::cerr << "Error: NULL subrun tree" << std::endl;
  }
  subrun_tree->SetBranchAddress("subruns", &subrun);

  for (int subrun_ind = 0; subrun_ind < subrun_tree->GetEntries(); subrun_ind++) {
    subrun_tree->GetEntry(subrun_ind);
    ProcessSubRun(subrun);
  }

  // process all the file meta-data
  f.GetObject("sbnfilemeta", filemeta_tree);
  if (filemeta_tree == NULL) {
    std::cerr << "Error: NULL filemeta tree" << std::endl;
  }
  filemeta_tree->SetBranchAddress("filemeta", &filemeta);
  for (int filemeta_ind = 0; filemeta_ind < filemeta_tree->GetEntries(); filemeta_ind++) {
    filemeta_tree->GetEntry(filemeta_ind);
    ProcessFileMeta(filemeta);
  } 

  FileCleanup(event_tree);
}

unsigned PostProcessorBase::WorkerID() {
  std::thread::id this_id = std::this_thread::get_id();
  for (unsigned index = 0; index < fThreadIDs.size(); index++) {
    if (this_id == fThreadIDs[index]) return index;
  }
  assert(false);
  return 0;
}

void PostProcessorBase::Run(std::vector<std::string> inputFiles) {
  // single threaded
  if (fNWorkers == 1) {
    fThreadIDs.push_back(std::this_thread::get_id());
    for (auto const& fname: inputFiles) {
      ProcessFile(fname);
    } 
  }
  // multi-threaded
  else {
    fThreadIDs = std::vector<std::thread::id>(fNWorkers);
    std::atomic<unsigned> f_index = 0;
    std::atomic<unsigned> w_index = 0;
    std::vector<std::thread> workers;
    for (unsigned i = 0; i < fNWorkers; i++) {
      workers.emplace_back([&]() {
        // setup the mapping to thread indices
        unsigned this_index = w_index.fetch_add(1);
        fThreadIDs[this_index] = std::this_thread::get_id();
        InitializeThread();
        // start fetching files
        unsigned file_index = f_index.fetch_add(1);
        while (file_index < inputFiles.size()) {
          ProcessFile(inputFiles[file_index]);
          file_index = f_index.fetch_add(1);
        } 
      });
    }
    for (std::thread &w: workers) w.join();
  }

  if (fProviderManager != NULL) {
    delete fProviderManager;
    fProviderManager = NULL;
  }

  // teardown
  Finalize();
}

}  // namespace core

