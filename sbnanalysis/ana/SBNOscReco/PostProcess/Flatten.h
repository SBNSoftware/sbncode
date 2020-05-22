#ifndef __sbnanalysis_ana_SBNOsc_Flatten__
#define __sbnanalysis_ana_SBNOsc_Flatten__

/**
 * \file Flatten.h
 */

#include "fhiclcpp/ParameterSet.h"
#include "core/PostProcessorBase.hh"

#include <string>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cassert>

#include <TFile.h>
#include <TNtuple.h>

#include "Cuts.h"

#include "core/Event.hh"
#include "../Data/FlatInteraction.h" 
#include "../Data/RecoEvent.h" 

namespace ana {
  namespace SBNOsc {

class Flatten: public core::PostProcessorBase {
public:
  // implementing PostProcessor
  //void FileCleanup(TTree *eventTree) {}
  void FileSetup(TFile *f, TTree *eventTree);
  void Initialize(fhicl::ParameterSet *config);
  void ProcessEvent(const event::Event *event);
  void InitializeThread();
  void ProcessSubRun(const SubRun *subrun);
  void ProcessFileMeta(const FileMeta *meta);
  void Finalize();

private:
  Cuts fCuts;
  std::vector<TFile *> fOutputFiles;
  std::vector<TNtuple *> fNtuples;
  std::vector<numu::flat::FlatInteraction> fInteractions;
  std::vector<numu::MCType> fMCTypes;
  std::vector<numu::RecoEvent *> fRecoEvents;
};

  }   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Flatten__
