#ifndef __sbnanalysis_ana_SBNOsc_PIDDataMaker__
#define __sbnanalysis_ana_SBNOsc_PIDDataMaker__

/**
 * \file PIDDataMaker.h
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
#include "../Data/PIDData.h"
#include "../Data/RecoEvent.h"

namespace ana {
  namespace SBNOsc {

class PIDDataMaker: public core::PostProcessorBase {
public:
  // implementing PostProcessor
  //void FileCleanup(TTree *eventTree) {}
  void Initialize(fhicl::ParameterSet *config);
  void ProcessEvent(const event::Event *event);
  void Finalize();

  void FileSetup(TFile *f, TTree *eventTree);

private:
  Cuts fCuts;
  TFile *fOutputFile;
  TTree *fPIDTree;
  numu::flat::PIDData *fPIDData;
  numu::RecoEvent *fRecoEvent;
};

  }   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_PIDDataMaker__
