#ifndef __sbnanalysis_ana_SBNOsc_Selection__
#define __sbnanalysis_ana_SBNOsc_Selection__

/**
 * \file Selection.h
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
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMatrixDSym.h>

#include "Histograms.h"
#include "Cuts.h"
#include "../Data/RecoEvent.h"

class TTree;

namespace ana {
  namespace SBNOsc {

class Selection: public core::PostProcessorBase {
public:
  // Constructor
  Selection() {}
  
  // implementing PostProcessor
  void FileCleanup(TTree *eventTree) {}
  void FileSetup(TTree *eventTree);
  void Initialize(fhicl::ParameterSet *config);
  void ProcessEvent(const event::Event *event);
  void Finalize();

private:
  Cuts fCuts;
  Histograms fHistograms;
  double fCRTHitDistance;
  numu::RecoEvent *fRecoEvent;
  TFile *fOutputFile;

};

  }   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Selection__
