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

#include "../Histograms/Histograms.h"
#include "../Histograms/TrajHistograms.h"
#include "../Histograms/Cuts.h"
#include "Normalize.h"
#include "ROC.h"
#include "../Data/RecoEvent.h"

class TTree;

namespace ana {
  namespace SBNOsc {

class Selection: public core::PostProcessorBase {
public:
  // Constructor
  Selection():
    fHistograms(),
    fNeutrinoHistograms("Neutrino"),
    fCosmicHistograms("Cosmic") {}
  
  // implementing PostProcessor
  void FileCleanup(TTree *eventTree) {}
  void FileSetup(TFile *f, TTree *eventTree);
  void Initialize(fhicl::ParameterSet *config);
  void ProcessEvent(const event::Event *event);
  void Finalize();

private:
  Cuts fCuts;
  ROC fROC;
  Histograms fHistograms;
  Histograms fNeutrinoHistograms;
  Histograms fCosmicHistograms;
  Histograms *fHistsToFill;
  TrajHistograms fTrajHistograms;
  Normalize fNormalize;
  fhicl::ParameterSet fCutConfig;
  double fCRTHitDistance;
  bool fDoNormalize;
  bool fFillAllTracks;
  double fGoalPOT;
  double fNCosmicData;

  numu::RecoEvent *fRecoEvent;
  TFile *fOutputFile;

  std::vector<std::string> fFileTypes;
  int fFileIndex;
  std::string fFileType;

};

  }   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Selection__
