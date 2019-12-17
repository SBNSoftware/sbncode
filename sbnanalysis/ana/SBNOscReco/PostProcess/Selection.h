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

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

#include "Histograms.h"
#include "Cuts.h"
#include "Normalize.h"
#include "ROC.h"

#include "../Histograms/TrajHistograms.h"
#include "../Histograms/CRTHisto.h"
#include "../Histograms/DynamicSelector.h"
#include "../Data/RecoEvent.h"
#include "../Data/MCType.h" 

class TTree;

namespace ana {
  namespace SBNOsc {

class Selection: public core::PostProcessorBase {
public:
  // implementing PostProcessor
  void FileCleanup(TTree *eventTree) {}
  void FileSetup(TFile *f, TTree *eventTree);
  void Initialize(fhicl::ParameterSet *config);
  void ProcessEvent(const event::Event *event);
  void ProcessSubRun(const SubRun *subrun);
  void ProcessFileMeta(const FileMeta *meta);
  void Finalize();

private:
  sbnd::CRTGeoAlg *fCRTGeo;
  Cuts fCuts;
  ROC fROC;
  Histograms fHistograms;
  Histograms fCosmicHistograms;
  Histograms *fHistsToFill;
  TrajHistograms fTrajHistograms;

  CRTHistos fCRTCosmicHistos;
  CRTHistos fCRTNeutrinoHistos;

  Normalize fNormalize;
  fhicl::ParameterSet fCutConfig;
  double fCRTHitDistance;
  bool fDoNormalize;
  bool fFillAllTracks;
  double fGoalPOT;
  double fNCosmicData;

  std::vector<numu::TrackSelector> fTrackSelectors;
  std::vector<std::string> fTrackSelectorNames;
  std::vector<std::string> fTrajHistoNames;

  std::vector<numu::ROOTValue> fTrackProfileValues;

  numu::RecoEvent *fRecoEvent;
  TFile *fOutputFile;

  int fFileIndex;
  numu::MCType fFileType;

};

  }   // namespace SBNOsc
}   // namespace ana

#endif// __sbnanalysis_ana_SBNOsc_Selection__
