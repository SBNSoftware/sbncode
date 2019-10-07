#include <string>
#include <vector>
#include <TChain.h>
#include <TTree.h>

#include <map>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>

#include <TFile.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TCanvas.h>
#include <core/Event.hh>

#include "Selection.h"

namespace ana {
namespace SBNOsc {
  void Selection::Initialize(fhicl::ParameterSet *config) {
    fOutputFile = new TFile("output.root", "CREATE");
    fOutputFile->cd();
    fCRTHitDistance = config->get<double>("CRTHitDistance", -1);

    fFileIndex = 0;

    fCutConfig = config->get<fhicl::ParameterSet>("Cuts");

    fDoNormalize = config->get<bool>("DoNormalize", false);
    fGoalPOT = config->get<double>("GoalPOT", 0.);
    if (fDoNormalize) {
      fNormalize.Initialize(config->get<fhicl::ParameterSet>("Normalize", {}));
      fFileTypes = config->get<std::vector<std::string>>("FileTypes");
    }
    fRecoEvent = NULL;
  }

  void Selection::FileSetup(TTree *eventTree) {
    eventTree->SetBranchAddress("reco_event", &fRecoEvent);
    fCuts.Initialize(fCutConfig, fProviderManager->GetGeometryProvider());
    if (fDoNormalize) {
      fFileType = fFileTypes[fFileIndex];
      if (fFileType == "cosmic") {
        fHistsToFill = &fCosmicHistograms;
      }
      else if (fFileType == "neutrino") {
        fHistsToFill = &fNeutrinoHistograms;
      }
      else assert(false);
    }
    else {
      fHistsToFill = &fHistograms;
    }
    fFileIndex += 1;

  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    std::cout << "Event: " << core_event->metadata.eventID << std::endl;
    fHistsToFill->Fill(*fRecoEvent, *core_event, fCuts);
    if (fDoNormalize) {
      if (fFileType == "cosmic") fNormalize.AddCosmicEvent(*core_event);
      else fNormalize.AddNeutrinoEvent(*core_event);
    }
  }

  void Selection::Finalize() {
    fOutputFile->cd();
    if (fDoNormalize) {
      fNeutrinoHistograms.Scale(fNormalize.ScaleNeutrino(fGoalPOT));
      fCosmicHistograms.Scale(fNormalize.ScaleCosmic(fGoalPOT));
      fHistograms.Add(fNeutrinoHistograms);
      fHistograms.Add(fCosmicHistograms);
    }
    fHistograms.Write();
  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

