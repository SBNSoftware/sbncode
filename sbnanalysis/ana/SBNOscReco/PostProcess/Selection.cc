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

#include "../NumuReco/PrimaryTrack.h"
#include "../NumuReco/TruthMatch.h"

namespace ana {
namespace SBNOsc {
  void Selection::Initialize(fhicl::ParameterSet *config) {
    fOutputFile = new TFile(config->get<std::string>("OutputFile", "output.root").c_str(), "CREATE");
    fOutputFile->cd();
    fCRTHitDistance = config->get<double>("CRTHitDistance", -1);

    fFileIndex = 0;

    fCutConfig = config->get<fhicl::ParameterSet>("Cuts");

    fDoNormalize = config->get<bool>("DoNormalize", false);
    fGoalPOT = config->get<double>("GoalPOT", 0.);
    fFillAllTracks = config->get<bool>("FillAllTracks", true);
    if (fDoNormalize) {
      fNormalize.Initialize(config->get<fhicl::ParameterSet>("Normalize", {}));
      fFileTypes = config->get<std::vector<std::string>>("FileTypes");
      fNCosmicData = config->get<double>("NCosmicData", 0.);
    }
    fTrajHistograms.Initialize();
    fROC.Initialize();
    fRecoEvent = NULL;
  }

  void Selection::FileSetup(TFile *f, TTree *eventTree) {
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
    
    TrajHistograms from_file;
    from_file.Get(*f);
    fTrajHistograms.Add(from_file);
  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    std::cout << "Event: " << core_event->metadata.eventID << std::endl;

    // update each reco Interaction to a smarter primary track selector
    unsigned i = 0;
    while (i < fRecoEvent->reco.size()) {
      int primary_track = numu::SelectLongestIDdMuon(fRecoEvent->reco_tracks, fRecoEvent->reco[i].slice);

      // remove vertices without a good primary track
      if (primary_track < 0) {
        fRecoEvent->reco.erase(fRecoEvent->reco.begin() + i); 
        continue;
      }

      fRecoEvent->reco[i].slice.primary_track_index = primary_track;
      fRecoEvent->reco[i].primary_track = fRecoEvent->reco_tracks.at(primary_track);
      // re-do truth matching
      fRecoEvent->reco[i].match = numu::InteractionTruthMatch(fRecoEvent->truth, fRecoEvent->reco_tracks, fRecoEvent->reco[i]);

      i += 1;
    }

    // make sure no two vertices match to the same true neutrino interaction
    numu::CorrectMultiMatches(*fRecoEvent, fRecoEvent->reco);

    fHistsToFill->Fill(*fRecoEvent, *core_event, fCuts, fFillAllTracks);
    if (fDoNormalize) {
      if (fFileType == "cosmic") fNormalize.AddCosmicEvent(*core_event);
      else fNormalize.AddNeutrinoEvent(*core_event);
    }
    fROC.Fill(fCuts, *fRecoEvent);
  }

  void Selection::Finalize() {
    fOutputFile->cd();
    if (fDoNormalize) {
      fNeutrinoHistograms.Scale(fNormalize.ScaleNeutrino(fGoalPOT));
      fCosmicHistograms.Scale(fNormalize.ScaleCosmic(fGoalPOT));
      fHistograms.Add(fNeutrinoHistograms);
      fHistograms.Add(fCosmicHistograms);
      fROC.Normalize(fNormalize.ScaleNeutrino(fGoalPOT), fNormalize.ScaleCosmic(fGoalPOT));  
    }
    fROC.BestCuts();
    fROC.Write();
    fHistograms.Write();
    fTrajHistograms.Write();
  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

