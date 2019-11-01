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
#include "../Histograms/DynamicSelector.h"

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

    std::vector<std::vector<std::string>> track_selector_strings = config->get<std::vector<std::vector<std::string>>>("TrackSelectors", {{""}});
    std::vector<std::vector<std::string>> track_selector_names = config->get<std::vector<std::vector<std::string>>>("TrackSelectorNames", {{""}});

    assert(track_selector_strings.size() == track_selector_names.size());

    fTrackSelectorNames = numu::MultiplyNames(track_selector_names);
    fTrackSelectors = numu::MultiplySelectors(track_selector_strings);

    assert(fTrackSelectorNames.size() == fTrackSelectors.size());

    fTrajHistoNames = numu::MultiplyNames(config->get<std::vector<std::vector<std::string>>>("TrajHistoNames", {{""}}));

    std::vector<std::string> track_profile_values = config->get<std::vector<std::string>>("TrackProfileValues", {});
    for (const std::string &name: track_profile_values) {
      fTrackProfileValues.push_back(numu::MakeROOTValue(name));
    }

    std::vector<std::string> track_profile_value_names = config->get<std::vector<std::string>>("TrackProfileValueNames", {});
    std::vector<std::tuple<unsigned, float, float>> track_profile_xranges =  config->get<std::vector<std::tuple<unsigned, float, float>>>("TrackProfileXRanges", {});

    assert(track_profile_value_names.size() == track_profile_xranges.size());
    assert(fTrackProfileValues.size() == track_profile_xranges.size());

    if (fDoNormalize) {
      fNormalize.Initialize(config->get<fhicl::ParameterSet>("Normalize", {}));
      fFileTypes = config->get<std::vector<std::string>>("FileTypes");
      fNCosmicData = config->get<double>("NCosmicData", 0.);
    }
    fTrajHistograms.Initialize(fTrajHistoNames);
    fROC.Initialize();
    fRecoEvent = NULL;

    fHistograms.Initialize("",  fTrackSelectorNames, track_profile_value_names, track_profile_xranges); 
    fNeutrinoHistograms.Initialize("Neutrino", fTrackSelectorNames, track_profile_value_names, track_profile_xranges);
    fCosmicHistograms.Initialize("Cosmic", fTrackSelectorNames, track_profile_value_names, track_profile_xranges);

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
    from_file.Get(*f, fTrajHistoNames);
    fTrajHistograms.Add(from_file);
  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    std::cout << "Event: " << core_event->metadata.eventID << std::endl;

    // update each reco Interaction to a smarter primary track selector
    unsigned i = 0;
    while (i < fRecoEvent->reco.size()) {
      // int primary_track = numu::SelectLongestIDdMuon(fRecoEvent->reco_tracks, fRecoEvent->reco[i].slice);
      int primary_track = numu::SelectLongestTrack(fRecoEvent->reco_tracks, fRecoEvent->reco[i].slice);

      // remove vertices without a good primary track
      if (primary_track < 0) {
        fRecoEvent->reco.erase(fRecoEvent->reco.begin() + i); 
        continue;
      }

      fRecoEvent->reco[i].slice.primary_track_index = primary_track;
      fRecoEvent->reco[i].primary_track = fRecoEvent->reco_tracks.at(primary_track);
      // re-do truth matching
      fRecoEvent->reco[i].match = numu::InteractionTruthMatch(fRecoEvent->truth, fRecoEvent->reco_tracks, fRecoEvent->reco[i]);

      numu::RecoTrack &track = fRecoEvent->reco_tracks.at(primary_track); 
      // refine CRT information on primary track
      if (!fCuts.HasCRTHitMatch(track)) {
        track.crt_match.hit.present = false;
      }
      if (!fCuts.HasCRTTrackMatch(track)) {
        track.crt_match.track.present = false;
      }

   
      fRecoEvent->reco[i].primary_track = fRecoEvent->reco_tracks.at(primary_track);
      i += 1;
    }

    for (const numu::RecoInteraction &reco: fRecoEvent->reco) {
      if (reco.primary_track.match.has_match && reco.primary_track.match.mctruth_origin==2/*kCosmicRay*/) {
        if (fRecoEvent->true_tracks.at(reco.primary_track.match.mcparticle_id).wall_enter == 1 /*wTop*/) {
          if (abs(reco.primary_track.start.Y()) < 50. && abs(reco.primary_track.end.Y()) < 50.) {
            std::cout << "Weird Cosmic track\n"; 
            std::cout << "vertex: " << reco.position.X() << " " << reco.position.Y() << " " << reco.position.Z() << std::endl;
            std::cout << "start: " << reco.primary_track.start.X() << " " << reco.primary_track.start.Y() << " " << reco.primary_track.start.Z() << std::endl;
            std::cout << "end: " << reco.primary_track.end.X() << " " << reco.primary_track.end.Y() << " " << reco.primary_track.end.Z() << std::endl;
            std::cout << "length: " << reco.primary_track.length << std::endl;
            std::cout << "True PDG: " << reco.primary_track.match.match_pdg << std::endl;
            std::cout << "True Completion: " << reco.primary_track.match.completion << std::endl; 
            const numu::RecoTrack &true_track = fRecoEvent->true_tracks.at(reco.primary_track.match.mcparticle_id);
            std::cout << "True start: " << true_track.start.X() << " " << true_track.start.Y() << " " << true_track.start.Z() << std::endl;
            std::cout << "True end: " << true_track.end.X() << " " << true_track.end.Y() << " " << true_track.end.Z() << std::endl;
            std::cout << "True length: " << true_track.length << std::endl;
          }
        }
      }
    }

    // make sure no two vertices match to the same true neutrino interaction
    numu::CorrectMultiMatches(*fRecoEvent, fRecoEvent->reco);

    fHistsToFill->Fill(*fRecoEvent, *core_event, fCuts, fTrackSelectors, fTrackProfileValues, fFillAllTracks);
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

