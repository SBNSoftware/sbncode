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
#include <TParameter.h>
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

    fUseCalorimetry = config->get<bool>("UseCalorimetry", true);

    fHistogramPostfix = config->get<std::string>("HistogramPostfix", "");

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
    }
    fTrajHistograms.Initialize(fTrajHistoNames);
    fROC.Initialize();
    fRecoEvent = NULL;

    fCRTGeo = new sbnd::CRTGeoAlg(fProviderManager->GetGeometryProvider(), fProviderManager->GetAuxDetGeometryProvider());
    if (fDoNormalize) {
      fCosmicHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, "Cosmic", fTrackSelectorNames, track_profile_value_names, track_profile_xranges);
    }
    fHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, fHistogramPostfix,  fTrackSelectorNames, track_profile_value_names, track_profile_xranges); 
  }

  void Selection::FileSetup(TFile *f, TTree *eventTree) {
    eventTree->SetBranchAddress("reco_event", &fRecoEvent);
    fCuts.Initialize(fCutConfig, fProviderManager->GetGeometryProvider());
    TParameter<int> *file_type = (TParameter<int> *)f->Get("MCType");
    fFileType = (numu::MCType) file_type->GetVal();

    if (fDoNormalize) {
      if (fFileType == numu::fIntimeCosmic) {
        fHistsToFill = &fCosmicHistograms;
      }
      else if (fFileType == numu::fOverlay) {
        fHistsToFill = &fHistograms;
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

    CRTHistos crt_from_file;
    crt_from_file.Get(*f, "_all");
    if (fFileType == numu::fIntimeCosmic) {
      fCRTCosmicHistos.Add(crt_from_file);
    }
    else if (fFileType == numu::fOverlay) {
      fCRTNeutrinoHistos.Add(crt_from_file);
    }
  }

  void Selection::ProcessSubRun(const SubRun *subrun) {
    if (fDoNormalize && fFileType == numu::fOverlay) {
      fNormalize.AddNeutrinoSubRun(*subrun);
    }
  }

  void Selection::ProcessFileMeta(const FileMeta *meta) {
    if (fDoNormalize && fFileType == numu::fIntimeCosmic) {
      fNormalize.AddCosmicFile(*meta);
    }
  }

  void Selection::ProcessEvent(const event::Event *core_event) {
    std::cout << "Event: " << core_event->metadata.eventID << std::endl;

    // update each reco Interaction to a smarter primary track selector
    unsigned i = 0;
    while (i < fRecoEvent->reco.size()) {
      for (size_t particle_ind: fRecoEvent->reco[i].slice.tracks) {
        numu::RecoTrack &track = fRecoEvent->tracks.at(particle_ind);

        // set containment
        if (fCuts.InCalorimetricContainment(track.start) && fCuts.InCalorimetricContainment(track.end)) {
          track.is_contained = true;
        }
        else {
          track.is_contained = false;
        }
      }

      int primary_track;
      if (fUseCalorimetry) {
        primary_track = numu::SelectLongestIDdMuon(fRecoEvent->tracks, fRecoEvent->reco[i].slice);
      }
      else {
        primary_track = numu::SelectLongestTrack(fRecoEvent->tracks, fRecoEvent->reco[i].slice);
      }

      // remove vertices without a good primary track
      if (primary_track < 0) {
        fRecoEvent->reco.erase(fRecoEvent->reco.begin() + i); 
        continue;
      }

      fRecoEvent->reco[i].slice.primary_track_index = primary_track;
      fRecoEvent->reco[i].primary_track = fRecoEvent->tracks.at(primary_track);

      // re-do truth matching
      fRecoEvent->reco[i].match = numu::InteractionTruthMatch(core_event->truth, fRecoEvent->tracks, fRecoEvent->reco[i]);

      // if this is an in-time cosmic file, update the cosmic mode
      if (fFileType == numu::fIntimeCosmic) {
        assert(fRecoEvent->reco[i].match.mode == numu::mCosmic || fRecoEvent->reco[i].match.mode == numu::mOther);
        fRecoEvent->reco[i].match.mode = numu::mIntimeCosmic;
      }

      i += 1;
    }

    fROC.Fill(fCuts, *fRecoEvent, fFileType == numu::fIntimeCosmic);

    // make sure no two vertices match to the same true neutrino interaction
    numu::CorrectMultiMatches(*fRecoEvent, fRecoEvent->reco);

    fHistsToFill->Fill(*fRecoEvent, *core_event, fCuts, fTrackSelectors, fTrackProfileValues, fFillAllTracks);
    if (fDoNormalize) {
      if (fFileType == numu::fIntimeCosmic) fNormalize.AddCosmicEvent(*core_event);
      else fNormalize.AddNeutrinoEvent(*core_event);
    }

  }

  void Selection::Finalize() {
    fOutputFile->cd();
    std::cout << "Finalizing!\n";
    if (fDoNormalize) {
      std::cout << "Normalizing!\n";
      // save scale per 1e20 POT
      TParameter<float> scale_neutrino("ScaleNeutrino", fNormalize.ScaleNeutrino(1.e20));
      TParameter<float> scale_cosmic("ScaleCosmic", fNormalize.ScaleCosmic(1.e20));
      scale_neutrino.Write();
      scale_cosmic.Write();
      std::cout << "Scale Neutrino: " << fNormalize.ScaleNeutrino(1.e20) << std::endl;
      std::cout << "Scale Cosmic: " << fNormalize.ScaleCosmic(1.e20) << std::endl;

      fHistograms.Scale(fNormalize.ScaleNeutrino(fGoalPOT));
      fCosmicHistograms.Scale(fNormalize.ScaleCosmic(fGoalPOT));

      fCRTCosmicHistos.Scale(fNormalize.ScaleNeutrino(fGoalPOT));
      fCRTNeutrinoHistos.Scale(fNormalize.ScaleCosmic(fGoalPOT));

      fHistograms.Add(fCosmicHistograms);
      fROC.Normalize(fNormalize.ScaleNeutrino(fGoalPOT), fNormalize.ScaleCosmic(fGoalPOT));  
    }
    fROC.BestCuts();
    fROC.Write();
    fHistograms.Write();
    fTrajHistograms.Write();
    fCRTCosmicHistos.Write();
    fCRTNeutrinoHistos.Write();

  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

