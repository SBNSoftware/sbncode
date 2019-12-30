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
    fRecoEventReaders = std::vector<TTreeReaderValue<numu::RecoEvent>*>(config->get<unsigned>("NThreads"), NULL);
    std::cerr << "N readers: " << fRecoEventReaders.size() << std::endl;

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

    fCRTGeo = new sbnd::CRTGeoAlg(fProviderManager->GetGeometryProvider(), fProviderManager->GetAuxDetGeometryProvider());
    if (fDoNormalize) {
      fCosmicHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, "Cosmic", fTrackSelectorNames, track_profile_value_names, track_profile_xranges);
    }
    fHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, fHistogramPostfix,  fTrackSelectorNames, track_profile_value_names, track_profile_xranges); 
  }

  void Selection::FileSetup(TFile *f) {
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

  void Selection::EventTreeSetup(TTreeReader &reader, unsigned thread_ind) {
    if (fRecoEventReaders[thread_ind] != NULL) delete fRecoEventReaders[thread_ind];
    fRecoEventReaders[thread_ind] = new TTreeReaderValue<numu::RecoEvent>(reader, "reco_event");
  }

  void Selection::ProcessSubRun(const SubRun *subrun) {
    if (fDoNormalize && fFileType == numu::fOverlay) {
      fNormalize.AddNeutrinoSubRun(*subrun);
    }
  }

  void Selection::ProcessFileMeta(const FileMeta *meta) {
    if (fDoNormalize) {
      if (fFileType == numu::fIntimeCosmic) {
        fNormalize.AddCosmicFile(*meta);
      }
      else {
        fNormalize.AddNeutrinoFile(*meta);
      }
    }
  }

  void Selection::ProcessEvent(const event::Event *core_event, unsigned thread_ind) {
    std::cout << "Event: " << core_event->metadata.eventID << std::endl;
    numu::RecoEvent *reco_event = fRecoEventReaders[thread_ind]->Get();

    // update each reco Interaction to a smarter primary track selector
    unsigned i = 0;
    while (i < reco_event->reco.size()) {
      for (size_t ind: reco_event->reco[i].slice.tracks) {
        if (ind >= reco_event->reco_tracks.size()) continue;
   
        // Set final momentum values
        numu::RecoTrack &track = reco_event->reco_tracks.at(ind);
        // TODO: apply calorimetry, set for non-primary tracks
        track.range_momentum = track.range_momentum_muon;
        // TODO: use forward/backward?
        track.mcs_momentum = track.mcs_muon.fwd_mcs_momentum;
        if (fCuts.InCalorimetricContainment(track.start) && fCuts.InCalorimetricContainment(track.end)) {
          track.momentum = track.range_momentum;
          track.is_contained = true;
        }
        else {
          track.momentum = track.mcs_momentum;
          track.is_contained = false;
        }
      }

      int primary_track;
      if (fUseCalorimetry) {
        primary_track = numu::SelectLongestIDdMuon(reco_event->reco_tracks, reco_event->reco[i].slice);
      }
      else {
        primary_track = numu::SelectLongestTrack(reco_event->reco_tracks, reco_event->reco[i].slice);
      }

      // remove vertices without a good primary track
      if (primary_track < 0) {
        reco_event->reco.erase(reco_event->reco.begin() + i); 
        continue;
      }

      reco_event->reco[i].slice.primary_track_index = primary_track;
      reco_event->reco[i].primary_track = reco_event->reco_tracks.at(primary_track);
      // re-do truth matching
      reco_event->reco[i].match = numu::InteractionTruthMatch(reco_event->truth, reco_event->reco_tracks, reco_event->reco[i]);

      // if this is an in-time cosmic file, update the cosmic mode
      if (fFileType == numu::fIntimeCosmic) {
        assert(reco_event->reco[i].match.mode == numu::mCosmic || reco_event->reco[i].match.mode == numu::mOther);
        reco_event->reco[i].match.mode = numu::mIntimeCosmic;
      }

      i += 1;
    }

    // if this is an in-time cosmic file, update the cosmic mode
    if (fFileType == numu::fIntimeCosmic) {
      assert(reco_event->truth.size() == 1);
      reco_event->truth[0].match.mode = numu::mIntimeCosmic;
    }

    fROC.Fill(fCuts, *reco_event, fFileType == numu::fIntimeCosmic);

    for (const numu::RecoInteraction &reco: reco_event->reco) {
      if (reco.primary_track.match.has_match && reco.primary_track.match.mctruth_origin==1) { // kBeamNeutrino
        if (reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).wall_enter == 0 &&
            reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).wall_exit == 0) { // wNone
          if (abs(reco.primary_track.end.Y() - 200.) < 2. || abs(reco.primary_track.start.Y() - 200.) < 2.) {
            std::cout << "Weird CC Track\n";
            std::cout << "start: " << reco.primary_track.start.X() << " " << reco.primary_track.start.Y() << " " << reco.primary_track.start.Z() << std::endl;
            std::cout << "end: " << reco.primary_track.end.X() << " " << reco.primary_track.end.Y() << " " << reco.primary_track.end.Z() << std::endl;
            std::cout << "true start: " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).start.X() << " " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).start.Y() << " " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).start.Z() << std::endl;
            std::cout << "true end: " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).end.X() << " " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).end.Y() << " " << reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).end.Z() << std::endl;
          }
        }
      }
      /*
      if (reco.primary_track.match.has_match && reco.primary_track.match.mctruth_origin==2) { // kComsicRay
        if (reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id).wall_enter == 1) { // wTop
          if (abs(reco.primary_track.start.Y()) < 50. && abs(reco.primary_track.end.Y()) < 50.) {
            std::cout << "Weird Cosmic track\n"; 
            std::cout << "vertex: " << reco.position.X() << " " << reco.position.Y() << " " << reco.position.Z() << std::endl;
            std::cout << "start: " << reco.primary_track.start.X() << " " << reco.primary_track.start.Y() << " " << reco.primary_track.start.Z() << std::endl;
            std::cout << "end: " << reco.primary_track.end.X() << " " << reco.primary_track.end.Y() << " " << reco.primary_track.end.Z() << std::endl;
            std::cout << "length: " << reco.primary_track.length << std::endl;
            std::cout << "True PDG: " << reco.primary_track.match.match_pdg << std::endl;
            std::cout << "True Completion: " << reco.primary_track.match.completion << std::endl; 
            const numu::RecoTrack &true_track = reco_event->true_tracks.at(reco.primary_track.match.mcparticle_id);
            std::cout << "True start: " << true_track.start.X() << " " << true_track.start.Y() << " " << true_track.start.Z() << std::endl;
            std::cout << "True end: " << true_track.end.X() << " " << true_track.end.Y() << " " << true_track.end.Z() << std::endl;
            std::cout << "True length: " << true_track.length << std::endl;
          }
        }
      }*/
    }

    // further refine reco stuff
    /*
    for (numu::RecoInteraction &reco: reco_event->reco) {
      numu::RecoTrack &track = reco_event->reco_tracks.at(reco.slice.primary_track_index); 
      // refine CRT information on primary track
      if (!fCuts.HasCRTHitMatch(track)) {
        track.crt_match.hit.present = false;
      }
      if (!fCuts.HasCRTTrackMatch(track)) {
        track.crt_match.track.present = false;
      }

      reco.primary_track = track;
    }*/

    // make sure no two vertices match to the same true neutrino interaction
    numu::CorrectMultiMatches(*reco_event, reco_event->reco);

    fHistsToFill->Fill(*reco_event, *core_event, fCuts, fTrackSelectors, fTrackProfileValues, fFillAllTracks);
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

