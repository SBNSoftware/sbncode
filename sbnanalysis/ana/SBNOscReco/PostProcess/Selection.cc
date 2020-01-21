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
#include "SetEvent.h"
#include "../uScript/api.h"
#include "../uScript/value.h"

namespace ana {
namespace SBNOsc {
  void Selection::Initialize(fhicl::ParameterSet *config) {
    std::cout << "Input Configuration:\n";
    std::cout << config->to_indented_string() << std::endl;
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

    fTrackSelectorNames = numu::MultiplyNames(track_selector_names, '/');
    fTrackSelectors = numu::MultiplyTrackSelectors(track_selector_strings);

    assert(fTrackSelectorNames.size() == fTrackSelectors.size());

    std::vector<std::string> track_profile_values = config->get<std::vector<std::string>>("TrackProfileValues", {});
    for (const std::string &source: track_profile_values) {
      fTrackProfileValues.push_back(uscript::compile<numu::RecoTrack, numu::TrueParticle, unsigned>("track", "particle", "mctype", source.c_str()));
    }

    std::vector<std::string> track_profile_value_names = config->get<std::vector<std::string>>("TrackProfileValueNames", {});
    std::vector<std::tuple<unsigned, float, float>> track_profile_xranges =  config->get<std::vector<std::tuple<unsigned, float, float>>>("TrackProfileXRanges", {});

    assert(track_profile_value_names.size() == track_profile_xranges.size());
    assert(fTrackProfileValues.size() == track_profile_xranges.size());

    if (fDoNormalize) {
      fNormalize.Initialize(config->get<fhicl::ParameterSet>("Normalize", {}));
    }
    fROC.Initialize();
    fRecoEvent = NULL;

    fCRTGeo = new sbnd::CRTGeoAlg(fProviderManager->GetGeometryProvider(), fProviderManager->GetAuxDetGeometryProvider());

    fCuts.Initialize(fCutConfig, fProviderManager->GetGeometryProvider());

    if (fDoNormalize) {
      fCosmicHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, fCuts, "Cosmic", fTrackSelectorNames, track_profile_value_names, track_profile_xranges);
    }
    fHistograms.Initialize(fProviderManager->GetGeometryProvider(), *fCRTGeo, fCuts, fHistogramPostfix,  fTrackSelectorNames, track_profile_value_names, track_profile_xranges); 
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

    // set stuff in the event
    SetEvent(*fRecoEvent, *core_event, fCuts, fFileType, fUseCalorimetry);

    for (const numu::RecoInteraction &reco: fRecoEvent->reco) {
      std::cout << "Event type: " << Histograms::mode2Str(reco.match.mode) << std::endl;
      const std::vector<size_t> &tracks = reco.slice.tracks;
      for (size_t ind: tracks) {
        const numu::RecoTrack &track = fRecoEvent->tracks.at(ind);
        const numu::RecoParticle &t_particle = reco.slice.particles.at(ind);
        std::cout << "Reco track: "; 
        std::cout << "length: " << track.length << " PID: " << track.match.match_pdg << " Chi2 m: " << track.chi2_muon << " Chi2 p: " << track.chi2_proton << " n daughter: " << t_particle.daughters.size() << std::endl;
        for (size_t d_ind: t_particle.daughters) {
          if (fRecoEvent->tracks.count(d_ind)) {
            std::cout << "Daughter track: ";
            const numu::RecoTrack &d_track = fRecoEvent->tracks.at(ind);
            const numu::RecoParticle &d_t_particle = reco.slice.particles.at(ind);
            std::cout << "length: " << d_track.length << " PID: " << d_track.match.match_pdg << " Chi2 m: " << d_track.chi2_muon << " Chi2 p: " << d_track.chi2_proton <<" n daughter: " << d_t_particle.daughters.size() << std::endl;
          }
          else {
            std::cout << "Shower Daughter\n";
          }
        }
      }
    }

    fROC.Fill(fCuts, *fRecoEvent, fFileType == numu::fIntimeCosmic);
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
    fCRTCosmicHistos.Write();
    fCRTNeutrinoHistos.Write();

  }
  }   // namespace SBNOsc
}   // namespace ana

DECLARE_SBN_POSTPROCESSOR(ana::SBNOsc::Selection);

