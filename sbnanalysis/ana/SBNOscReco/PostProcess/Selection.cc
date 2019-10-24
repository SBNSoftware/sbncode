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

    // check for no two reco matched to same truth
    std::map<int, unsigned> matched_vertices;
    for (unsigned reco_i = 0; reco_i < fRecoEvent->reco.size(); reco_i++) {
      numu::RecoInteraction &reco = fRecoEvent->reco[reco_i];
      unsigned set_i = reco_i; 

      if (reco.primary_track.match.has_match && reco.primary_track.match.is_primary && reco.match.event_track_id >= 0 && reco.match.event_track_id < fRecoEvent->truth.size()-1) {
        if (matched_vertices.count(reco.match.event_track_id)) {
          std::cout << "BAAAAAAAAD: two reco matched to same truth." << std::endl;
          std::cout << "This id: " << reco.match.event_track_id << " N ids: " << fRecoEvent->reco.size() << " N truth: " << core_event->truth.size() << std::endl;
          for (const numu::RecoInteraction &reco: fRecoEvent->reco) {
            std::cout << "Interaction x: " << reco.position.X() << " match: " << reco.match.event_track_id << " primary track pdg: " << reco.primary_track.match.match_pdg << std::endl;
          }

          numu::RecoInteraction &other = fRecoEvent->reco[matched_vertices[reco.match.event_track_id]];

          float dist_reco = std::min( (reco.primary_track.start - fRecoEvent->true_tracks[reco.primary_track.match.mcparticle_id].start).Mag(),
                                      (reco.primary_track.end   - fRecoEvent->true_tracks[reco.primary_track.match.mcparticle_id].start).Mag());

          float dist_other = std::min( (other.primary_track.start - fRecoEvent->true_tracks[other.primary_track.match.mcparticle_id].start).Mag(),
                                      (other.primary_track.end   - fRecoEvent->true_tracks[other.primary_track.match.mcparticle_id].start).Mag());

          // Handle if one of the matched particles in not a muon
          if (abs(reco.primary_track.match.match_pdg) != 13 && abs(other.primary_track.match.match_pdg) == 13) {
            set_i = matched_vertices[reco.match.event_track_id];
            reco.primary_track.match.is_primary = false; 
            std::cout << "Corrected -- matched track is non-muon\n";
          } 
          else if (abs(reco.primary_track.match.match_pdg) == 13 && abs(other.primary_track.match.match_pdg) != 13) {
            other.primary_track.match.is_primary = false; 
            std::cout << "Corrected -- matched track is non-muon\n";
          } 
          // Handle both two different particles neither of which are muons
          else if (reco.primary_track.match.mcparticle_id != other.primary_track.match.mcparticle_id) {
            std::cout << "Reco track 1 length: " << reco.primary_track.length << " pdg: " << reco.primary_track.match.match_pdg << std::endl; 
            std::cout << "Reco track 2 length: " << other.primary_track.length << " pdg: " << other.primary_track.match.match_pdg << std::endl; 
            if (reco.primary_track.length > other.primary_track.length) {
              other.primary_track.match.is_primary = false;
              std::cout << "Corrected -- used longer track. Track 1 is primary.\n";
            }
            else {
              reco.primary_track.match.is_primary = false;
              std::cout << "Corrected -- used longer track. Track 2 is primary.\n";
            }
          }
          // Handle if the priamry track was split in two
          else if (reco.primary_track.match.mcparticle_id == other.primary_track.match.mcparticle_id) {
            std::cout << "Reco track 1 pos: " << reco.primary_track.start.X() << " " << reco.primary_track.start.Y() << " " << reco.primary_track.start.Z() << " to: " << reco.primary_track.end.X() << " " << reco.primary_track.end.Y() << " " << reco.primary_track.end.Z() << std::endl;
            std::cout << "Reco track 2 pos: " << other.primary_track.start.X() << " " << other.primary_track.start.Y() << " " << other.primary_track.start.Z() << " to: " << other.primary_track.end.X() << " " << other.primary_track.end.Y() << " " << other.primary_track.end.Z() << std::endl;
            TVector3 match_start = fRecoEvent->true_tracks[other.primary_track.match.mcparticle_id].start;
            std::cout << "Match start pos: " << match_start.X() << " " << match_start.Y() << " " << match_start.Z() << std::endl;
            if (dist_reco <= dist_other) {
              std::cout << "Corrected -- split muon. Track 1 is primary.\n";
              other.primary_track.match.is_primary = false; 
            }
            else {
              std::cout << "Corrected -- split muon. Track 2 is primary.\n";
              set_i = matched_vertices[reco.match.event_track_id];
              reco.primary_track.match.is_primary = false; 
            }
          }
          else { 
            std::cout << "Unable to correct!\n";
            std::cerr << "Exiting on failure.\n";
            assert(false);
          }
        }
        matched_vertices[reco.match.event_track_id] = set_i;
      } 
    }

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

