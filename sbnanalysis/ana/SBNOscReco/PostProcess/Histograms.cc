#include <TH1D.h>
#include <TH2D.h>

#include "nusimdata/SimulationBase/MCTruth.h"

#include "Histograms.h"
#include "../Histograms/Derived.h"

namespace ana {
  namespace SBNOsc {

void Histograms::Fill(
  const numu::RecoEvent &event, 
  const event::Event &core, 
  const Cuts &cutmaker, 
  const std::vector<numu::TrackSelector> &selectors, 
  const std::vector<numu::ROOTValue> &xvalues, 
  bool fill_all_tracks) { 

  if (fill_all_tracks) {
    for (const auto &track_pair: event.reco_tracks) {
      const numu::RecoTrack &track = track_pair.second;
      for (unsigned i = 0; i < fAllTracks.size(); i++) {
        bool select = selectors[i](track, event);
        if (select) {
          fAllTracks[i].Fill(track, event.true_tracks);
        }
      }
    }
  }

  for (unsigned i = 0; i < event.reco.size(); i++) {
    std::array<bool, Cuts::nCuts> cuts = cutmaker.ProcessRecoCuts(event, i);
    const numu::RecoInteraction &interaction = event.reco[i];

    if (event.reco_tracks.size() > (unsigned)interaction.slice.primary_track_index) {
      const numu::RecoTrack &track = event.reco_tracks.at(interaction.slice.primary_track_index);
      for (unsigned j = 0; j < fPrimaryTracks.size(); j++) { 
        bool select = selectors[j](track, event);
        if (select) {
          for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
            if (cuts[cut_i]) {
              fPrimaryTracks[j][cut_i].Fill(track, event.true_tracks);
              for (unsigned k = 0; k < xvalues.size(); k++) {
                fPrimaryTrackProfiles[j][k][cut_i].Fill(xvalues[k], track, event);
              }
            }
          }
        }
      }
    }

    // fill histos
    for (size_t cut_i=0; cut_i < Cuts::nCuts; cut_i++) {
      int mode = interaction.match.mode; 
      if (cuts[cut_i]) {
        fInteraction[cut_i+Cuts::nTruthCuts][mode].Fill(i, false, event, core.truth);
        fInteraction[cut_i+Cuts::nTruthCuts][numu::mAll].Fill(i, false, event, core.truth);
      }
    }
  }

  for (unsigned i = 0; i < event.truth.size(); i++) {
    const numu::RecoInteraction &truth = event.truth[i];
    std::array<bool, Cuts::nTruthCuts> cuts = cutmaker.ProcessTruthCuts(event, i); 
    for (int cut_i = 0; cut_i < Cuts::nTruthCuts; cut_i++) {
      int mode = truth.match.mode;
      if (cuts[cut_i]) {
        fInteraction[cut_i][mode].Fill(i, true, event, core.truth);
        fInteraction[cut_i][numu::mAll].Fill(i, true, event, core.truth);
      }
    }
  }
}

void Histograms::Initialize(
  const std::string &prefix, 
  const std::vector<std::string> &track_histo_types, 
  const std::vector<std::string> &track_profile_types,
  const std::vector<std::tuple<unsigned, float, float>> &track_profile_xranges) {

  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      std::string postfix = mode2Str(mode) + prefix + histoNames[i];
      fInteraction[i][mode].Initialize(postfix); 
    }
  }
  fAllTracks.reserve(track_histo_types.size());
  fPrimaryTracks.reserve(track_histo_types.size());

  for (unsigned i = 0; i < track_histo_types.size(); i++) {
    fAllTracks.emplace_back();
    fPrimaryTracks.emplace_back();
    fAllTracks[i].Initialize(prefix + "All_" + track_histo_types[i]);
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      fPrimaryTracks[i][j].Initialize(prefix + "Primary_" + track_histo_types[i] + "_" + std::string(Histograms::histoNames[Cuts::nTruthCuts+j]));

      for (unsigned k = 0; k < track_profile_types.size(); k++) {
        unsigned n_bin;
        float xlo, xhi;
        std::tie(n_bin, xlo, xhi) = track_profile_xranges[k];
        fPrimaryTrackProfiles[i][k][j].Initialize(
          prefix + "Primary_" + track_histo_types[i] + "_" + track_profile_types[k] + "_" + std::string(Histograms::histoNames[Cuts::nTruthCuts+j]),
          n_bin, xlo, xhi);
      } 

    }
  }

  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      fAllHistos.insert(fAllHistos.end(), fInteraction[i][mode].fAllHistos.begin(), fInteraction[i][mode].fAllHistos.end());
    }
  }

  for (unsigned i = 0; i < fAllTracks.size(); i++) {
    fAllHistos.insert(fAllHistos.end(), fAllTracks[i].fAllHistos.begin(), fAllTracks[i].fAllHistos.end());
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      fAllHistos.insert(fAllHistos.end(), fPrimaryTracks[i][j].fAllHistos.begin(), fPrimaryTracks[i][j].fAllHistos.end());
    } 
  }

}

  } // namespace SBNOsc
} // namespace ana
