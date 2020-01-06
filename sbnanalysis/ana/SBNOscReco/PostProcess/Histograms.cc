#include <TH1D.h>
#include <TH2D.h>

#include "nusimdata/SimulationBase/MCTruth.h"

#include "Histograms.h"
#include "../Histograms/Derived.h"
#include "../RecoUtils/GeoUtil.h"
#include "../NumuReco/TruthMatch.h"

namespace ana {
  namespace SBNOsc {

void Histograms::Fill( const numu::RecoEvent &event, 
		       const event::Event &core, 
		       const Cuts &cutmaker, 
		       const std::vector<numu::TrackSelector> &selectors, 
                       const std::vector<numu::TrackFunction> &xfunctions,
		       bool fill_all_tracks) { 

  if (event.type == numu::fOverlay) {
    fCosmic[0].Fill(event.particles);
    if (cutmaker.PassFlashTrigger(event)) {
      fCosmic[1].Fill(event.particles);
    }
  }
  else {
    std::cout << "Filling Cosmic!\n";
    fCosmic[2].Fill(event.particles);
    if (cutmaker.PassFlashTrigger(event)) {
      fCosmic[3].Fill(event.particles);
    }
  }

  numu::TrueParticle bad;

  if (fill_all_tracks) {
    for (const auto &track_pair: event.tracks) {
      const numu::RecoTrack &track = track_pair.second;
      for (unsigned i = 0; i < fAllTracks.size(); i++) {
        const numu::TrueParticle &part = (track.match.has_match) ? event.particles.at(track.match.mcparticle_id) : bad;
        bool select = selectors[i](track, part, event.type);
        if (select) {
          fAllTracks[i].Fill(track, event.particles);
        }
      }
    }
  }

  for (unsigned i = 0; i < event.reco.size(); i++) {
    std::array<bool, Cuts::nCuts> cuts = cutmaker.ProcessRecoCuts(event, i);
    const numu::RecoInteraction &interaction = event.reco[i];


    if (event.tracks.size() > (unsigned)interaction.slice.primary_track_index) {
      const numu::RecoTrack &track = event.tracks.at(interaction.slice.primary_track_index);

      for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
        if (cuts[cut_i] && cutmaker.HasCRTHitMatch(track)) {
          std::cout << "Filling CRT Histo!\n";
          fCRTs[cut_i].Fill(track.crt_match.hit);
        }
      }

      for (unsigned j = 0; j < fPrimaryTracks.size(); j++) { 
        const numu::TrueParticle &part = (track.match.has_match) ? event.particles.at(track.match.mcparticle_id) : bad;
        bool select = selectors[i](track, part, event.type);
        if (select) {
          for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
            if (cuts[cut_i]) {
              fPrimaryTracks[j][cut_i].Fill(track, event.particles);
              for (unsigned k = 0; k < xfunctions.size(); k++) {
                const numu::TrueParticle &part = (track.match.has_match) ? event.particles.at(track.match.mcparticle_id) : bad;
                uscript::Value x = xfunctions[k](&track, &part, (unsigned*)&event.type);
                assert(IS_NUMBER(x));
                fPrimaryTrackProfiles[j][k][cut_i].Fill(AS_NUMBER(x), 
							track, event);
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
        fInteraction[cut_i+Cuts::nTruthCuts][mode].Fill(event.reco[i], 
							event, 
							core.truth);
        fInteraction[cut_i+Cuts::nTruthCuts][numu::mAll].Fill(event.reco[i], 
							      event, 
							      core.truth);
      }
    }
  }

  for (unsigned i = 0; i < core.truth.size(); i++) {
    std::array<bool, Cuts::nTruthCuts> cuts = cutmaker.ProcessTruthCuts(event, core, i); 
    for (int cut_i = 0; cut_i < Cuts::nTruthCuts; cut_i++) {
      int mode = numu::GetMode(core.truth[i]);
      if (cuts[cut_i]) {
        fInteraction[cut_i][mode].Fill(core.truth[i], i, event);
        fInteraction[cut_i][numu::mAll].Fill(core.truth[i], i, event);
      }
    }
  }
}

void Histograms::Initialize(
  const geo::GeometryCore *geometry,
  const sbnd::CRTGeoAlg &crt_geo,
  const std::string &prefix, 
  const std::vector<std::string> &track_histo_types, 
  const std::vector<std::string> &track_profile_types,
  const std::vector<std::tuple<unsigned, float, float>> &track_profile_xranges) {

  double max_length = SBNRecoUtils::MaxLength(geometry);
  geo::BoxBoundedGeo detector = SBNRecoUtils::DetectorVolume(geometry);
  std::vector<double> tagger_volume = crt_geo.CRTLimits();

  std::cout << "Limits: " << tagger_volume[0] << " " << tagger_volume[3] << " " << tagger_volume[1] << " " << tagger_volume[4] << " " << tagger_volume[2] << " " << tagger_volume[5] << std::endl;

  fCosmic[0].Initialize(prefix + "_Cosmic", detector);
  fCosmic[1].Initialize(prefix + "_CosmicTrig", detector);
  fCosmic[2].Initialize(prefix + "_IntimeCosmic", detector);
  fCosmic[3].Initialize(prefix + "_IntimeCosmicTrig", detector);

  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      std::string postfix = mode2Str(mode) + prefix + histoNames[i];
      fInteraction[i][mode].Initialize(postfix, detector, tagger_volume); 
    }
  }

  for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
    fCRTs[cut_i].Initialize(prefix + "_crt_cut_" + std::string(Histograms::histoNames[Cuts::nTruthCuts+cut_i]), tagger_volume); 
  }

  fAllTracks.reserve(track_histo_types.size());
  fPrimaryTracks.reserve(track_histo_types.size());

  for (unsigned i = 0; i < track_histo_types.size(); i++) {
    fAllTracks.emplace_back();
    fPrimaryTracks.emplace_back();
    fPrimaryTrackProfiles.emplace_back();
    fAllTracks[i].Initialize(prefix + "All_" + 
			     track_histo_types[i], 
			     detector, max_length);
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      fPrimaryTracks[i][j].Initialize(prefix + "Primary_" + 
				      track_histo_types[i] + 
				      "_" + 
				      std::string(Histograms::histoNames[Cuts::nTruthCuts+j]), detector, max_length);

      for (unsigned k = 0; k < track_profile_types.size(); k++) {
        if (j == 0) fPrimaryTrackProfiles[i].emplace_back();
        unsigned n_bin;
        float xlo, xhi;
        std::tie(n_bin, xlo, xhi) = track_profile_xranges[k];
        fPrimaryTrackProfiles[i][k][j].Initialize(prefix + 
						  "Primary_" + 
						  track_histo_types[i] + 
						  "_" + 
						  track_profile_types[k] + 
						  "_" + 
						  std::string(Histograms::histoNames[Cuts::nTruthCuts+j]),
						  n_bin, xlo, xhi);
      } 

    }
  }


  for (unsigned i = 0; i < 4; i++) {
    fAllHistos.insert(fAllHistos.end(), 
		      fCosmic[i].fAllHistos.begin(), 
		      fCosmic[i].fAllHistos.end());

  }
  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      fAllHistos.insert(fAllHistos.end(), 
			fInteraction[i][mode].fAllHistos.begin(), 
			fInteraction[i][mode].fAllHistos.end());
    }
  }

  for (unsigned i = 0; i < Cuts::nCuts; i++) {
    fAllHistos.insert(fAllHistos.end(), fCRTs[i].fAllHistos.begin(), fCRTs[i].fAllHistos.end());
  }

  for (unsigned i = 0; i < fAllTracks.size(); i++) {
    fAllHistos.insert(fAllHistos.end(), 
		      fAllTracks[i].fAllHistos.begin(), 
		      fAllTracks[i].fAllHistos.end());
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      fAllHistos.insert(fAllHistos.end(), 
			fPrimaryTracks[i][j].fAllHistos.begin(), 
			fPrimaryTracks[i][j].fAllHistos.end());
      for (unsigned k = 0; k < track_profile_types.size(); k++) {
        fAllHistos.insert(fAllHistos.end(), 
			  fPrimaryTrackProfiles[i][k][j].fAllHistos.begin(), 
			  fPrimaryTrackProfiles[i][k][j].fAllHistos.end());
      } 
    } 
  }

}

  } // namespace SBNOsc
} // namespace ana
