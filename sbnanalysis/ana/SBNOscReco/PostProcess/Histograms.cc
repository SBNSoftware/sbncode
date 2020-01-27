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
  // Fill the particle histos (only with CC events)
  for (const event::Interaction &truth: core.truth) {
    if (truth.neutrino.iscc) {
      for (const event::FinalStateParticle &p: truth.finalstate) {
        if (event.particles.count(p.G4ID)) {
          const numu::TrueParticle &particle = event.particles.at(p.G4ID);

          bool is_primary = particle.start_process == numu::primary;
          if (!is_primary) continue;

          int pdg_ind = -1;
          if (abs(particle.pdgid) == 13) pdg_ind = 0;
          else if (abs(particle.pdgid) == 211) pdg_ind = 1;
          else if (abs(particle.pdgid) == 2212) pdg_ind = 2;
          else if (abs(particle.pdgid) == 321) pdg_ind = 3;
          else continue;

          bool has_reco = false;
          for (auto const &track_pair: event.tracks) {
            if (track_pair.second.match.mcparticle_id == p.G4ID) {
              has_reco = true;
              break;
            }
          }

          bool cc_reco = false;
          for (const numu::RecoInteraction &reco: event.reco) {
            if (reco.slice.match.mode != numu::mCC) continue; 
            for (size_t ind: reco.slice.tracks) {
              if (event.tracks.at(ind).match.mcparticle_id == p.G4ID) {
                cc_reco = true;
                break;
              }
            }
            if (cc_reco) break;
          }
          
          bool cc_np_reco = false;
          for (const numu::RecoInteraction &reco: event.reco) {
            if (reco.slice.match.mode != numu::mCCNonPrimary) continue; 
            for (size_t ind: reco.slice.tracks) {
              if (event.tracks.at(ind).match.mcparticle_id == p.G4ID) {
                cc_np_reco = true;
                break;
              }
            }
            if (cc_np_reco) break;
          }

          bool cosmic_reco = false;
          for (const numu::RecoInteraction &reco: event.reco) {
            if (reco.slice.match.mode != numu::mCosmic) continue; 
            for (size_t ind: reco.slice.tracks) {
              if (event.tracks.at(ind).match.mcparticle_id == p.G4ID) {
                cosmic_reco = true;
                break;
              }
            }
            if (cosmic_reco) break;
          }

          fParticles[pdg_ind][0].Fill(particle); 
          if (has_reco) {
            fParticles[pdg_ind][1].Fill(particle); 
          }
          if (cc_reco) {
            fParticles[pdg_ind][2].Fill(particle); 
          }
          else if (cc_np_reco) {
            fParticles[pdg_ind][3].Fill(particle); 
          }
          else if (cosmic_reco) {
            fParticles[pdg_ind][4].Fill(particle); 
          }

        }
      }
    }
  }

  numu::TrueParticle bad;

  if (fill_all_tracks) {
    for (const auto &track_pair: event.tracks) {
      const numu::RecoTrack &track = track_pair.second;
      const numu::TrueParticle &part = (track.match.has_match) ? event.particles.at(track.match.mcparticle_id) : bad;

      for (unsigned i = 0; i < fAllTracks.size(); i++) {
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


    if (event.tracks.size() > (unsigned)interaction.primary_track_index) {
      const numu::RecoTrack &track = event.tracks.at(interaction.primary_track_index);
      const numu::TrueParticle &part = (track.match.has_match) ? event.particles.at(track.match.mcparticle_id) : bad;

      for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
        if (cuts[cut_i] && cutmaker.HasCRTHitMatch(track)) {
          fCRTs[cut_i].Fill(track.crt_match.hit);
        }
      }

      for (unsigned j = 0; j < fPrimaryTracks.size(); j++) { 
        bool select = selectors[j](track, part, event.type);
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
      int mode = interaction.slice.match.mode; 
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
  const Cuts &cuts,
  const std::string &prefix, 
  const std::vector<std::string> &track_histo_types, 
  const std::vector<std::string> &track_profile_types,
  const std::vector<std::tuple<unsigned, float, float>> &track_profile_xranges) {

  double max_length = SBNRecoUtils::MaxLength(geometry);
  geo::BoxBoundedGeo detector = SBNRecoUtils::DetectorVolume(geometry);
  std::vector<double> tagger_volume = crt_geo.CRTLimits();

  std::cout << "Limits: " << tagger_volume[0] << " " << tagger_volume[3] << " " << tagger_volume[1] << " " << tagger_volume[4] << " " << tagger_volume[2] << " " << tagger_volume[5] << std::endl;

  // make new directory for histograms
  TDirectory *d_top = gDirectory->mkdir("histograms");

  if (prefix.size() != 0) {
    d_top = d_top->mkdir(prefix.c_str());
    d_top->cd();
  }

  d_top->mkdir("cosmic");

  d_top->mkdir("cosmic/outtime");
  d_top->cd("cosmic/outtime");
  fCosmic[0].Initialize("", detector);
  d_top->cd();

  d_top->mkdir("cosmic/outtime_trig");
  d_top->cd("cosmic/outtime_trig");
  fCosmic[1].Initialize("", detector);
  d_top->cd();

  d_top->mkdir("cosmic/intime");
  d_top->cd("cosmic/intime");
  fCosmic[2].Initialize("", detector);
  d_top->cd();

  d_top->mkdir("cosmic/intime_trig");
  d_top->cd("cosmic/intime_trig");
  fCosmic[3].Initialize("", detector);
  d_top->cd();

  // true-particle histograms
  d_top->mkdir("particle");
  for (unsigned i = 0; i < Histograms::nPIDs; i++) {
    std::string all = "particle/" + std::string(Histograms::allPIDs[i]) + "/all";
    std::string has_reco = "particle/" + std::string(Histograms::allPIDs[i]) + "/has_reco";
    std::string cc_reco = "particle/" + std::string(Histograms::allPIDs[i]) + "/cc_reco";
    std::string cc_np_reco = "particle/" + std::string(Histograms::allPIDs[i]) + "/cc_np_reco";
    std::string cosmic_reco = "particle/" + std::string(Histograms::allPIDs[i]) + "/cosmic_reco";

    d_top->mkdir(all.c_str());
    d_top->cd(all.c_str());
    fParticles[i][0].Initialize("");
    d_top->cd();

    d_top->mkdir(has_reco.c_str());
    d_top->cd(has_reco.c_str());
    fParticles[i][1].Initialize("");
    d_top->cd();

    d_top->mkdir(cc_reco.c_str());
    d_top->cd(cc_reco.c_str());
    fParticles[i][2].Initialize("");
    d_top->cd();

    d_top->mkdir(cc_np_reco.c_str());
    d_top->cd(cc_np_reco.c_str());
    fParticles[i][3].Initialize("");
    d_top->cd();

    d_top->mkdir(cosmic_reco.c_str());
    d_top->cd(cosmic_reco.c_str());
    fParticles[i][4].Initialize("");
    d_top->cd();
  }

  std::vector<std::string> cut_order = cuts.CutOrder();
  std::vector<std::string> truth_cut_order = cuts.TruthCutOrder();

  d_top->mkdir("interaction");
  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      std::string cut_name = (i < truth_cut_order.size()) ? truth_cut_order[i] : cut_order[i - truth_cut_order.size()];
      std::string postfix = mode2Str(mode) + prefix + cut_name;
      std::string dirname = "interaction/" + mode2Str(mode) + "/" + cut_name;
      d_top->mkdir(dirname.c_str());
      d_top->cd(dirname.c_str());
      fInteraction[i][mode].Initialize("", detector, tagger_volume); 
      d_top->cd();
    }
  }

  d_top->mkdir("crt");
  for (unsigned cut_i = 0; cut_i < Cuts::nCuts; cut_i++) {
    std::string dirname = "crt/" + cut_order[cut_i];
    d_top->mkdir(dirname.c_str());
    d_top->cd(dirname.c_str());
    fCRTs[cut_i].Initialize("", tagger_volume);
    d_top->cd();
  }

  fAllTracks.reserve(track_histo_types.size());
  fPrimaryTracks.reserve(track_histo_types.size());

  d_top->mkdir("ptrack");
  d_top->mkdir("alltrack");
  for (unsigned i = 0; i < track_histo_types.size(); i++) {
    fAllTracks.emplace_back();
    fPrimaryTracks.emplace_back();
    fPrimaryTrackProfiles.emplace_back();

    std::string dirname = track_histo_types[i];
    std::string all_dirname = "alltrack/" + dirname;
    
    d_top->mkdir(all_dirname.c_str());
    d_top->cd(all_dirname.c_str());
    fAllTracks[i].Initialize("", detector, max_length);
    d_top->cd();
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      std::string p_dirname = "ptrack/" + dirname + cut_order[j];
      d_top->mkdir(p_dirname.c_str());
      d_top->cd(p_dirname.c_str());
      fPrimaryTracks[i][j].Initialize("", detector, max_length);
      d_top->cd();

      for (unsigned k = 0; k < track_profile_types.size(); k++) {
        if (j == 0) fPrimaryTrackProfiles[i].emplace_back();
        unsigned n_bin;
        float xlo, xhi;
        std::tie(n_bin, xlo, xhi) = track_profile_xranges[k];
        std::string p_profile_dirname = "ptrack/" + dirname + cut_order[j] + "/profile_" + track_profile_types[k];
        d_top->mkdir(p_profile_dirname.c_str());
        d_top->cd(p_profile_dirname.c_str());
        fPrimaryTrackProfiles[i][k][j].Initialize("", n_bin, xlo, xhi);
        d_top->cd();
      } 
    }
  }


  for (unsigned i = 0; i < 4; i++) {
    Merge(fCosmic[i]);
  }
  for (unsigned i = 0; i < Histograms::nPIDs; i++) {
    for (unsigned j = 0; j < 5; j ++) {
      Merge(fParticles[i][j]);
    }
  }
  for (unsigned i = 0; i < Histograms::nHistos; i++) {
    for (const auto mode: Histograms::allModes) {
      Merge(fInteraction[i][mode]);
    }
  }

  for (unsigned i = 0; i < Cuts::nCuts; i++) {
    Merge(fCRTs[i]);
  }

  for (unsigned i = 0; i < fAllTracks.size(); i++) {
    Merge(fAllTracks[i]);
    for (unsigned j = 0; j < Cuts::nCuts; j++) {
      Merge(fPrimaryTracks[i][j]);
      for (unsigned k = 0; k < track_profile_types.size(); k++) {
        Merge(fPrimaryTrackProfiles[i][k][j]);
      } 
    } 
  }

}

  } // namespace SBNOsc
} // namespace ana
