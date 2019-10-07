#include "Cuts.h"
#include "../RecoUtils/GeoUtil.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry) {
  fConfig.active_volumes.clear();
  fConfig.fiducial_volumes.clear();

  fConfig.trackMatchCompletionCut = cfg.get<double>("trackMatchCompletionCut", -1);
  fConfig.active_volumes = SBNRecoUtils::ActiveVolumes(geometry); 

  {
    fhicl::ParameterSet dFV = \
     cfg.get<fhicl::ParameterSet>("fiducial_volume_inset");
    double dx = dFV.get<double>("x");
    double dy = dFV.get<double>("y");
    double zfront = dFV.get<double>("zfront");
    double zback = dFV.get<double>("zback");
    for (const geo::BoxBoundedGeo &geo: fConfig.active_volumes) {
      fConfig.fiducial_volumes.emplace_back(geo.MinX() + dx, geo.MaxX() - dx, geo.MinY() + dy, geo.MaxY() - dy, geo.MinZ() + zfront, geo.MaxZ() - zback);
    }
  }

}

std::array<bool, 2> Cuts::ProcessTruthCuts(const numu::RecoEvent &event, unsigned truth_vertex_index) const {
  bool is_neutrino = event.truth[truth_vertex_index].match.mode == numu::mCC || event.truth[truth_vertex_index].match.mode == numu::mNC;
  return {is_neutrino, is_neutrino &&  containedInFV(event.truth[truth_vertex_index].position)};
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, unsigned reco_vertex_index) const {
  bool is_reco = true;

  // require fiducial
  bool fiducial = containedInFV(event.reco[reco_vertex_index].position);

  // require close to truth
  bool v_quality = event.reco[reco_vertex_index].match.event_vertex_id >= 0 && event.reco[reco_vertex_index].match.truth_vertex_distance < 10.;

  bool t_quality = false;
  int t_mcparticle_id = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).match.mcparticle_id;
  for (const numu::RecoInteraction &truth: event.truth) {
    if (truth.slice.primary_track_index >= 0 && 
      event.reco_tracks.at(truth.slice.primary_track_index).match.mcparticle_id == t_mcparticle_id) {
      t_quality = true;
      break;
    } 
  }
  // require completion
  t_quality = t_quality && (fConfig.trackMatchCompletionCut < 0 ||
    event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).match.completion > fConfig.trackMatchCompletionCut);

  bool is_contained = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index).is_contained;

  return {
    is_reco,
    fiducial,
    v_quality,
    t_quality && v_quality,
    is_contained
  };
}

bool Cuts::containedInFV(const geo::Point_t &v) const {
  for (auto const& FV: fConfig.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

bool Cuts::containedInFV(const TVector3 &v) const {
  for (auto const& FV: fConfig.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}


/*
bool Cuts::SelectReco(std::array<bool, Cuts::nCuts> &cuts) {
  return 
    cuts[0] && 
    (cuts[1] || !fConfig.requireTrack) &&
    (cuts[3] || !fConfig.requireMatched) &&
    (cuts[4] || !fConfig.requireContained);
}
*/
  }
}
