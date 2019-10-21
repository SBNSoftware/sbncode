#include "Cuts.h"
#include "Derived.h"
#include "../RecoUtils/GeoUtil.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry) {
  fConfig.active_volumes.clear();
  fConfig.fiducial_volumes.clear();

  fConfig.trackMatchCompletionCut = cfg.get<double>("trackMatchCompletionCut", -1);
  fConfig.active_volumes = SBNRecoUtils::ActiveVolumes(geometry); 
  fConfig.TruthCompletion = cfg.get<float>("TruthCompletion", 0.5);
  fConfig.TruthMatchDist = cfg.get<float>("TruthMatchDist", 5.);

  fConfig.TrackLength = cfg.get<float>("TrackLength", 50.);

  fConfig.CRTHitDist = cfg.get<float>("CRTHitDist", 35.);
  fConfig.CRTHitTimeRange = cfg.get<std::array<float, 2>>("CRTHitTimeRange", {-0.4, 2.0});
  fConfig.CRTTrackAngle = cfg.get<float>("CRTTrackAngle", 0.4);

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

  {
    fhicl::ParameterSet dFV = \
     cfg.get<fhicl::ParameterSet>("containment_volume_inset");
    double dx = dFV.get<double>("x");
    double dy = dFV.get<double>("y");
    double zfront = dFV.get<double>("zfront");
    double zback = dFV.get<double>("zback");
    for (const geo::BoxBoundedGeo &geo: fConfig.active_volumes) {
      fConfig.containment_volumes.emplace_back(geo.MinX() + dx, geo.MaxX() - dx, geo.MinY() + dy, geo.MaxY() - dy, geo.MinZ() + zfront, geo.MaxZ() - zback);
    }
  }

}

std::array<bool, Cuts::nTruthCuts> Cuts::ProcessTruthCuts(const numu::RecoEvent &event, unsigned truth_vertex_index) const {
  bool is_neutrino = event.truth[truth_vertex_index].match.mode == numu::mCC || event.truth[truth_vertex_index].match.mode == numu::mNC;
  bool is_fiducial = InFV(event.truth[truth_vertex_index].position) && is_neutrino;
  bool is_matched = (fConfig.TruthMatchDist < 0. || dist2Match(event.truth[truth_vertex_index], event.reco) < fConfig.TruthMatchDist) 
                    && is_fiducial;
  bool is_completed = (fConfig.TruthCompletion < 0. || trackMatchCompletion(truth_vertex_index, event) > fConfig.TruthCompletion)
                      && is_matched;

  bool has_reco = false; 
  for (unsigned i = 0; i < event.reco.size(); i++) {
    if (event.reco[i].match.has_match && event.reco[i].match.event_vertex_id == truth_vertex_index) {
      has_reco = true;
      break;
    }
  }

  return {is_neutrino, is_fiducial, is_matched, is_completed, has_reco};
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, unsigned reco_vertex_index) const {
  bool is_reco = true;

  // require fiducial
  bool fiducial = InFV(event.reco[reco_vertex_index].position);

  const numu::RecoTrack &primary_track = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index);
  
  bool pass_crt_track = !HasCRTTrackMatch(primary_track) && fiducial;
  bool pass_crt_hit = (!HasCRTHitMatch(primary_track) || TimeInSpill(CRTMatchTime(primary_track)))
    && fiducial && pass_crt_track;

  bool pass_length = fConfig.TrackLength < 0. || primary_track.length > fConfig.TrackLength
    && pass_crt_hit;

  bool is_contained = primary_track.is_contained && fiducial && pass_crt_track && pass_crt_hit && pass_length;

  return {
    is_reco,
    fiducial,
    pass_crt_track,
    pass_crt_hit,
    pass_length,
    is_contained
  };
}

bool Cuts::HasCRTTrackMatch(const numu::RecoTrack &track) const {
  if (track.crt_match.size() == 0) return false;
  return track.crt_match[0].track.present && (fConfig.CRTTrackAngle < 0. || track.crt_match[0].track.angle < fConfig.CRTTrackAngle);
}

bool Cuts::HasCRTHitMatch(const numu::RecoTrack &track) const {
  if (track.crt_match.size() == 0) return false;
  return track.crt_match[0].hit.present && (fConfig.CRTHitDist < 0. || track.crt_match[0].hit.distance < fConfig.CRTHitDist);
}

float Cuts::CRTMatchTime(const numu::RecoTrack &track) const {
  if (HasCRTTrackMatch(track)) return track.crt_match[0].track.time;
  if (HasCRTHitMatch(track)) return track.crt_match[0].hit.time;
  return -99999;
}

bool Cuts::TimeInSpill(float time) const {
  return time > fConfig.CRTHitTimeRange[0] && time < fConfig.CRTHitTimeRange[1]; 
}

bool Cuts::InFV(const geo::Point_t &v) const {
  for (auto const& FV: fConfig.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

bool Cuts::InFV(const TVector3 &v) const {
  for (auto const& FV: fConfig.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

bool Cuts::InContainment(const TVector3 &v) const {
  for (auto const& FV: fConfig.containment_volumes) {
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
