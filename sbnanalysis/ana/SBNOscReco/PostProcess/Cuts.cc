#include "Cuts.h"
#include "../Histograms/Derived.h"
#include "../RecoUtils/GeoUtil.h"
#include "../TriggerEmulator/PMTTrigger.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry) {
  fConfig.active_volumes.clear();
  fConfig.fiducial_volumes.clear();

  fConfig.trackMatchCompletionCut = cfg.get<double>("trackMatchCompletionCut", -1);
  fConfig.active_volumes = SBNRecoUtils::ActiveVolumes(geometry); 
  fConfig.TruthCompletion = cfg.get<float>("TruthCompletion", 0.5);
  fConfig.TruthMatchDist = cfg.get<float>("TruthMatchDist", 5.);

  fConfig.MCSTrackLength = cfg.get<float>("MCSTrackLength", 100.);
  fConfig.TrackLength = cfg.get<float>("TrackLength", 50.);

  fConfig.CRTHitDist = cfg.get<float>("CRTHitDist", 35.);
  fConfig.CRTHitTimeRange = cfg.get<std::array<float, 2>>("CRTHitTimeRange", {-0.2, 1.8});
  fConfig.CRTActivityTimeRange = cfg.get<std::array<float, 2>>("CRTActivityTimeRange", {-0.2, 1.8});
  fConfig.CRTTrackAngle = cfg.get<float>("CRTTrackAngle", 0.4);

  fConfig.TruthFlashMatch = cfg.get<bool>("TruthFlashMatch", true);

  fConfig.PMTTriggerTreshold = cfg.get<int>("PMTTriggerTreshold", 7950);
  fConfig.PMTNAboveThreshold = cfg.get<unsigned>("PMTNAboveThreshold", 5);

  fConfig.CRTActivityPEThreshold = cfg.get<float>("CRTActivityPEThreshold", 100.);

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
     cfg.get<fhicl::ParameterSet>("calorimetric_containment_volume_inset");
    double dx = dFV.get<double>("x");
    double dy = dFV.get<double>("y");
    double zfront = dFV.get<double>("zfront");
    double zback = dFV.get<double>("zback");
    for (const geo::BoxBoundedGeo &geo: fConfig.active_volumes) {
      fConfig.calorimetric_containment_volumes.emplace_back(geo.MinX() + dx, 
							    geo.MaxX() - dx, 
							    geo.MinY() + dy, 
							    geo.MaxY() - dy, 
							    geo.MinZ() + 
							    zfront, 
							    geo.MaxZ() -zback);
    }
  }

  {
    fhicl::ParameterSet dFV = \
     cfg.get<fhicl::ParameterSet>("cosmic_containment_volume_inset");
    double ytop = dFV.get<double>("ytop");
    double ybottom = dFV.get<double>("ybottom");
    double zfront = dFV.get<double>("zfront");
    double zback = dFV.get<double>("zback");
    for (const geo::BoxBoundedGeo &geo: fConfig.active_volumes) {
      VolYZ contain;
      contain.Y = {geo.MinY() + ybottom, geo.MaxY() - ytop};
      contain.Z = {geo.MinZ() + zfront, geo.MaxZ() - zback};
      fConfig.cosmic_containment_volumes.emplace_back(contain);
    }
  }

}

std::array<bool, Cuts::nTruthCuts> Cuts::ProcessTruthCuts(const numu::RecoEvent &event, unsigned truth_vertex_index) const {
  bool is_truth = true;

  bool is_neutrino = event.truth[truth_vertex_index].match.mode == numu::mCC ||
                     event.truth[truth_vertex_index].match.mode == numu::mNC;
  bool is_fiducial = InFV(event.truth[truth_vertex_index].position) &&
                     is_neutrino;
  bool is_matched = (fConfig.TruthMatchDist < 0. || 
                     dist2Match(event.truth[truth_vertex_index], event.reco) < fConfig.TruthMatchDist) 
                     && is_fiducial;
  bool is_completed = (fConfig.TruthCompletion < 0. || 
		       trackMatchCompletion(truth_vertex_index, event) > fConfig.TruthCompletion)
                      && is_matched;

  bool has_reco = false;
 
  for (unsigned i = 0; i < event.reco.size(); i++) {
    if (event.reco[i].match.has_match && 
	event.reco[i].match.event_track_id == truth_vertex_index && 
	event.reco[i].primary_track.match.is_primary) {
      has_reco = true;
      break;
    }
  }

  return {is_truth, is_fiducial, is_matched, is_completed, has_reco};
}

bool Cuts::PassFlashTrigger(const numu::RecoEvent &event) const {
  return numu::HasTrigger(event.flash_trigger_primitives, fConfig.PMTTriggerTreshold, fConfig.PMTNAboveThreshold);
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, 
						    unsigned reco_vertex_index) const {
  bool is_reco = true;

  bool has_trigger = PassFlashTrigger(event);

  // require fiducial
  bool fiducial = InFV(event.reco[reco_vertex_index].position) && has_trigger;

  const numu::RecoTrack &primary_track = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index);

  bool good_mcs = (InCalorimetricContainment(primary_track.start) && 
		   InCalorimetricContainment(primary_track.end) )  || /* use range momentum */
    ( primary_track.mcs_momentum < 7. /* garbage value */ && 
      (fConfig.MCSTrackLength <0. || 
       primary_track.length > fConfig.MCSTrackLength) ) /* use MCS*/&& 
    fiducial;
  
  bool pass_crt_track = !HasCRTTrackMatch(primary_track) && good_mcs;
  bool pass_crt_hit = ( !HasCRTHitMatch(primary_track) || 
			TimeInSpill(CRTMatchTime(primary_track) )) && 
                        pass_crt_track;

  // for now, use truth-based flash matching
  bool time_in_spill = false;
  if (primary_track.match.has_match) {
    time_in_spill = TimeInSpill(event.true_tracks.at(primary_track.match.mcparticle_id).start_time);
  }
  bool flash_match = (!fConfig.TruthFlashMatch || time_in_spill) && pass_crt_hit; 

  bool crt_activity = flash_match;
  for (const numu::CRTHit &crt_hit: event.in_time_crt_hits) {
     if (TimeInCRTActiveSpill(crt_hit.time) && crt_hit.pes > fConfig.CRTActivityPEThreshold) {
       crt_activity = false;
     }
     if (!crt_activity) break;
  }

  bool is_contained = InCosmicContainment(primary_track.start) && 
                      InCosmicContainment(primary_track.end) && 
                      crt_activity;

  bool pass_length = fConfig.TrackLength < 0. || 
                     primary_track.length > fConfig.TrackLength
                     && is_contained;
  return {
    is_reco,
    has_trigger,
    fiducial,
    good_mcs,
    pass_crt_track,
    pass_crt_hit,
    flash_match,
    crt_activity,
    is_contained,
    pass_length
  };
}

bool Cuts::HasCRTTrackMatch(const numu::RecoTrack &track) const {
  return track.crt_match.track.present && 
         (fConfig.CRTTrackAngle < 0. || 
	  track.crt_match.track.angle < fConfig.CRTTrackAngle);
}

bool Cuts::HasCRTHitMatch(const numu::RecoTrack &track) const {
  return track.crt_match.hit_match.present && 
         (fConfig.CRTHitDist < 0. || 
	  track.crt_match.hit_match.distance < fConfig.CRTHitDist);
}

float Cuts::CRTMatchTime(const numu::RecoTrack &track) const {
  if (HasCRTTrackMatch(track)) return track.crt_match.track.time;
  if (HasCRTHitMatch(track)) return track.crt_match.hit_match.time;
  return -99999;
}

bool Cuts::TimeInSpill(float time) const {
  return time > fConfig.CRTHitTimeRange[0] && 
         time < fConfig.CRTHitTimeRange[1]; 
}

bool Cuts::TimeInCRTActiveSpill(float time) const {
  return time > fConfig.CRTActivityTimeRange[0] && 
         time < fConfig.CRTActivityTimeRange[1]; 
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

bool Cuts::InCalorimetricContainment(const TVector3 &v) const {
  for (auto const& FV: fConfig.calorimetric_containment_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

bool Cuts::InCosmicContainment(const TVector3 &v) const {
  for (auto const& CV: fConfig.cosmic_containment_volumes) {
    if (CV.Z[0] < v.Z() 
     && CV.Z[1] > v.Z() 
     && CV.Y[0] < v.Y()
     && CV.Y[1] > v.Y()) {
      return true;
    }
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
