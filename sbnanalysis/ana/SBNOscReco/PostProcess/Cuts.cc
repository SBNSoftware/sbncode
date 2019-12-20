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
  fConfig.CRTHitTimeRange = cfg.get<std::array<float, 2>>("CRTHitTimeRange", {-0.2, 2.0});
  fConfig.CRTActivityTimeRange = cfg.get<std::array<float, 2>>("CRTActivityTimeRange", {-0.2, 2.0});
  fConfig.CRTTrackAngle = cfg.get<float>("CRTTrackAngle", 0.4);

  fConfig.TruthFlashMatch = cfg.get<bool>("TruthFlashMatch", false);
  fConfig.FlashMatchScore = cfg.get<float>("FlashMatchScore", 10.);

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
  bool has_trigger = PassFlashTrigger(event) && is_fiducial;

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

  return {is_truth, is_fiducial, has_trigger, is_matched, is_completed, has_reco};
}

bool Cuts::PassFlashTrigger(const numu::RecoEvent &event) const {
  return numu::HasTrigger(event.flash_trigger_primitives, fConfig.PMTTriggerTreshold, fConfig.PMTNAboveThreshold);
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, 
						    unsigned reco_vertex_index, 
						    bool fSequentialCuts) const {
  bool is_reco = true;

  bool has_trigger = PassFlashTrigger(event);

  // require an in-time flash time
  bool has_intime_flash = false;
  for (const numu::RecoInteraction &reco: event.reco) {
    if (reco.slice.flash_match.present && reco.slice.flash_match.time > 0. /*intime*/) {
      has_intime_flash = true;
      break;
    }
  }

  // require fiducial
  bool fiducial = InFV(event.reco[reco_vertex_index].position) && has_trigger;

  const numu::RecoTrack &primary_track = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index);

  //   good_mcs = track contain OR (range momentum OR mcs trk length)
  bool good_mcs = ( ( InCalorimetricContainment(primary_track.start) && 
		      InCalorimetricContainment(primary_track.end) )  || 
		    ( primary_track.mcs_momentum < 7. /*garbage value*/&& 
		      (fConfig.MCSTrackLength <0. || 
		       primary_track.length > fConfig.MCSTrackLength) )
		    );

  // allow truth-based flash matching or reco based
  bool time_in_spill = false;
  if (primary_track.match.has_match) {
    time_in_spill = TimeInSpill(event.true_tracks.at(primary_track.match.mcparticle_id).start_time);
  }

  bool flashmatch;
  if (fConfig.TruthFlashMatch) flashmatch = time_in_spill;
  else flashmatch = fConfig.FlashMatchScore < 0. || event.reco[reco_vertex_index].slice.flash_match.score < fConfig.FlashMatchScore;

  bool pass_crt_track = !HasCRTTrackMatch(primary_track);

  bool pass_crt_hit = ( !HasCRTHitMatch(primary_track) || 
			TimeInSpill(CRTMatchTime(primary_track)) );

  bool pass_length = fConfig.TrackLength < 0. || 
                     primary_track.length > fConfig.TrackLength;

  bool is_contained = ( InCosmicContainment(primary_track.start) && 
			InCosmicContainment(primary_track.end) );


  bool no_crt_activity = true; //is_contained;
  for (const numu::CRTHit &crt_hit: event.in_time_crt_hits) {
     if ( TimeInSpill(crt_hit.time)  && 
	  crt_hit.pes > fConfig.CRTActivityPEThreshold ) 
       no_crt_activity = false;
     if ( !no_crt_activity )  break;
  }

  // Sequential Cuts
  bool has_trigger_seq      = is_reco              && has_trigger;
  bool has_intime_flash_seq = has_trigger_seq      && has_intime_flash;
  bool fiducial_seq         = has_intime_flash_seq && fiducial;
  bool good_mcs_seq         = fiducial_seq         && good_mcs;
  bool flashmatch_seq       = good_mcs_seq         && flashmatch;
  bool pass_crt_track_seq   = flashmatch_seq       && pass_crt_track;
  bool pass_crt_hit_seq     = pass_crt_track_seq   && pass_crt_hit;
  bool no_crt_activity_seq  = pass_crt_hit_seq     && no_crt_activity;
  bool is_contained_seq     = no_crt_activity_seq  && is_contained;
  bool pass_length_seq      = is_contained_seq     && pass_length;

  if ( fSequentialCuts ){
    has_trigger = has_trigger_seq;
    has_intime_flash = has_intime_flash_seq;
    fiducial = fiducial_seq;
    good_mcs = good_mcs_seq;
    flashmatch  = flashmatch_seq;
    pass_crt_track = pass_crt_track_seq;
    pass_crt_hit = pass_crt_hit_seq;
    no_crt_activity = no_crt_activity_seq;
    is_contained = is_contained_seq;
    pass_length = pass_length_seq;
  }

  return {
    is_reco,
    has_trigger,
    has_intime_flash,
    fiducial,
    good_mcs,
    flashmatch,
    pass_crt_track,
    pass_crt_hit,
    no_crt_activity,
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
