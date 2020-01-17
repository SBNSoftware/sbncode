#include "Cuts.h"
#include "../Histograms/Derived.h"
#include "../RecoUtils/GeoUtil.h"
#include "../TriggerEmulator/PMTTrigger.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry) {
  fConfig.active_volumes.clear();
  fConfig.fiducial_volumes.clear();

  fConfig.UseTrueVertex = cfg.get<bool>("UseTrueVertex", true);
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

  fConfig.CutOrder = cfg.get<std::vector<std::string>>("CutOrder");
  fConfig.TruthCutOrder = cfg.get<std::vector<std::string>>("TruthCutOrder");

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

std::array<bool, Cuts::nTruthCuts> Cuts::ProcessTruthCuts(const numu::RecoEvent &event, unsigned truth_vertex_index, bool SequentialCuts) const {
  std::map<std::string, bool> cuts;
  cuts["Truth"] = true;

  bool is_neutrino = event.truth[truth_vertex_index].match.mode == numu::mCC ||
                     event.truth[truth_vertex_index].match.mode == numu::mNC;
  cuts["T_fid"]  = InFV(event.truth[truth_vertex_index].position) &&
                     is_neutrino;
  cuts["T_trig"] = PassFlashTrigger(event);

  cuts["T_vqual"]  = (fConfig.TruthMatchDist < 0. || 
                     dist2Match(event.truth[truth_vertex_index], event.reco) < fConfig.TruthMatchDist); 
  cuts["T_tqual"]  = (fConfig.TruthCompletion < 0. || 
		       trackMatchCompletion(truth_vertex_index, event) > fConfig.TruthCompletion);

  bool has_reco = false;
  for (unsigned i = 0; i < event.reco.size(); i++) {
    if (event.reco[i].match.has_match && 
	event.reco[i].match.event_track_id == truth_vertex_index && 
	event.reco[i].primary_track.match.is_primary) {
      has_reco = true;
      break;
    }
  }
  cuts["T_reco"] = has_reco;

  std::array<bool, Cuts::nTruthCuts> ret;
  for (unsigned i = 0; i < fConfig.TruthCutOrder.size(); i++) {
    if (SequentialCuts && i > 0) ret[i] = ret[i-1] && cuts.at(fConfig.TruthCutOrder[i]);
    else ret[i] = cuts.at(fConfig.TruthCutOrder[i]);
  }
  return ret;

}

bool Cuts::PassFlashTrigger(const numu::RecoEvent &event) const {
  return numu::HasTrigger(event.flash_trigger_primitives, fConfig.PMTTriggerTreshold, fConfig.PMTNAboveThreshold);
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, 
						    unsigned reco_vertex_index, 
						    bool fSequentialCuts) const {
  std::map<std::string, bool> cuts;
  cuts["Reco"] = true;

  cuts["R_trig"] = PassFlashTrigger(event);

  // require an in-time flash time
  bool has_intime_flash = false;
  for (const numu::RecoInteraction &reco: event.reco) {
    if (reco.slice.flash_match.present && reco.slice.flash_match.time > 0. /*intime*/) {
      has_intime_flash = true;
      break;
    }
  }
  cuts["R_flashtime"] = has_intime_flash;

  // require fiducial
  if (fConfig.UseTrueVertex) {
    cuts["R_fid"] = event.reco[reco_vertex_index].match.event_track_id >= 0 && InFV(event.truth.at(event.reco[reco_vertex_index].match.event_track_id).position);
  }
  else {
    cuts["R_fid"] = InFV(event.reco[reco_vertex_index].position);
  }

  const numu::RecoTrack &primary_track = event.reco_tracks.at(event.reco[reco_vertex_index].slice.primary_track_index);

  //   good_mcs = track contain OR (range momentum OR mcs trk length)
  cuts["R_goodmcs"] = ( ( InCalorimetricContainment(primary_track.start) && 
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

  if (fConfig.TruthFlashMatch) cuts["R_flashmatch"] = time_in_spill;
  else cuts["R_flashmatch"] = fConfig.FlashMatchScore < 0. || (event.reco[reco_vertex_index].slice.flash_match.present && event.reco[reco_vertex_index].slice.flash_match.score < fConfig.FlashMatchScore);

  cuts["R_crttrack"] = !HasCRTTrackMatch(primary_track);

  cuts["R_crthit"] = ( !HasCRTHitMatch(primary_track) || 
			TimeInSpill(CRTMatchTime(primary_track)) );

  cuts["R_length"]  = fConfig.TrackLength < 0. || 
                     primary_track.length > fConfig.TrackLength;

  cuts["R_contained"] = ( InCosmicContainment(primary_track.start) && 
			InCosmicContainment(primary_track.end) );


  cuts["R_crtactive"]  = true; //is_contained;
  for (const numu::CRTHit &crt_hit: event.in_time_crt_hits) {
     if ( TimeInSpill(crt_hit.time)  && 
	  crt_hit.pes > fConfig.CRTActivityPEThreshold )  {
       cuts["R_crtactive"]  = false;
       break;
     }
  }

  std::array<bool, Cuts::nCuts> ret;
  for (unsigned i = 0; i < fConfig.CutOrder.size(); i++) {
    if (fSequentialCuts && i > 0) ret[i] = ret[i-1] && cuts.at(fConfig.CutOrder[i]);
    else ret[i] = cuts.at(fConfig.CutOrder[i]);
  }
  return ret;
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
