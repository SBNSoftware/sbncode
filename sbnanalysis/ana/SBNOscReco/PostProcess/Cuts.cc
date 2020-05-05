#include "Cuts.h"
#include "../Histograms/Derived.h"
#include "../RecoUtils/GeoUtil.h"
#include "../TriggerEmulator/PMTTrigger.h"
#include "../NumuReco/TrackAlgo.h"

namespace ana {
 namespace SBNOsc {

void Cuts::Initialize(const fhicl::ParameterSet &cfg, const geo::GeometryCore *geometry) {
  fConfig.active_volumes.clear();
  fConfig.fiducial_volumes.clear();

  fConfig.UseTrueVertex = cfg.get<bool>("UseTrueVertex", false);
  fConfig.RequireTrueVertex = cfg.get<bool>("RequireTrueVertex", false);
  fConfig.trackMatchCompletionCut = cfg.get<double>("trackMatchCompletionCut", -1);
  fConfig.active_volumes = SBNRecoUtils::ActiveVolumes(geometry); 
  fConfig.TruthCompletion = cfg.get<float>("TruthCompletion", 0.5);
  fConfig.TruthMatchDist = cfg.get<float>("TruthMatchDist", 5.);

  fConfig.MCSTrackLength = cfg.get<float>("MCSTrackLength", 100.);
  fConfig.TrackLength = cfg.get<float>("TrackLength", 50.);

  fConfig.CRTHitDist = cfg.get<float>("CRTHitDist", 35.);
  fConfig.CRTHitTimeRange = cfg.get<std::array<float, 2>>("CRTHitTimeRange", {-0.2, 2.0});
  fConfig.CRTActivityTimeRange = cfg.get<std::array<float, 2>>("CRTActivityTimeRange", {-0.2, 2.0});
  fConfig.CRTTrackAngle = cfg.get<float>("CRTTrackAngle", -1.);

  fConfig.TruthFlashMatch = cfg.get<bool>("TruthFlashMatch", false);
  fConfig.FlashMatchScore = cfg.get<float>("FlashMatchScore", 10.);

  fConfig.PMTTriggerTreshold = cfg.get<int>("PMTTriggerTreshold", 7950);
  fConfig.PMTNAboveThreshold = cfg.get<unsigned>("PMTNAboveThreshold", 5);

  fConfig.CRTActivityPEThreshold = cfg.get<float>("CRTActivityPEThreshold", 100.);

  fConfig.CutOrder = cfg.get<std::vector<std::string>>("CutOrder");
  fConfig.TruthCutOrder = cfg.get<std::vector<std::string>>("TruthCutOrder");

  fConfig.SelectMaxNuScore = cfg.get<bool>("SelectMaxNuScore", false);
  fConfig.NuScoreCut = cfg.get<float>("NuScoreCut", -2.);

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

  if (cfg.has_key("XGBoostPIDModel")) {
    std::string model_file = cfg.get<std::string>("XGBoostPIDModel");
    fXGBPID.SetModelFile(model_file.c_str());
  }
  fConfig.ProtonMuonScore = cfg.get<float>("ProtonMuonScore", 0.5);

  fConfig.RequireIdentCCMuon = cfg.get<bool>("RequireIdentCCMuon", true);
  fConfig.RequireRecoCCMuon = cfg.get<bool>("RequireRecoCCMuon", true);

}

std::array<bool, Cuts::nTruthCuts> Cuts::ProcessTruthCuts(const numu::RecoEvent &event, const event::Event &core, unsigned truth_vertex_index, bool SequentialCuts) const {
  std::map<std::string, bool> cuts;
  cuts["Truth"] = true;


  cuts["T_fid"]  = InFV(core.truth[truth_vertex_index].neutrino.position);

  cuts["T_trig"] = PassFlashTrigger(event);

  cuts["T_vqual"]  = (fConfig.TruthMatchDist < 0. || 
                     dist2Match(core.truth[truth_vertex_index], event.reco) < fConfig.TruthMatchDist); 
  cuts["T_tqual"]  = (fConfig.TruthCompletion < 0. || 
		       trackMatchCompletion(truth_vertex_index, event) > fConfig.TruthCompletion);

  bool has_reco = false;
  for (unsigned i = 0; i < event.reco.size(); i++) {
    const numu::RecoTrack &primary_track = event.tracks.at(event.reco[i].primary_track_index);
    if (event.reco[i].slice.truth.interaction_id == truth_vertex_index && 
        event.reco[i].slice.truth.IsPrimary()) {
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

bool Cuts::IdentCCMuon(const numu::RecoEvent &event, const event::Event &core, unsigned reco_vertex_index) const {
  // check if we match to a CC interaction
  bool iscc = event.reco[reco_vertex_index].slice.truth.interaction_id >= 0 &&
         core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].neutrino.iscc;
  if (!iscc) return false;

  // now grab the muon
  int muon_id = core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].lepton.G4ID;
  if (!event.particles.count(muon_id)) return false;
  const numu::TrueParticle &muon = event.particles.at(muon_id);
  for (auto const &track_pair: event.tracks) {
    const numu::RecoTrack &track = track_pair.second;
    if (track.truth.matches.size() > 0 && track.truth.matches[0].G4ID == muon_id) {
      if (track.truth.matches[0].energy / muon.deposited_energy > 0.5) {
        return true;
      }
    } 
  }
  return false;
}

bool Cuts::RecoCCMuon(const numu::RecoEvent &event, const event::Event &core, unsigned reco_vertex_index) const {
  // check if we match to a CC interaction
  bool iscc = event.reco[reco_vertex_index].slice.truth.interaction_id >= 0 &&
         core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].neutrino.iscc;
  if (!iscc) return false;

  // now grab the muon
  int muon_id = core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].lepton.G4ID;
  if (!event.particles.count(muon_id)) return false;
  const numu::TrueParticle &muon = event.particles.at(muon_id);
  for (auto const &particle_pair: event.reco[reco_vertex_index].slice.particles) {
    if (!event.tracks.count(particle_pair.second.ID)) continue;

    const numu::RecoTrack &track = event.tracks.at(particle_pair.second.ID);
    if (track.truth.matches.size() > 0 && track.truth.matches[0].G4ID == muon_id) {
      if (track.truth.matches[0].energy / muon.deposited_energy > 0.5) {
        return true;
      }
    } 
  }
  return false;
}

std::array<bool, Cuts::nCuts> Cuts::ProcessRecoCuts(const numu::RecoEvent &event, 
                                                    const event::Event &core,
						    unsigned reco_vertex_index, 
						    bool fSequentialCuts) const {
  std::map<std::string, bool> cuts;
  cuts["Reco"] = true;

  cuts["R_goodreco"] = (!fConfig.RequireIdentCCMuon || IdentCCMuon(event, core, reco_vertex_index)) &&
                       (!fConfig.RequireRecoCCMuon  || RecoCCMuon(event, core, reco_vertex_index));

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
    cuts["R_fid"] = event.reco[reco_vertex_index].slice.truth.interaction_id >= 0 &&
                      InFV(core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].neutrino.position);
  }
  else {
    cuts["R_fid"] = InFV(event.reco[reco_vertex_index].position);
  }

  if (fConfig.RequireTrueVertex) {
    cuts["R_fid"] = cuts["R_fid"] && 
	    event.reco[reco_vertex_index].slice.truth.interaction_id >= 0 &&
	    InFV(core.truth[event.reco[reco_vertex_index].slice.truth.interaction_id].neutrino.position);
  }

  const numu::RecoTrack &primary_track = event.tracks.at(event.reco[reco_vertex_index].primary_track_index);

  //   good_mcs = track contain OR (range momentum OR mcs trk length)
  cuts["R_goodmcs"] = ( ( InCalorimetricContainment(primary_track.start) && 
		      InCalorimetricContainment(primary_track.end) )  || 
		    ( numu::MCSMomentum(primary_track) < 7. /*garbage value*/&& 
		      (fConfig.MCSTrackLength <0. || 
		       primary_track.length > fConfig.MCSTrackLength) )
		    );

  // allow truth-based flash matching or reco based
  bool time_in_spill = false;
  int ptrack_id = primary_track.truth.GetPrimaryMatchID();
  if (event.particles.count(ptrack_id)) {
    time_in_spill = TimeInSpill(event.particles.at(ptrack_id).start_time);
  }

  if (fConfig.TruthFlashMatch) cuts["R_flashmatch"] = time_in_spill;
  else cuts["R_flashmatch"] = fConfig.FlashMatchScore < 0. || (event.reco[reco_vertex_index].slice.flash_match.present && event.reco[reco_vertex_index].slice.flash_match.score < fConfig.FlashMatchScore);

  cuts["R_crttrack"] = !HasCRTTrackMatch(primary_track);

  cuts["R_crthit"] = ( !HasCRTHitMatch(primary_track) || 
			TimeInSpill(CRTMatchTime(primary_track)) );

  cuts["R_length"]  = fConfig.TrackLength < 0. || 
                     primary_track.length > fConfig.TrackLength;

  float nu_score = event.reco[reco_vertex_index].slice.particles.at(event.reco[reco_vertex_index].slice.primary_index).p_nu_score;
  bool is_max_nu_score = true;
  for (unsigned i = 0; i < event.reco.size(); i++) {
    float this_nu_score = event.reco[i].slice.particles.at(event.reco[i].slice.primary_index).p_nu_score;
    if (i != reco_vertex_index && this_nu_score > nu_score) {
      is_max_nu_score = false;
      break;
    } 
  } 

  cuts["R_maxnuscore"] = !fConfig.SelectMaxNuScore || is_max_nu_score;

  cuts["R_nuscore"] = nu_score > fConfig.NuScoreCut;

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
         fConfig.CRTTrackAngle >= 0. &&
	  track.crt_match.track.angle < fConfig.CRTTrackAngle;
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

float Cuts::PredictTrack(const numu::RecoTrack &track) const {
  if (!fXGBPID.Ready()) return 1.;

  float chi2_proton = (track.chi2_proton > 0) ? track.chi2_proton : NAN;
  float chi2_muon = (track.chi2_muon > 0) ? track.chi2_muon : NAN;
  float chi2_diff = (track.chi2_muon > 0 && track.chi2_proton > 0) ? track.chi2_proton - track.chi2_muon : NAN;
  float crt_hit_distance = (track.crt_match.hit_match.present) ? track.crt_match.hit_match.distance : NAN;

  std::map<std::string, float> data;

  data["theta"] = track.theta;
  data["phi"] = track.phi;
  data["length"] = track.length;
  data["contained"] = track.is_contained;
  data["chi2_proton"] = chi2_proton;
  data["chi2_muon"] = chi2_muon;
  data["chi2_diff"] = chi2_diff;
  data["crt_hit_distance"] = crt_hit_distance;

  return fXGBPID.PredictOne(data);
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
