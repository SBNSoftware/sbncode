//////////////////////////////////////////////////////////////////////
// \file    FillReco.cxx
// \brief   Fill reco SR branches 
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillReco.h"
#include "RecoUtils/RecoUtils.h"

// declare helpers
caf::SRTrackTruth MatchTrack2Truth(const std::vector<art::Ptr<recob::Hit>> &hits);

namespace caf
{

  //......................................................................
  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic) {
    return slice.is_clear_cosmic || !cut_clear_cosmic;
  }

  //......................................................................
  void FillSliceVars(const recob::Slice& slice,
                     const SliceData sdata,
                     caf::SRSlice& srslice,
                     bool allowEmpty)
  {

    srslice.id           = slice.ID() + sdata.slice_id_offset;
    srslice.charge       = slice.Charge();

    std::cout << "Slice: " << srslice.id << " has primary particle: " << (sdata.primary != NULL) << std::endl;

    assert(sdata.primary != NULL);
    assert(sdata.primary_meta != NULL);
    
    // collect the properties
    auto const &properties = sdata.primary_meta->GetPropertiesMap();
    if (properties.count("IsClearCosmic")) {
      assert(!properties.count("IsNeutrino"));
      srslice.is_clear_cosmic = false;
    }
    else {
      assert(properties.count("IsNeutrino"));
      srslice.is_clear_cosmic = true;
    }
    if (properties.count("NuScore")) {
      srslice.nu_score = properties.at("NuScore"); 
    } 
    else {
      srslice.nu_score = -1;
    }

    if (sdata.fmatch != NULL) {
      srslice.fmatch.present = true;
      srslice.fmatch.time = sdata.fmatch->Time();
      srslice.fmatch.score = sdata.fmatch->TriggerConfidence();
      srslice.fmatch.pe = sdata.fmatch->TriggerType();
    }
    else {
      srslice.fmatch.present = false;
    }

    // get the priamry tracks/showers
    for (unsigned id: sdata.primary->Daughters()) {
      srslice.primary.push_back(id + sdata.particle_id_offset);
    }

  }
  //......................................................................

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle &particle,
                     const caf::TrackData tdata,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {

    std::cout << "Track with ID: " << track.ID() << " len: " << track.Length() << std::endl;

    srtrack.npts = track.CountValidPoints();
    srtrack.len  = track.Length();
    srtrack.costh = track.StartDirection().Z() / sqrt(track.StartDirection().Mag2());

    srtrack.ID = particle.Self() + tdata.particle_index_offset;

    // set the daughters in the particle flow
    for (unsigned id: particle.Daughters()) {
      srtrack.daughters.push_back(id + tdata.particle_index_offset);
    }

    // calculate MCS fits
    recob::MCSFitResult mcs_fit_muon= tdata.mcs_calculator->fitMcs(track, 13);
    srtrack.mcs_muon.fwd_momentum = mcs_fit_muon.fwdMomentum();
    srtrack.mcs_muon.fwd_momentum_err = mcs_fit_muon.fwdMomUncertainty();
    srtrack.mcs_muon.bwd_momentum = mcs_fit_muon.bwdMomentum();
    srtrack.mcs_muon.bwd_momentum_err = mcs_fit_muon.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_proton= tdata.mcs_calculator->fitMcs(track, 2212);
    srtrack.mcs_proton.fwd_momentum = mcs_fit_proton.fwdMomentum();
    srtrack.mcs_proton.fwd_momentum_err = mcs_fit_proton.fwdMomUncertainty();
    srtrack.mcs_proton.bwd_momentum = mcs_fit_proton.bwdMomentum();
    srtrack.mcs_proton.bwd_momentum_err = mcs_fit_proton.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_pion= tdata.mcs_calculator->fitMcs(track, 211);
    srtrack.mcs_pion.fwd_momentum = mcs_fit_pion.fwdMomentum();
    srtrack.mcs_pion.fwd_momentum_err = mcs_fit_pion.fwdMomUncertainty();
    srtrack.mcs_pion.bwd_momentum = mcs_fit_pion.bwdMomentum();
    srtrack.mcs_pion.bwd_momentum_err = mcs_fit_pion.bwdMomUncertainty();

    recob::MCSFitResult mcs_fit_kaon= tdata.mcs_calculator->fitMcs(track, 321);
    srtrack.mcs_kaon.fwd_momentum = mcs_fit_kaon.fwdMomentum();
    srtrack.mcs_kaon.fwd_momentum_err = mcs_fit_kaon.fwdMomUncertainty();
    srtrack.mcs_kaon.bwd_momentum = mcs_fit_kaon.bwdMomentum();
    srtrack.mcs_kaon.bwd_momentum_err = mcs_fit_kaon.bwdMomUncertainty();

    // calculate range momentum
    srtrack.range_momentum_muon = tdata.range_calculator->GetTrackMomentum(track.Length(), 13);
    srtrack.range_momentum_proton = tdata.range_calculator->GetTrackMomentum(track.Length(), 2212);

    // get the particle ID's
    //
    // iterate over the planes -- use the conduction plane to get the particle ID
    srtrack.pid_ndof = 0;
    for (unsigned i = 0; i < 3; i++) { 
      const anab::ParticleID &particle_id = *tdata.particleIDs[i];
      if (particle_id.PlaneID() && tdata.geom->SignalType(particle_id.PlaneID()) == geo::kCollection) {
        srtrack.chi2_muon = particle_id.Chi2Muon();
        srtrack.chi2_pion = particle_id.Chi2Kaon();
        srtrack.chi2_kaon = particle_id.Chi2Pion();
        srtrack.chi2_proton = particle_id.Chi2Proton();
        srtrack.pid_ndof = particle_id.Ndf();
      }
    }

    // bad particle ID -- set chi2 to -1
    if (srtrack.pid_ndof == 0) {
      srtrack.chi2_muon = -1;
      srtrack.chi2_proton = -1;
      srtrack.chi2_kaon = -1;
      srtrack.chi2_pion = -1;
    }

    // TODO: what to do with calorimetry
    
    // TODO: crt matching

    // Truth matching
    srtrack.truth = MatchTrack2Truth(tdata.hits);    

    std::cout << "Track matched to particle: " << srtrack.truth.GetPrimaryMatchID() << std::endl;
  }
  //......................................................................
} // end namespace 

// define helpers
caf::SRTrackTruth MatchTrack2Truth(const std::vector<art::Ptr<recob::Hit>> &hits) {

  // this id is the same as the mcparticle ID as long as we got it from geant4
  std::vector<std::pair<int, float>> matches = CAFRecoUtils::AllTrueParticleIDEnergyMatches(hits, true);
  float total_energy = CAFRecoUtils::TotalHitEnergy(hits);

  caf::SRTrackTruth ret;

  ret.total_deposited_energy = total_energy / 1000. /* MeV -> GeV */;

  // setup the matches
  for (auto const &pair: matches) {
    caf::SRTrackTruth::ParticleMatch match;
    match.G4ID = pair.first;
    match.energy = pair.second / 1000. /* MeV -> GeV */;
    ret.matches.push_back(match);
  }

  // sort highest energy match to lowest
  std::sort(ret.matches.begin(), ret.matches.end(), 
    [](const caf::SRTrackTruth::ParticleMatch &a, const caf::SRTrackTruth::ParticleMatch &b) {
      return a.energy > b.energy;
    }
  );

  return ret;
}






