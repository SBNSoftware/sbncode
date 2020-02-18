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
    return (slice.is_clear_cosmic || !cut_clear_cosmic) // No clear cosmics
           && slice.primary.size() > 0; // must have primary tracks/showers
  }

  //......................................................................
  void FillSliceVars(const recob::Slice& slice,
                     const SliceData sdata,
                     caf::SRSlice& srslice,
                     bool allowEmpty)
  {

    srslice.charge       = slice.Charge();

    assert(sdata.primary != NULL);
    assert(sdata.primary_meta != NULL);
 
    // default values   
    srslice.nu_score = -1;
    srslice.is_clear_cosmic = true;
 
    // collect the properties
    if (sdata.primary != NULL) {
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
    if (sdata.primary != NULL) {
      for (unsigned id: sdata.primary->Daughters()) {
        srslice.primary.push_back(id);
      }
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

    srtrack.ID = particle.Self();

    // set the daughters in the particle flow
    for (unsigned id: particle.Daughters()) {
      srtrack.daughters.push_back(id);
    }

    // calculate MCS fits
    recob::MCSFitResult mcs_fit_muon= tdata.mcs_calculator->fitMcs(track, 13);
    recob::MCSFitResult mcs_fit_proton= tdata.mcs_calculator->fitMcs(track, 2212);
    recob::MCSFitResult mcs_fit_pion= tdata.mcs_calculator->fitMcs(track, 211);
    recob::MCSFitResult mcs_fit_kaon= tdata.mcs_calculator->fitMcs(track, 321);

    srtrack.mcsP.fwdP_muon     = mcs_fit_muon.fwdMomentum();
    srtrack.mcsP.fwdP_err_muon = mcs_fit_muon.fwdMomUncertainty();
    srtrack.mcsP.bwdP_muon     = mcs_fit_muon.bwdMomentum();
    srtrack.mcsP.bwdP_err_muon = mcs_fit_muon.bwdMomUncertainty();

    srtrack.mcsP.fwdP_proton     = mcs_fit_proton.fwdMomentum();
    srtrack.mcsP.fwdP_err_proton = mcs_fit_proton.fwdMomUncertainty();
    srtrack.mcsP.bwdP_proton     = mcs_fit_proton.bwdMomentum();
    srtrack.mcsP.bwdP_err_proton = mcs_fit_proton.bwdMomUncertainty();

    srtrack.mcsP.fwdP_pion     = mcs_fit_pion.fwdMomentum();
    srtrack.mcsP.fwdP_err_pion = mcs_fit_pion.fwdMomUncertainty();
    srtrack.mcsP.bwdP_pion     = mcs_fit_pion.bwdMomentum();
    srtrack.mcsP.bwdP_err_pion = mcs_fit_pion.bwdMomUncertainty();

    srtrack.mcsP.fwdP_kaon     = mcs_fit_kaon.fwdMomentum();
    srtrack.mcsP.fwdP_err_kaon = mcs_fit_kaon.fwdMomUncertainty();
    srtrack.mcsP.bwdP_kaon     = mcs_fit_kaon.bwdMomentum();
    srtrack.mcsP.bwdP_err_kaon = mcs_fit_kaon.bwdMomUncertainty();

    // calculate range momentum
    srtrack.rangeP.p_muon = tdata.range_calculator->GetTrackMomentum(track.Length(), 13);
    srtrack.rangeP.p_proton = tdata.range_calculator->GetTrackMomentum(track.Length(), 2212);

    // get the particle ID's
    //
    // iterate over the planes -- use the conduction plane to get the particle ID
    srtrack.chi2pid.pid_ndof = 0;
    assert(tdata.particleIDs.size() == 0 || tdata.particleIDs == 3);
    for (unsigned i = 0; i < tdata.particleIDs.size(); i++) { 
      const anab::ParticleID &particle_id = *tdata.particleIDs[i];
      if (particle_id.PlaneID() && tdata.geom->SignalType(particle_id.PlaneID()) == geo::kCollection) {
        srtrack.chi2pid.chi2_muon = particle_id.Chi2Muon();
        srtrack.chi2pid.chi2_pion = particle_id.Chi2Kaon();
        srtrack.chi2pid.chi2_kaon = particle_id.Chi2Pion();
        srtrack.chi2pid.chi2_proton = particle_id.Chi2Proton();
        srtrack.chi2pid.pid_ndof = particle_id.Ndf();
      }
    }

    // bad particle ID -- set chi2 to -1
    if (srtrack.chi2pid.pid_ndof == 0) {
      srtrack.chi2pid.chi2_muon = -1;
      srtrack.chi2pid.chi2_proton = -1;
      srtrack.chi2pid.chi2_kaon = -1;
      srtrack.chi2pid.chi2_pion = -1;
    }

    // TO DO: Fill pdg value from Chi2ParticePID
    // TO DO: Make a FillTrackChi2PID

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
