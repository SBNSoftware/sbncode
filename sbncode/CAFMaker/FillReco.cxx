//////////////////////////////////////////////////////////////////////
// \file    FillReco.cxx
// \brief   Fill reco SR branches 
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#include "FillReco.h"
#include "RecoUtils/RecoUtils.h"


namespace caf
{

  //......................................................................
  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic) {
    return (slice.is_clear_cosmic || !cut_clear_cosmic) // No clear cosmics
           && slice.primary.size() > 0; // must have primary tracks/showers
  }

  //......................................................................
  void FillShowerVars(const recob::Shower& shower, 
                            caf::SRShower &srshower,
                            bool allowEmpty)
  {

    srshower.dir    = SRVector3D( shower.Direction() );
    srshower.start  = SRVector3D( shower.ShowerStart() );
    srshower.dEdx     = shower.dEdx();
    srshower.energy   = shower.Energy();

    // TO DO: work out conversion gap
    // It's sth like this but not quite. And will need to pass a simb::MCtruth object vtx position anyway.
    // srshower.conversion_gap = (shower.ShowerStart() - vertex.Position()).Mag();

    if(shower.best_plane() != -999){
      srshower.bestplane        = shower.best_plane();
      srshower.bestplane_dEdx   = shower.dEdx().at(shower.best_plane());
      srshower.bestplane_energy = shower.Energy().at(shower.best_plane());
    }

    if(shower.Length() > 0) {
      srshower.len = shower.Length();
      if(shower.best_plane() != -999){
        srshower.density = shower.Energy().at(shower.best_plane()) / shower.Length();
      }
    }

    // if(shower.has_open_angle()) {
    srshower.open_angle = shower.OpenAngle();
    // }

  }

  void FillSliceVars(const recob::Slice& slice,
                     const recob::PFParticle *primary /* can be null */,
                     caf::SRSlice &srslice,
                     bool allowEmpty)
  {

    srslice.charge       = slice.Charge();

    // get the primary tracks/showers
    if (primary != NULL) {
      for (unsigned id: primary->Daughters()) {
        srslice.primary.push_back(id);
      }
    }
  }

  void FillSliceMetadata(const larpandoraobj::PFParticleMetadata *primary_meta,
                        caf::SRSlice &srslice,
                        bool allowEmpty)
  {
    // default values   
    srslice.nu_score = -1;
    srslice.is_clear_cosmic = true;
 
    // collect the properties
    if (primary_meta != NULL) {
      auto const &properties = primary_meta->GetPropertiesMap();
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

  }

  void FillSliceFlashMatch(const anab::T0 *fmatch /* can be NULL */,
                           caf::SRSlice &srslice,
                           bool allowEmpty)
  {
    if (fmatch != NULL) {
      srslice.fmatch.present = true;
      srslice.fmatch.time = fmatch->Time();
      srslice.fmatch.score = fmatch->TriggerConfidence();
      srslice.fmatch.pe = fmatch->TriggerType();
    }
    else {
      srslice.fmatch.present = false;
    }
  }
  void FillSliceVertex(const recob::Vertex *vertex,
                       caf::SRSlice& slice,
                       bool allowEmpty) {
    if (vertex != NULL) {
      slice.vertex.x = vertex->position().X();
      slice.vertex.y = vertex->position().Y();
      slice.vertex.z = vertex->position().Z();
    }
  }


  //......................................................................

  void FillTrackCRTHit(const std::vector<art::Ptr<sbn::crt::CRTHit>> &hitmatch, 
                       const std::vector<const anab::T0*> &t0match, 
                       caf::SRTrack &srtrack,
                       bool allowEmpty)
  {
    if (hitmatch.size()) {
      assert(hitmatch.size() == 1);
      assert(t0match.size() == 1);
      srtrack.crthit.distance = t0match[0]->fTriggerConfidence; 
      srtrack.crthit.hit.time = t0match[0]->fTime;
      srtrack.crthit.hit.position.x = hitmatch[0]->x_pos;
      srtrack.crthit.hit.position.y = hitmatch[0]->y_pos;
      srtrack.crthit.hit.position.z = hitmatch[0]->z_pos;
      srtrack.crthit.hit.position_err.x = hitmatch[0]->x_err;
      srtrack.crthit.hit.position_err.y = hitmatch[0]->y_err;
      srtrack.crthit.hit.position_err.z = hitmatch[0]->z_err;

    }
  }

  void FillTrackMCS(const recob::Track& track,
                    const std::array<std::vector<art::Ptr<recob::MCSFitResult>>, 4> &mcs_results,
                    caf::SRTrack& srtrack,
                    bool allowEmpty)
  {
    // gather MCS fits
    if (mcs_results[0].size()) {
      recob::MCSFitResult mcs_fit_muon = *mcs_results[0][0];

      srtrack.mcsP.fwdP_muon     = mcs_fit_muon.fwdMomentum();
      srtrack.mcsP.fwdP_err_muon = mcs_fit_muon.fwdMomUncertainty();
      srtrack.mcsP.bwdP_muon     = mcs_fit_muon.bwdMomentum();
      srtrack.mcsP.bwdP_err_muon = mcs_fit_muon.bwdMomUncertainty();
    }

    if (mcs_results[1].size()) {
      recob::MCSFitResult mcs_fit_proton = *mcs_results[1][0];

      srtrack.mcsP.fwdP_proton     = mcs_fit_proton.fwdMomentum();
      srtrack.mcsP.fwdP_err_proton = mcs_fit_proton.fwdMomUncertainty();
      srtrack.mcsP.bwdP_proton     = mcs_fit_proton.bwdMomentum();
      srtrack.mcsP.bwdP_err_proton = mcs_fit_proton.bwdMomUncertainty();
    }

    if (mcs_results[2].size()) {
      recob::MCSFitResult mcs_fit_pion = *mcs_results[2][0];

      srtrack.mcsP.fwdP_pion     = mcs_fit_pion.fwdMomentum();
      srtrack.mcsP.fwdP_err_pion = mcs_fit_pion.fwdMomUncertainty();
      srtrack.mcsP.bwdP_pion     = mcs_fit_pion.bwdMomentum();
      srtrack.mcsP.bwdP_err_pion = mcs_fit_pion.bwdMomUncertainty();
    }

    if (mcs_results[3].size()) {
      recob::MCSFitResult mcs_fit_kaon = *mcs_results[3][0];

      srtrack.mcsP.fwdP_kaon     = mcs_fit_kaon.fwdMomentum();
      srtrack.mcsP.fwdP_err_kaon = mcs_fit_kaon.fwdMomUncertainty();
      srtrack.mcsP.bwdP_kaon     = mcs_fit_kaon.bwdMomentum();
      srtrack.mcsP.bwdP_err_kaon = mcs_fit_kaon.bwdMomUncertainty();
    }
  }

  void FillTrackRangeP(const recob::Track& track,
                       const std::array<std::vector<art::Ptr<sbn::RangeP>>, 2> &range_results,
                       caf::SRTrack& srtrack,
                       bool allowEmpty)
  {
    // calculate range momentum
    if (range_results[0].size()) {
      srtrack.rangeP.p_muon = range_results[0][0]->range_p; 
      assert(track.ID() == range_results[0][0]->trackID);
    }

    if (range_results[1].size()) {
      srtrack.rangeP.p_proton = range_results[1][0]->range_p; 
      assert(track.ID() == range_results[1][0]->trackID);
    }
  }

  void FillTrackChi2PID(const std::vector<art::Ptr<anab::ParticleID>> particleIDs,
                        const geo::GeometryCore *geom,
                        caf::SRTrack& srtrack,
                        bool allowEmpty)
  {
    // get the particle ID's
    //
    // iterate over the planes -- use the conduction plane to get the particle ID
    srtrack.chi2pid.pid_ndof = 0;
    assert(particleIDs.size() == 0 || particleIDs == 3);
    for (unsigned i = 0; i < particleIDs.size(); i++) { 
      const anab::ParticleID &particle_id = *particleIDs[i];
      if (particle_id.PlaneID() && geom->SignalType(particle_id.PlaneID()) == geo::kCollection) {
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

  }

  void FillTrackCalo(const std::vector<art::Ptr<anab::Calorimetry>> &calos,
                     const geo::GeometryCore *geom,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {
    // TODO: what to do with calorimetry
  }

  // TODO: crt matching

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle &particle,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {

    srtrack.npts = track.CountValidPoints();
    srtrack.len  = track.Length();
    srtrack.costh = track.StartDirection().Z() / sqrt(track.StartDirection().Mag2());

    srtrack.ID = particle.Self();

    // set the daughters in the particle flow
    for (unsigned id: particle.Daughters()) {
      srtrack.daughters.push_back(id);
    }

  }
  //......................................................................
  
  // TODO: implement
  void SetNuMuCCPrimary(std::vector<caf::StandardRecord> &recs,
                        std::vector<caf::SRTrueInteraction> &srneutrinos) {}

} // end namespace 

