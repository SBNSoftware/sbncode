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
    return (!slice.is_clear_cosmic || !cut_clear_cosmic) // No clear cosmics
           && slice.primary.size() > 0; // must have primary tracks/showers
  }

  void FillCRTHit(const sbn::crt::CRTHit &hit,
                  bool use_ts0,
                  caf::SRCRTHit &srhit,
                  bool allowEmpty) {
    srhit.time = (use_ts0 ? (float)hit.ts0_ns : hit.ts1_ns) / 1000.;

    srhit.position.x = hit.x_pos;
    srhit.position.y = hit.y_pos;
    srhit.position.z = hit.z_pos;

    srhit.position_err.x = hit.x_err;
    srhit.position_err.y = hit.y_err;
    srhit.position_err.z = hit.z_err;

    srhit.pe = hit.peshit;
    srhit.plane = hit.plane;

  }

  void FillCRTTrack(const sbn::crt::CRTTrack &track,
                  bool use_ts0,
                  caf::SRCRTTrack &srtrack,
                  bool allowEmpty) {

    srtrack.time = (use_ts0 ? (float)track.ts0_ns : track.ts1_ns) / 1000.;

    srtrack.hita.position.x = track.x1_pos;
    srtrack.hita.position.y = track.y1_pos;
    srtrack.hita.position.z = track.z1_pos;

    srtrack.hita.position_err.x = track.x1_err;
    srtrack.hita.position_err.y = track.y1_err;
    srtrack.hita.position_err.z = track.z1_err;

    srtrack.hita.plane = track.plane1;

    srtrack.hitb.position.x = track.x2_pos;
    srtrack.hitb.position.y = track.y2_pos;
    srtrack.hitb.position.z = track.z2_pos;

    srtrack.hitb.position_err.x = track.x2_err;
    srtrack.hitb.position_err.y = track.y2_err;
    srtrack.hitb.position_err.z = track.z2_err;

    srtrack.hitb.plane = track.plane2;
  }


  std::vector<float> double_to_float_vector(const std::vector<double>& v)
  {
    std::vector<float> ret;
    ret.reserve(v.size());
    for(double x: v) ret.push_back(x);
    return ret;
  }

  //......................................................................
  void FillShowerVars(const recob::Shower& shower,
                      const recob::Vertex* vertex,
                      const std::vector<art::Ptr<recob::Hit>> &hits,
                      const geo::GeometryCore *geom,
                      unsigned producer,
                      caf::SRShower &srshower,
                      bool allowEmpty)
  {

    srshower.producer = producer;
    srshower.dEdx_plane0 = (shower.dEdx())[0];
    srshower.dEdx_plane1 = (shower.dEdx())[1];
    srshower.dEdx_plane2 = (shower.dEdx())[2];
    srshower.energy_plane0 = (shower.Energy())[0];
    srshower.energy_plane1 = (shower.Energy())[1];
    srshower.energy_plane2 = (shower.Energy())[2];
    srshower.dEdx   = double_to_float_vector( shower.dEdx() );
    srshower.energy = double_to_float_vector( shower.Energy() );
    srshower.dir    = SRVector3D( shower.Direction() );
    srshower.start  = SRVector3D( shower.ShowerStart() );

    // TO DO: work out conversion gap
    // It's sth like this but not quite. And will need to pass a simb::MCtruth object vtx position anyway.
    // srshower.conversion_gap = (shower.ShowerStart() - vertex.Position()).Mag();

    if(shower.best_plane() != -999){
      srshower.bestplane        = shower.best_plane();
      srshower.bestplane_dEdx   = shower.dEdx().at(shower.best_plane());
      srshower.bestplane_energy = shower.Energy().at(shower.best_plane());
    }

    if(shower.has_open_angle())
      srshower.open_angle = shower.OpenAngle();
    if(shower.has_length())
      srshower.len = shower.Length();

    if(srshower.len > std::numeric_limits<float>::epsilon() && srshower.bestplane_energy > 0)
        srshower.density = srshower.bestplane_energy / srshower.len;

    if (vertex && shower.ShowerStart().Z()>-990) {
      // Need to do some rearranging to make consistent types
      const geo::Point_t vertexPos(vertex->position());
      const TVector3 vertexTVec3{vertexPos.X(), vertexPos.Y(), vertexPos.Z()};

      srshower.conversion_gap = (shower.ShowerStart() - vertexTVec3).Mag();
    }

    if (shower.Direction().Z()>-990 && shower.ShowerStart().Z()>-990 && shower.Length()>0) {
      srshower.end = shower.ShowerStart()+ (shower.Length() * shower.Direction());
    }

    for (auto const& hit:hits) {
      switch (hit->WireID().Plane) {
        case(0):
          srshower.nHits_plane0++;
        case(1):
          srshower.nHits_plane1++;
        case(2):
          srshower.nHits_plane2++;
      }
    }

    for (geo::PlaneGeo const& plane: geom->IteratePlanes()) {

      const double angleToVert(geom->WireAngleToVertical(plane.View(), plane.ID()) - 0.5*M_PI);
      const double cosgamma(std::abs(std::sin(angleToVert)*shower.Direction().Y()+std::cos(angleToVert)*shower.Direction().Z()));

      switch (plane.ID().Plane) {
        case(0):
          srshower.wirePitch_plane0 = plane.WirePitch()/cosgamma;
        case(1):
          srshower.wirePitch_plane1 = plane.WirePitch()/cosgamma;
        case(2):
          srshower.wirePitch_plane2 = plane.WirePitch()/cosgamma;
      }
    }
  }

  void FillShowerMVAPID(const art::Ptr<sbn::MVAPID> mvaPID,
      caf::SRShower& srshower,
      bool allowEmpty)
  {
    srshower.mvaPID.electronScore = mvaPID->mMVAScoreMap.at(11);
    srshower.mvaPID.photonScore = mvaPID->mMVAScoreMap.at(22);
    srshower.mvaPID.otherScore = mvaPID->mMVAScoreMap.at(0);

    srshower.mvaPID.pdg = mvaPID->BestPDG();
    srshower.mvaPID.bestScore = mvaPID->BestScore();
  }


  void FillShowerCosmicCylinder(const std::vector<art::Ptr<float> >& cosmicCylinderVec,
                      caf::SRShower& srshower)
  {
      if (cosmicCylinderVec.size() != 1)
        return;
      srshower.cosmicCylinder = *cosmicCylinderVec.front();
  }

  void FillShowerResiduals(const std::vector<art::Ptr<float> >& residuals,
                      caf::SRShower& srshower)
  {
    for (auto const& res: residuals) {
      srshower.selVars.showerResiduals.push_back(*res);
    }
  }

  void FillShowerTrackFit(const sbn::ShowerTrackFit& trackFit,
                      caf::SRShower& srshower)
  {
    srshower.selVars.trackLength = trackFit.mTrackLength;
    srshower.selVars.trackWidth  = trackFit.mTrackWidth;
  }

  void FillShowerDensityFit(const sbn::ShowerDensityFit& densityFit,
                      caf::SRShower& srshower)
  {
    srshower.selVars.densityGradient      = densityFit.mDensityGrad;
    srshower.selVars.densityGradientPower = densityFit.mDensityPow;
  }

  void FillSliceVars(const recob::Slice& slice,
                     const recob::PFParticle *primary /* can be null */,
                     unsigned producer,
                     caf::SRSlice &srslice,
                     bool allowEmpty)
  {

    srslice.producer = producer;
    srslice.charge       = slice.Charge();

    // get the primary tracks/showers
    if (primary != NULL) {
      for (unsigned id: primary->Daughters()) {
        srslice.primary.push_back(id);
      }
      srslice.self = primary->Self();
      srslice.nu_pdg = primary->PdgCode();
    }
    else {
      srslice.self = -1;
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
        srslice.is_clear_cosmic = true;
      }
      else {
        assert(properties.count("IsNeutrino"));
        srslice.is_clear_cosmic = false;
      }
      if (properties.count("NuScore")) {
        srslice.nu_score = properties.at("NuScore");
      }
      else {
        srslice.nu_score = -1;
      }
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

  void FillTrackCRTHit(const std::vector<art::Ptr<anab::T0>> &t0match,
                       caf::SRTrack &srtrack,
                       bool allowEmpty)
  {
    if (t0match.size()) {
      assert(t0match.size() == 1);
      srtrack.crthit.distance = t0match[0]->fTriggerConfidence;
      srtrack.crthit.hit.time = t0match[0]->fTime / 1e3; /* ns -> us */

      // TODO/FIXME: FILL THESE ONCE WE HAVE THE CRT HIT!!!
      // srtrack.crthit.hit.position.x = hitmatch[0]->x_pos;
      // srtrack.crthit.hit.position.y = hitmatch[0]->y_pos;
      // srtrack.crthit.hit.position.z = hitmatch[0]->z_pos;
      // srtrack.crthit.hit.position_err.x = hitmatch[0]->x_err;
      // srtrack.crthit.hit.position_err.y = hitmatch[0]->y_err;
      // srtrack.crthit.hit.position_err.z = hitmatch[0]->z_err;

    }
  }

  void FillTrackCRTTrack(const std::vector<art::Ptr<anab::T0>> &t0match,
                       caf::SRTrack &srtrack,
                       bool allowEmpty)
  {
    if (t0match.size()) {
      assert(t0match.size() == 1);
      srtrack.crttrack.angle = t0match[0]->fTriggerConfidence;
      srtrack.crttrack.time = t0match[0]->fTime / 1e3; /* ns -> us */

      // TODO/FIXME: FILL MORE ONCE WE HAVE THE CRT HIT!!!

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
                       const std::array<std::vector<art::Ptr<sbn::RangeP>>, 3> &range_results,
                       caf::SRTrack& srtrack,
                       bool allowEmpty)
  {
    if (range_results[0].size()) {
      srtrack.rangeP.p_muon = range_results[0][0]->range_p;
      assert(track.ID() == range_results[0][0]->trackID);
    }

    if (range_results[1].size()) {
      srtrack.rangeP.p_pion = range_results[1][0]->range_p;
      assert(track.ID() == range_results[1][0]->trackID);
    }

    if (range_results[2].size()) {
      srtrack.rangeP.p_proton = range_results[1][0]->range_p;
      assert(track.ID() == range_results[1][0]->trackID);
    }
  }

  void FillPlaneChi2PID(const anab::ParticleID &particle_id, caf::SRTrkChi2PID &srpid) {
    srpid.chi2_muon = particle_id.Chi2Muon();
    srpid.chi2_pion = particle_id.Chi2Pion();
    srpid.chi2_kaon = particle_id.Chi2Kaon();
    srpid.chi2_proton = particle_id.Chi2Proton();
    srpid.pid_ndof = particle_id.Ndf();
    srpid.pida = particle_id.PIDA();
  }

  void FillTrackChi2PID(const std::vector<art::Ptr<anab::ParticleID>> particleIDs,
                        const geo::GeometryCore *geom,
                        caf::SRTrack& srtrack,
                        bool allowEmpty)
  {
    // get the particle ID's
    for (unsigned i = 0; i < particleIDs.size(); i++) {
      const anab::ParticleID &particle_id = *particleIDs[i];
      if (particle_id.PlaneID()) {
        unsigned plane_id  = particle_id.PlaneID().Plane;
        assert(plane_id < 3);
        caf::SRTrkChi2PID &this_pid = (plane_id == 0) ? srtrack.chi2pid0 : ((plane_id == 1) ? srtrack.chi2pid1 : srtrack.chi2pid2);
        FillPlaneChi2PID(particle_id, this_pid);
      }
    }
  }

  void FillTrackPlaneCalo(const anab::Calorimetry &calo, float constant, caf::SRTrackCalo &srcalo) {
    const std::vector<float> &dqdx = calo.dQdx();
    const std::vector<float> &dedx = calo.dEdx();
    const std::vector<float> &pitch = calo.TrkPitchVec();
    srcalo.charge = 0.;
    srcalo.ke = 0.;
    srcalo.nhit = 0;
    for (unsigned i = 0; i < dedx.size(); i++) {
      if (dedx[i] > 1000.) continue;
      srcalo.nhit ++;
      srcalo.charge += dqdx[i] * pitch[i] / constant; /* convert ADC*tick to electrons */
      srcalo.ke += dedx[i] * pitch[i];
    }
  }

  void FillTrackScatterDCA(const art::Ptr<sbn::ScatterDCA> dca,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.scatterDCA.mean = dca->mMean;
    srtrack.scatterDCA.stdDev = dca->mStdDev;
    srtrack.scatterDCA.max = dca->mMax;
  }

  void FillTrackStoppingChi2Fit(const art::Ptr<sbn::StoppingChi2Fit> stoppingChi2,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.stoppingChi2Fit.pol0Chi2 = stoppingChi2->mPol0Chi2;
    srtrack.stoppingChi2Fit.expChi2 = stoppingChi2->mExpChi2;
    srtrack.stoppingChi2Fit.pol0Fit = stoppingChi2->mPol0Fit;
  }

  void FillTrackMVAPID(const art::Ptr<sbn::MVAPID> mvaPID,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.mvaPID.muonScore = mvaPID->mMVAScoreMap.at(13);
    srtrack.mvaPID.pionScore = mvaPID->mMVAScoreMap.at(211);
    srtrack.mvaPID.protonScore = mvaPID->mMVAScoreMap.at(2212);
    srtrack.mvaPID.otherScore = mvaPID->mMVAScoreMap.at(0);

    srtrack.mvaPID.pdg = mvaPID->BestPDG();
    srtrack.mvaPID.bestScore = mvaPID->BestScore();
  }

  void FillTrackCalo(const std::vector<art::Ptr<anab::Calorimetry>> &calos,
                     const geo::GeometryCore *geom,
                     const std::array<float, 3> &calo_constants,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {
    // count up the kinetic energy on each plane --
    // ignore any charge with a deposition > 1000 MeV/cm
    // TODO: ignore first and last hit???
    //    assert(calos.size() == 0 || calos == 3);
    for (unsigned i = 0; i < calos.size(); i++) {
      const anab::Calorimetry &calo = *calos[i];
      if (calo.PlaneID()) {
        unsigned plane_id = calo.PlaneID().Plane;
        assert(plane_id < 3);
        caf::SRTrackCalo &this_calo = (plane_id == 0) ? srtrack.calo0 : ((plane_id == 1) ? srtrack.calo1 : srtrack.calo2);
        FillTrackPlaneCalo(calo, calo_constants[plane_id], this_calo);
      }
    }

    // Set the plane with the most hits
    //
    // We expect the noise to be lowest at planes 2 -> 0 -> 1, so use this to break ties
    caf::Plane_t bestplane = caf::kUnknown;
    if (srtrack.calo2.nhit > 0 && srtrack.calo2.nhit >= srtrack.calo0.nhit && srtrack.calo2.nhit >=  srtrack.calo1.nhit) {
      bestplane = (caf::Plane_t)2;
    }
    else if (srtrack.calo0.nhit > 0 && srtrack.calo0.nhit >= srtrack.calo2.nhit && srtrack.calo0.nhit >=  srtrack.calo1.nhit) {
      bestplane = (caf::Plane_t)0;
    }
    else if (srtrack.calo1.nhit > 1) {
      bestplane = (caf::Plane_t)1;
    }
    srtrack.bestplane = bestplane;

  }

  // TODO: crt matching

  void FillTrackVars(const recob::Track& track,
                     unsigned producer,
                     caf::SRTrack& srtrack,
                     bool allowEmpty)
  {

    srtrack.producer = producer;
    srtrack.npts = track.CountValidPoints();
    srtrack.len  = track.Length();
    srtrack.costh = track.StartDirection().Z() / sqrt(track.StartDirection().Mag2());
    srtrack.phi = track.StartDirection().Phi();

    srtrack.dir_end.x = track.EndDirection().X();
    srtrack.dir_end.y = track.EndDirection().Y();
    srtrack.dir_end.z = track.EndDirection().Z();

    srtrack.dir.x = track.StartDirection().X();
    srtrack.dir.y = track.StartDirection().Y();
    srtrack.dir.z = track.StartDirection().Z();

    srtrack.start.x = track.Start().X();
    srtrack.start.y = track.Start().Y();
    srtrack.start.z = track.Start().Z();

    srtrack.end.x = track.End().X();
    srtrack.end.y = track.End().Y();
    srtrack.end.z = track.End().Z();

  }

  void FillPFPVars(const recob::PFParticle &particle,
                   const recob::PFParticle *primary,
                   const larpandoraobj::PFParticleMetadata *pfpMeta,
                   caf::SRPFP& srpfp,
                   bool allowEmpty)
  {
    srpfp.ID = particle.Self();
    srpfp.slcID = (primary) ? primary->Self() : -1;

    // set the daughters in the particle flow
    for (unsigned id: particle.Daughters()) {
      srpfp.daughters.push_back(id);
    }
    srpfp.ndaughters = srpfp.daughters.size();

    srpfp.parent = particle.Parent();
    srpfp.parent_is_primary = (particle.Parent() == recob::PFParticle::kPFParticlePrimary) \
      || (primary && particle.Parent() == primary->Self());

    if (pfpMeta) {
      auto const &propertiesMap (pfpMeta->GetPropertiesMap());
      auto const &pfpTrackScoreIter(propertiesMap.find("TrackScore"));
      srpfp.trackScore = (pfpTrackScoreIter == propertiesMap.end()) ? -5.f : pfpTrackScoreIter->second;
    }
  }

  //......................................................................

  void SetNuMuCCPrimary(std::vector<caf::StandardRecord> &recs,
                        std::vector<caf::SRTrueInteraction> &srneutrinos) {
  //   // set is_primary to true by default
  //   for (caf::StandardRecord &rec: recs) {
  //     rec.slc.tmatch.is_numucc_primary = true;
  //   }

  //   for (unsigned i = 0; i < srneutrinos.size(); i++) {
  //     ApplyNumuCCMatching(recs, srneutrinos, i);
  //   }
  }

  void ApplyNumuCCMatching(std::vector<caf::StandardRecord> &recs,
                           const std::vector<caf::SRTrueInteraction> &srneutrinos,
                           unsigned truth_ind) {

  //   std::vector<unsigned> matches_truth;
  //   for (unsigned i = 0; i < recs.size(); i++) {
  //     if (recs[i].slc.tmatch.index == (int)truth_ind) {
  //       matches_truth.push_back(i);
  //     }
  //   }

  //   // first -- remove any cases where most of the slice
  //   // matches to non-primary particles of the neutrino
  //   unsigned ind = 0;
  //   std::vector<float> matching_primary_energy;
  //   while (ind < matches_truth.size()) {
  //     const caf::SRSliceRecoBranch &reco = recs[matches_truth[ind]].reco;
  //     const caf::SRSlice &slice = recs[matches_truth[ind]].slc;

  //     caf::SRVector3D vertex = slice.vertex;

  //     float primary_energy = 0.;
  //     float total_energy = 0.;

  //     // check the primary tracks of the slice
  //     for (const caf::SRTrack &track: reco.trk) {
  //       caf::SRVector3D start = track.start;
  //       float dist = sqrt((start.x - vertex.x) * (start.x - vertex.x) +
  //                         (start.y - vertex.y) * (start.y - vertex.y) +
  //                         (start.z - vertex.z) * (start.z - vertex.z));

  //       if (track.parent == slice.self && dist < 10.) {
  //         for (const caf::SRTrackTruth::ParticleMatch &pmatch: track.truth.matches) {
  //           total_energy += pmatch.energy;
  //           for (unsigned i_part = 0; i_part < recs[0].true_particles.size(); i_part++) {
  //             const caf::SRTrueParticle &particle = recs[0].true_particles[i_part];
  //             if (particle.G4ID == pmatch.G4ID) {
  //               if (particle.start_process == caf::kG4primary) {
  //                 primary_energy += pmatch.energy;
  //               }
  //               break;
  //             }
  //           }
  //         }
  //       }
  //     }
  //     if (primary_energy / total_energy < 0.5) {
  //       recs[matches_truth[ind]].slc.tmatch.is_numucc_primary = false;
  //       matches_truth.erase(matches_truth.begin()+ind);
  //     }
  //     else {
  //       matching_primary_energy.push_back(primary_energy);
  //       ind ++;
  //     }
  //   }

  //   // less than two matches! All good
  //   if (matches_truth.size() < 2) return;

  //   // If this is a numu CC interaction, break
  //   // tie by matching the muon
  //   // Whoever has a track matching closer to the
  //   // start of the muon wins
  //   if (abs(srneutrinos[truth_ind].pdg == 14) && srneutrinos[truth_ind].iscc) {
  //     const caf::SRTrueParticle &muon = srneutrinos[truth_ind].prim[0];
  //     float closest_dist = -1;
  //     int best_index = -1;
  //     for (unsigned ind = 0; ind < matches_truth.size(); ind++) {
  //       const caf::SRSliceRecoBranch &reco = recs[matches_truth[ind]].reco;
  //       const caf::SRSlice &slice = recs[matches_truth[ind]].slc;

  //       caf::SRVector3D vertex = slice.vertex;

  //       for (const caf::SRTrack &track: reco.trk) {
  //         caf::SRVector3D start = track.start;
  //         float dist = sqrt((start.x - vertex.x) * (start.x - vertex.x) +
  //                           (start.y - vertex.y) * (start.y - vertex.y) +
  //                           (start.z - vertex.z) * (start.z - vertex.z));

  //         if (track.parent == slice.self && dist < 10. && track.truth.matches.size()) {
  //           const caf::SRTrackTruth::ParticleMatch &pmatch = track.truth.matches[0];
  //           if (pmatch.energy / muon.planeVisE > 0.05 && pmatch.G4ID == muon.G4ID) {
  //              caf::SRVector3D start = track.start;
  //              caf::SRVector3D end = track.end;
  //              float start_dist = sqrt((start.x - muon.start.x) * (start.x - muon.start.x) +
  //                                      (start.y - muon.start.y) * (start.y - muon.start.y) +
  //                                      (start.z - muon.start.z) * (start.z - muon.start.z));
  //              float end_dist = sqrt((end.x - muon.start.x) * (end.x - muon.start.x) +
  //                                    (end.y - muon.start.y) * (end.y - muon.start.y) +
  //                                    (end.z - muon.start.z) * (end.z - muon.start.z));
  //              float this_dist = std::min(start_dist, end_dist);
  //              if (closest_dist < 0. || this_dist < closest_dist) {
  //                closest_dist = this_dist;
  //                best_index = ind;
  //              }
  //           }
  //         }
  //       }
  //     }

  //     // found a match!
  //     if (best_index >= 0) {
  //       for (unsigned i = 0; i < matches_truth.size(); i++) {
  //         if ((int)i == best_index) recs[matches_truth[i]].slc.tmatch.is_numucc_primary = true;
  //         else                 recs[matches_truth[i]].slc.tmatch.is_numucc_primary = false;
  //       }
  //       return;
  //     }
  //     // no match :( fallback on non numu-CC matching
  //     else {}
  //   }

  //   // Otherwise, take the most energetic one
  //   unsigned best_index = std::distance(matching_primary_energy.begin(),
  //                                       std::max_element(matching_primary_energy.begin(), matching_primary_energy.end()));

  //   for (unsigned i = 0; i < matches_truth.size(); i++) {
  //     if (i == best_index) recs[matches_truth[i]].slc.tmatch.is_numucc_primary = true;
  //     else                 recs[matches_truth[i]].slc.tmatch.is_numucc_primary = false;
  //   }
    return;
  }


} // end namespace
