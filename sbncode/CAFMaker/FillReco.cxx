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

  void FillStubVars(const sbn::Stub &stub,
                    const art::Ptr<recob::PFParticle> stubpfp,
                    caf::SRStub &srstub,
                    bool allowEmpty) {

    // allowEmpty does not (yet) matter here
    (void) allowEmpty;

    // Copy the stub object over
    srstub.vtx.x = stub.vtx.x();
    srstub.vtx.y = stub.vtx.y();
    srstub.vtx.z = stub.vtx.z();

    srstub.end.x = stub.end.x();
    srstub.end.y = stub.end.y();
    srstub.end.z = stub.end.z();

    srstub.efield_vtx = stub.efield_vtx;
    srstub.efield_end = stub.efield_end;

    for (unsigned ip = 0; ip < stub.plane.size(); ip++) {
      caf::SRStubPlane plane;
      plane.p = (caf::Plane_t)stub.plane[ip].Plane;
      plane.pitch = stub.pitch[ip];
      plane.trkpitch = stub.trkpitch[ip];
      plane.vtx_w = stub.vtx_w[ip];
      plane.hit_w = stub.hit_w[ip];
      for (unsigned ih = 0; ih < stub.hits[ip].size(); ih++) {
        caf::SRStubHit hit;
        hit.charge = stub.hits[ip][ih].charge;
        hit.wire = stub.hits[ip][ih].wire;
        hit.ontrack = stub.hits[ip][ih].ontrack;

        plane.hits.push_back(hit);
      }

      srstub.planes.push_back(plane);     
    }

    // If there is an overlaid PFParticle, save its ID
    if (stubpfp) srstub.pfpid = stubpfp->Self();

  }

  void FillCRTHit(const sbn::crt::CRTHit &hit,
                  bool use_ts0,
                  int64_t CRT_T0_reference_time, // ns, signed
                  double CRT_T1_reference_time, // us
                  caf::SRCRTHit &srhit,
                  bool allowEmpty) {

    srhit.t0 = ( (long long)(hit.ts0()) /*u_int64_t to int64_t*/ + CRT_T0_reference_time )/1000.;
    srhit.t1 = hit.ts1()/1000.+CRT_T1_reference_time; // ns -> us
    srhit.time = use_ts0 ? srhit.t0 : srhit.t1;

    srhit.position.x = hit.x_pos;
    srhit.position.y = hit.y_pos;
    srhit.position.z = hit.z_pos;

    srhit.position_err.x = hit.x_err;
    srhit.position_err.y = hit.y_err;
    srhit.position_err.z = hit.z_err;

    srhit.pe = hit.peshit;
    srhit.plane = hit.plane;

    srhit.flag = (hit.plane>34 && hit.plane<50 && hit.ts0_s_corr<7) ? hit.ts0_s_corr : 0 ;
    
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

  void FillCRTPMTMatch(const sbn::crt::CRTPMTMatching &match,
		       caf::SRCRTPMTMatch &srmatch,
		       bool allowEmpty){
    // allowEmpty does not (yet) matter here                                                           
    (void) allowEmpty;
    //srmatch.setDefault();
    srmatch.flashID = match.flashID;
    srmatch.flashTime_us = match.flashTime;
    srmatch.flashGateTime = match.flashGateTime;
    srmatch.firstOpHitPeakTime = match.firstOpHitPeakTime;
    srmatch.firstOpHitStartTime = match.firstOpHitStartTime;
    srmatch.flashInGate = match.flashInGate;
    srmatch.flashInBeam = match.flashInBeam;
    srmatch.flashPE = match.flashPE;
    srmatch.flashPosition = SRVector3D (match.flashPosition.X(), match.flashPosition.Y(), match.flashPosition.Z());
    srmatch.flashYWidth = match.flashYWidth;
    srmatch.flashZWidth = match.flashZWidth;
    srmatch.flashClassification = static_cast<int>(match.flashClassification);
    for(const auto& matchedCRTHit : match.matchedCRTHits){
      std::cout << "CRTPMTTimeDiff = "<< matchedCRTHit.PMTTimeDiff << "\n";
      caf::SRMatchedCRT matchedCRT;
      matchedCRT.PMTTimeDiff = matchedCRTHit.PMTTimeDiff; 
      matchedCRT.time = matchedCRTHit.time;
      matchedCRT.sys = matchedCRTHit.sys;
      matchedCRT.region = matchedCRTHit.region;
      matchedCRT.position = SRVector3D(matchedCRTHit.position.X(), matchedCRTHit.position.Y(), matchedCRTHit.position.Z());
      srmatch.matchedCRTHits.push_back(matchedCRT);
    }
  }


  void FillOpFlash(const recob::OpFlash &flash,
                  std::vector<recob::OpHit const*> const& hits,
                  int cryo, 
                  caf::SROpFlash &srflash,
                  bool allowEmpty) {

    srflash.setDefault();

    srflash.time = flash.Time();
    srflash.timewidth = flash.TimeWidth();

    double firstTime = std::numeric_limits<double>::max();
    for(const auto& hit: hits){
      double const hitTime = hit->HasStartTime()? hit->StartTime(): hit->PeakTime();
      if (firstTime > hitTime)
        firstTime = hitTime;
    }
    srflash.firsttime = firstTime;

    srflash.cryo = cryo; // 0 in SBND, 0/1 for E/W in ICARUS

    // Sum over each wall, not very SBND-compliant
    float sumEast = 0.;
    float sumWest = 0.;
    int countingOffset = 0;
    if ( cryo == 1 ) countingOffset += 180;
    for ( int PMT = 0 ; PMT < 180 ; PMT++ ) {
      if ( PMT <= 89 ) sumEast += flash.PEs().at(PMT + countingOffset);
      else sumWest += flash.PEs().at(PMT + countingOffset);
    }
    srflash.peperwall[0] = sumEast;
    srflash.peperwall[1] = sumWest;

    srflash.totalpe = flash.TotalPE();
    srflash.fasttototal = flash.FastToTotal();
    srflash.onbeamtime = flash.OnBeamTime();

    srflash.center.SetXYZ( -9999.f, flash.YCenter(), flash.ZCenter() );
    srflash.width.SetXYZ( -9999.f, flash.YWidth(), flash.ZWidth() );

    // Checks if ( recob::OpFlash.XCenter() != std::numeric_limits<double>::max() )
    // See LArSoft OpFlash.h at https://nusoft.fnal.gov/larsoft/doxsvn/html/OpFlash_8h_source.html
    if ( flash.hasXCenter() ) {
      srflash.center.SetX( flash.XCenter() );
      srflash.width.SetX( flash.XWidth() );
    }
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

    // We need to convert the energy from MeV to GeV
    // Also convert -999 -> -5 for consistency with other defaults in the CAFs
    for(int i = 0; i < 3; ++i){
      const float e = shower.Energy()[i];
      srshower.plane[i].energy = e > 0 ? e / 1000.f : -5.f;
      srshower.plane[i].dEdx = shower.dEdx()[i];
    }

    srshower.dir    = SRVector3D( shower.Direction() );
    srshower.start  = SRVector3D( shower.ShowerStart() );

    // TO DO: work out conversion gap
    // It's sth like this but not quite. And will need to pass a simb::MCtruth object vtx position anyway.
    // srshower.conversion_gap = (shower.ShowerStart() - vertex.Position()).Mag();

    if(shower.best_plane() != -999){
      srshower.bestplane        = shower.best_plane();
      srshower.bestplane_dEdx   = srshower.plane[shower.best_plane()].dEdx;
      srshower.bestplane_energy = srshower.plane[shower.best_plane()].energy;
    }

    if(shower.has_open_angle())
      srshower.open_angle = shower.OpenAngle();
    if(shower.has_length())
      srshower.len = shower.Length();

    // We want density to be in MeV/cm so need to convert the energy back to MeV from GeV
    if(srshower.len > std::numeric_limits<float>::epsilon() && srshower.bestplane_energy > 0)
        srshower.density = 1000.f * srshower.bestplane_energy / srshower.len;

    if (vertex && shower.ShowerStart().Z()>-990) {
      // Need to do some rearranging to make consistent types
      const geo::Point_t vertexPos(vertex->position());
      const TVector3 vertexTVec3{vertexPos.X(), vertexPos.Y(), vertexPos.Z()};

      srshower.conversion_gap = (shower.ShowerStart() - vertexTVec3).Mag();
    }

    if (shower.Direction().Z()>-990 && shower.ShowerStart().Z()>-990 && shower.Length()>0) {
      srshower.end = shower.ShowerStart()+ (shower.Length() * shower.Direction());
    }

    for(int p = 0; p < 3; ++p) srshower.plane[p].nHits = 0;
    for (auto const& hit:hits) ++srshower.plane[hit->WireID().Plane].nHits;

    for (geo::PlaneGeo const& plane: geom->Iterate<geo::PlaneGeo>()) {

      const double angleToVert(geom->WireAngleToVertical(plane.View(), plane.ID()) - 0.5*M_PI);
      const double cosgamma(std::abs(std::sin(angleToVert)*shower.Direction().Y()+std::cos(angleToVert)*shower.Direction().Z()));

      srshower.plane[plane.ID().Plane].wirePitch = plane.WirePitch()/cosgamma;
    }
  }

  void FillShowerRazzle(const art::Ptr<sbn::MVAPID> razzle,
      caf::SRShower& srshower,
      bool allowEmpty)
  {
    srshower.razzle.electronScore = razzle->mvaScoreMap.at(11);
    srshower.razzle.photonScore = razzle->mvaScoreMap.at(22);
    srshower.razzle.otherScore = razzle->mvaScoreMap.at(0);

    srshower.razzle.pdg = razzle->BestPDG();
    srshower.razzle.bestScore = razzle->BestScore();
  }


  void FillShowerCosmicDist(const std::vector<art::Ptr<float> >& cosmicDistVec,
                      caf::SRShower& srshower)
  {
      if (cosmicDistVec.size() != 1)
        return;
      srshower.cosmicDist = *cosmicDistVec.front();
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
    srslice.nuid.setDefault();

    // collect the properties
    if (primary_meta != NULL) {
      auto const &properties = primary_meta->GetPropertiesMap();
      if (properties.count("IsClearCosmic") || (properties.count("NuScore") && properties.at("NuScore") < 0)) {
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
      // NeutrinoID (SliceID) features
      CopyPropertyIfSet(properties, "NuNFinalStatePfos",        srslice.nuid.nufspfos);
      CopyPropertyIfSet(properties, "NuNHitsTotal",             srslice.nuid.nutothits);
      CopyPropertyIfSet(properties, "NuVertexY",                srslice.nuid.nuvtxy);
      CopyPropertyIfSet(properties, "NuWeightedDirZ",           srslice.nuid.nuwgtdirz);
      CopyPropertyIfSet(properties, "NuNSpacePointsInSphere",   srslice.nuid.nusps);
      CopyPropertyIfSet(properties, "NuEigenRatioInSphere",     srslice.nuid.nueigen);
      CopyPropertyIfSet(properties, "CRLongestTrackDirY",       srslice.nuid.crlongtrkdiry);
      CopyPropertyIfSet(properties, "CRLongestTrackDeflection", srslice.nuid.crlongtrkdef);
      CopyPropertyIfSet(properties, "CRFracHitsInLongestTrack", srslice.nuid.crlongtrkhitfrac);
      CopyPropertyIfSet(properties, "CRNHitsMax",               srslice.nuid.crmaxhits);
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


  void FillSliceCRUMBS(const sbn::CRUMBSResult *crumbs,
                       caf::SRSlice& slice,
                       bool allowEmpty) {
    if (crumbs != nullptr) {
      slice.crumbs_result.score = crumbs->score;
      slice.crumbs_result.ccnumuscore = crumbs->ccnumuscore;
      slice.crumbs_result.ccnuescore = crumbs->ccnuescore;
      slice.crumbs_result.ncscore = crumbs->ncscore;
      slice.crumbs_result.bestscore = crumbs->bestscore;
      slice.crumbs_result.bestid = crumbs->bestid;
      slice.crumbs_result.tpc.crlongtrackhitfrac = crumbs->tpc_CRFracHitsInLongestTrack;
      slice.crumbs_result.tpc.crlongtrackdefl = crumbs->tpc_CRLongestTrackDeflection;
      slice.crumbs_result.tpc.crlongtrackdiry = crumbs->tpc_CRLongestTrackDirY;
      slice.crumbs_result.tpc.crnhitsmax = crumbs->tpc_CRNHitsMax;
      slice.crumbs_result.tpc.nusphereeigenratio = crumbs->tpc_NuEigenRatioInSphere;
      slice.crumbs_result.tpc.nufinalstatepfos = crumbs->tpc_NuNFinalStatePfos;
      slice.crumbs_result.tpc.nutotalhits = crumbs->tpc_NuNHitsTotal;
      slice.crumbs_result.tpc.nuspherespacepoints = crumbs->tpc_NuNSpacePointsInSphere;
      slice.crumbs_result.tpc.nuvertexy = crumbs->tpc_NuVertexY;
      slice.crumbs_result.tpc.nuwgtdirz = crumbs->tpc_NuWeightedDirZ;
      slice.crumbs_result.tpc.stoppingchi2ratio = crumbs->tpc_StoppingChi2CosmicRatio;
      slice.crumbs_result.pds.fmtotalscore = crumbs->pds_FMTotalScore;
      slice.crumbs_result.pds.fmpe = crumbs->pds_FMPE;
      slice.crumbs_result.pds.fmtime = crumbs->pds_FMTime;
      slice.crumbs_result.crt.trackscore = crumbs->crt_TrackScore;
      slice.crumbs_result.crt.hitscore = crumbs->crt_HitScore;
      slice.crumbs_result.crt.tracktime = crumbs->crt_TrackTime;
      slice.crumbs_result.crt.hittime = crumbs->crt_HitTime;
    }
  }

  void FillSliceOpT0Finder(const std::vector<art::Ptr<sbn::OpT0Finder>> &opt0_v,
                           caf::SRSlice &slice)
  {
    if (opt0_v.empty()==false){
      unsigned int nopt0 = opt0_v.size();
      double max_score=-1.; // score of the opt0 object with the highest score 
      double sec_score=-1.; // score of the opt0 object with the 2nd highest score 

      unsigned int max_idx = 0;
      unsigned int sec_idx = 0;

      // fill the default, which is the maximum 
      for (unsigned int i = 0; i < nopt0; i++ ) {
        const sbn::OpT0Finder &thisOpT0 = *opt0_v[i];
        if (thisOpT0.score > max_score){
          max_score = thisOpT0.score;
          max_idx = i;
        }
      }

      const sbn::OpT0Finder &maxOpT0 = *opt0_v[max_idx];
      slice.opt0.tpc    = maxOpT0.tpc;
      slice.opt0.time   = maxOpT0.time;
      slice.opt0.score  = maxOpT0.score;
      slice.opt0.measPE = maxOpT0.measPE;
      slice.opt0.hypoPE = maxOpT0.hypoPE;
      
      // in case there are more matches, find the opt0 object with the second highest score
      // usually this is filled for a slice that is split across two tpcs
      if (nopt0>1){    
        for (unsigned int i = 0; i < nopt0; i++ ) {
          if (i == max_idx) continue;
          const sbn::OpT0Finder &thisOpT0 = *opt0_v[i];
          if (thisOpT0.score > sec_score){
            sec_score = thisOpT0.score;
            sec_idx = i;
          }
        }
        const sbn::OpT0Finder &secOpT0 = *opt0_v[sec_idx];
        slice.opt0_sec.tpc    = secOpT0.tpc;
        slice.opt0_sec.time   = secOpT0.time;
        slice.opt0_sec.score  = secOpT0.score;
        slice.opt0_sec.measPE = secOpT0.measPE;
        slice.opt0_sec.hypoPE = secOpT0.hypoPE;
      }
    }
  }

  void FillSliceBarycenter(const std::vector<art::Ptr<recob::Hit>> &inputHits,
                           const std::vector<art::Ptr<recob::SpacePoint>> &inputPoints,
                           caf::SRSlice &slice)
  {
    unsigned int nHits = inputHits.size();
    double sumCharge = 0.;
    double sumX = 0.; double sumY = 0.; double sumZ = 0.;
    double sumXX = 0.; double sumYY = 0.; double sumZZ = 0.;
    double thisHitCharge, thisPointXYZ[3], chargeCenter[3], chargeWidth[3];

    for ( unsigned int i = 0; i < nHits; i++ ) {
      const recob::Hit &thisHit = *inputHits[i];
      if ( thisHit.SignalType() != geo::kCollection ) continue;
      art::Ptr<recob::SpacePoint> const& thisPoint = inputPoints.at(i);
      if ( !thisPoint ) continue;

      thisHitCharge = thisHit.Integral();
      thisPointXYZ[0] = thisPoint->XYZ()[0];
      thisPointXYZ[1] = thisPoint->XYZ()[1];
      thisPointXYZ[2] = thisPoint->XYZ()[2];

      sumCharge += thisHitCharge;
      sumX += thisPointXYZ[0] * thisHitCharge;
      sumY += thisPointXYZ[1] * thisHitCharge;
      sumZ += thisPointXYZ[2] * thisHitCharge;
      sumXX += thisPointXYZ[0] * thisPointXYZ[0] * thisHitCharge;
      sumYY += thisPointXYZ[1] * thisPointXYZ[1] * thisHitCharge;
      sumZZ += thisPointXYZ[2] * thisPointXYZ[2] * thisHitCharge;
    }

    if( sumCharge != 0 ) {
      chargeCenter[0] = sumX / sumCharge;
      chargeCenter[1] = sumY / sumCharge;
      chargeCenter[2] = sumZ / sumCharge;
      chargeWidth[0] = pow( (sumXX/sumCharge - pow(sumX/sumCharge, 2)), .5);
      chargeWidth[1] = pow( (sumYY/sumCharge - pow(sumY/sumCharge, 2)), .5);
      chargeWidth[2] = pow( (sumZZ/sumCharge - pow(sumZ/sumCharge, 2)), .5);

      slice.charge_center.SetXYZ( chargeCenter[0], chargeCenter[1], chargeCenter[2] ); 
      slice.charge_width.SetXYZ( chargeWidth[0], chargeWidth[1], chargeWidth[2] );
    }
  }

  //......................................................................

  void FillTrackCRTHit(const std::vector<art::Ptr<anab::T0>> &t0match,
                       const std::vector<art::Ptr<sbn::crt::CRTHit>> &hitmatch,
                       bool use_ts0,
                       int64_t CRT_T0_reference_time, // ns, signed
                       double CRT_T1_reference_time, // us
                       caf::SRTrack &srtrack,
                       bool allowEmpty)
  {
    if (t0match.size()) {
      assert(t0match.size() == 1);
      srtrack.crthit.distance = t0match[0]->fTriggerConfidence;
      srtrack.crthit.hit.time = t0match[0]->fTime / 1e3; /* ns -> us */
    }
    if (hitmatch.size()) {
      FillCRTHit(*hitmatch[0], use_ts0, CRT_T0_reference_time, CRT_T1_reference_time, srtrack.crthit.hit, allowEmpty);
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
      srtrack.rangeP.p_proton = range_results[2][0]->range_p;
      assert(track.ID() == range_results[2][0]->trackID);
    }
  }

  void FillPlaneChi2PID(const anab::ParticleID &particle_id, caf::SRTrkChi2PID &srpid) {

    // Assign dummy values.

    srpid.chi2_muon = 0.;
    srpid.chi2_pion = 0.;
    srpid.chi2_kaon = 0.;
    srpid.chi2_proton = 0.;
    srpid.pid_ndof = 0;
    srpid.pida = 0.;

    // Loop over algorithm scores and extract the ones we want.
    // Get the ndof from any chi2 algorithm

    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = particle_id.ParticleIDAlgScores();
    for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){
      anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
      if (AlgScore.fAlgName == "Chi2"){
        if (TMath::Abs(AlgScore.fAssumedPdg) == 13) { // chi2mu
          srpid.chi2_muon = AlgScore.fValue;
          srpid.pid_ndof = AlgScore.fNdf;
        }
        else if (TMath::Abs(AlgScore.fAssumedPdg) == 211) { // chi2pi
          srpid.chi2_pion = AlgScore.fValue;
          srpid.pid_ndof = AlgScore.fNdf;
        }
        else if (TMath::Abs(AlgScore.fAssumedPdg) == 321) { // chi2ka
          srpid.chi2_kaon = AlgScore.fValue;
          srpid.pid_ndof = AlgScore.fNdf;
        }
        else if (TMath::Abs(AlgScore.fAssumedPdg) == 2212) { // chi2pr
          srpid.chi2_proton = AlgScore.fValue;
          srpid.pid_ndof = AlgScore.fNdf;
        }
      }
      else if (AlgScore.fVariableType==anab::kPIDA){
        srpid.pida = AlgScore.fValue;
      }
    }
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
        FillPlaneChi2PID(particle_id, srtrack.chi2pid[plane_id]);
      }
    }
  }

  void FillTrackPlaneCalo(const anab::Calorimetry &calo, 
        const std::vector<art::Ptr<recob::Hit>> &hits,
        bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend, 
        const detinfo::DetectorPropertiesData &dprop,
        caf::SRTrackCalo &srcalo) {

    // Collect info from Calorimetry
    const std::vector<float> &dqdx = calo.dQdx();
    const std::vector<float> &dedx = calo.dEdx();
    const std::vector<float> &pitch = calo.TrkPitchVec();
    const std::vector<float> &rr = calo.ResidualRange();
    const std::vector<geo::Point_t> &xyz = calo.XYZ();
    const std::vector<size_t> &tps = calo.TpIndices();

    srcalo.charge = 0.;
    srcalo.ke = 0.;
    srcalo.nhit = 0;

    float rrmax = !rr.empty() ? *std::max_element(rr.begin(), rr.end()) : 0.;

    for (unsigned i = 0; i < dedx.size(); i++) {
      // Save the points we need to
      if (fill_calo_points && (
          (rrmax - rr[i]) < fillhit_rrstart || // near start
          rr[i] < fillhit_rrend)) { // near end

        // Point information
        caf::SRCaloPoint p;
        p.rr = rr[i];
        p.dqdx = dqdx[i];
        p.dedx = dedx[i];
        p.pitch = pitch[i];
        p.x = xyz[i].x();
        p.y = xyz[i].y();
        p.z = xyz[i].z();

        // lookup the wire -- the Calorimery object makes this
        // __way__ harder than it should be
        for (const art::Ptr<recob::Hit> &h: hits) {
          if (h.key() == tps[i]) {
            p.wire = h->WireID().Wire;
            p.tpc = h->WireID().TPC;
            p.channel = h->Channel();
            p.sumadc = h->SummedADC();
            p.integral = h->Integral();
            p.t = h->PeakTime();
            p.width = h->RMS();
            p.mult = h->Multiplicity();
            p.start = h->StartTick();
            p.end = h->EndTick();
          }
        }

        // Save
        srcalo.points.push_back(p);
      }

      if (dedx[i] > 1000.) continue;
      srcalo.nhit ++;
      srcalo.charge += dqdx[i] * pitch[i]; // ADC
      srcalo.ke += dedx[i] * pitch[i];
    }

    // Sort the points by residual range hi->lo
    std::sort(srcalo.points.begin(), srcalo.points.end(),
      [](const caf::SRCaloPoint &lhs, const caf::SRCaloPoint &rhs) {
        return lhs.rr > rhs.rr;
    });

  }

  void FillTrackScatterClosestApproach(const art::Ptr<sbn::ScatterClosestApproach> closestapproach,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.scatterClosestApproach.mean = closestapproach->mean;
    srtrack.scatterClosestApproach.stdDev = closestapproach->stdDev;
    srtrack.scatterClosestApproach.max = closestapproach->max;
  }

  void FillTrackStoppingChi2Fit(const art::Ptr<sbn::StoppingChi2Fit> stoppingChi2,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.stoppingChi2Fit.pol0Chi2 = stoppingChi2->pol0Chi2;
    srtrack.stoppingChi2Fit.expChi2  = stoppingChi2->expChi2;
    srtrack.stoppingChi2Fit.pol0Fit  = stoppingChi2->pol0Fit;
  }

  void FillTrackDazzle(const art::Ptr<sbn::MVAPID> dazzle,
      caf::SRTrack& srtrack,
      bool allowEmpty)
  {
    srtrack.dazzle.muonScore = dazzle->mvaScoreMap.at(13);
    srtrack.dazzle.pionScore = dazzle->mvaScoreMap.at(211);
    srtrack.dazzle.protonScore = dazzle->mvaScoreMap.at(2212);
    srtrack.dazzle.otherScore = dazzle->mvaScoreMap.at(0);

    srtrack.dazzle.pdg = dazzle->BestPDG();
    srtrack.dazzle.bestScore = dazzle->BestScore();
  }

  void FillTrackCalo(const std::vector<art::Ptr<anab::Calorimetry>> &calos,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend,
                     const geo::GeometryCore *geom, const detinfo::DetectorPropertiesData &dprop,
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
        FillTrackPlaneCalo(calo, hits, fill_calo_points, fillhit_rrstart, fillhit_rrend, dprop, srtrack.calo[plane_id]);
      }
    }

    // Set the plane with the most hits
    //
    // We expect the noise to be lowest at planes 2 -> 0 -> 1, so use this to break ties
    caf::Plane_t bestplane = caf::kUnknown;
    int bestnhit = -1;
    for(int plane: {2, 0, 1}){
      if(srtrack.calo[plane].nhit > bestnhit){
        bestplane = caf::Plane_t(plane);
        bestnhit = srtrack.calo[plane].nhit;
      }
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
                   const art::Ptr<anab::T0> t0,
                   caf::SRPFP& srpfp,
                   bool allowEmpty)
  {
    srpfp.id = particle.Self();
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

      // Pfo Characterisation features
      srpfp.pfochar.setDefault();

      CopyPropertyIfSet(propertiesMap, "LArThreeDChargeFeatureTool_EndFraction",             srpfp.pfochar.chgendfrac);
      CopyPropertyIfSet(propertiesMap, "LArThreeDChargeFeatureTool_FractionalSpread",        srpfp.pfochar.chgfracspread);
      CopyPropertyIfSet(propertiesMap, "LArThreeDLinearFitFeatureTool_DiffStraightLineMean", srpfp.pfochar.linfitdiff);
      CopyPropertyIfSet(propertiesMap, "LArThreeDLinearFitFeatureTool_Length",               srpfp.pfochar.linfitlen);
      CopyPropertyIfSet(propertiesMap, "LArThreeDLinearFitFeatureTool_MaxFitGapLength",      srpfp.pfochar.linfitgaplen);
      CopyPropertyIfSet(propertiesMap, "LArThreeDLinearFitFeatureTool_SlidingLinearFitRMS",  srpfp.pfochar.linfitrms);
      CopyPropertyIfSet(propertiesMap, "LArThreeDOpeningAngleFeatureTool_AngleDiff",         srpfp.pfochar.openanglediff);
      CopyPropertyIfSet(propertiesMap, "LArThreeDPCAFeatureTool_SecondaryPCARatio",          srpfp.pfochar.pca2ratio);
      CopyPropertyIfSet(propertiesMap, "LArThreeDPCAFeatureTool_TertiaryPCARatio",           srpfp.pfochar.pca3ratio);
      CopyPropertyIfSet(propertiesMap, "LArThreeDVertexDistanceFeatureTool_VertexDistance",  srpfp.pfochar.vtxdist);
    }
    if (t0) {
      srpfp.t0 = t0->Time() / 1e3; /* ns -> us */
    }
  }

  void FillHitVars(const recob::Hit& hit,
                   unsigned producer,
                   const recob::SpacePoint& spacepoint,
                   const recob::PFParticle& particle,
                   caf::SRHit& srhit,
                   bool allowEmpty)
  {
    srhit.setDefault();

    srhit.peakTime = hit.PeakTime();
    srhit.RMS = hit.RMS();

    srhit.peakAmplitude = hit.PeakAmplitude();
    srhit.integral = hit.Integral();

    const geo::WireID wire = hit.WireID();
    srhit.cryoID = wire.Cryostat;
    srhit.tpcID = wire.TPC;
    srhit.planeID = wire.Plane;
    srhit.wireID = wire.Wire;
    srhit.spacepoint.XYZ = SRVector3D (spacepoint.XYZ());
    srhit.spacepoint.chisq = spacepoint.Chisq();
    srhit.spacepoint.pfpID = particle.Self();
    srhit.spacepoint.ID = spacepoint.ID();
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

  //......................................................................
  template<class T, class U>
  void CopyPropertyIfSet( const std::map<std::string, T>& props, const std::string& search, U& value )
  {
    auto it = props.find(search);
    if ( it != props.end() ) value = it->second;
  }

} // end namespace
