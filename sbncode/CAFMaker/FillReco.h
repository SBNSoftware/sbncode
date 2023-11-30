
#ifndef CAF_FILLRECO_H
#define CAF_FILLRECO_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbnobj/Common/Reco/Stub.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbnanaobj/StandardRecord/SRSlice.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"


namespace caf
{

  void FillStubVars(const sbn::Stub &stub,
                    const art::Ptr<recob::PFParticle> stubpfp,
                    caf::SRStub &srstub,
                    bool allowEmpty = false);

  void FillShowerVars(const recob::Shower& shower,
                      const recob::Vertex* vertex,
                      const std::vector<art::Ptr<recob::Hit>> &hits,
                      const geo::GeometryCore *geom,
                      unsigned producer,
                      caf::SRShower& srshower,
                      bool allowEmpty = false);

  void FillShowerRazzle(const art::Ptr<sbn::MVAPID> razzle,
                        caf::SRShower& srshower,
                        bool allowEmpty = false);

  void FillShowerCosmicDist(const std::vector<art::Ptr<float> >& cosmicDistVec,
                      caf::SRShower& srshower);

  void FillShowerResiduals(const std::vector<art::Ptr<float> >& residuals,
                      caf::SRShower& srshower);

  void FillShowerTrackFit(const sbn::ShowerTrackFit& trackFit,
                      caf::SRShower& srshower);

  void FillShowerDensityFit(const sbn::ShowerDensityFit& densityFit,
                      caf::SRShower& srshower);

  void FillSliceVars(const recob::Slice& slice,
                     const recob::PFParticle *primary,
                     unsigned producer,
                     caf::SRSlice& srslice,
                     bool allowEmpty = false);

  void FillSliceMetadata(const larpandoraobj::PFParticleMetadata *primary_meta,
                        caf::SRSlice &srslice,
                        bool allowEmpty = false);

  void FillSliceVertex(const recob::Vertex *vertex,
                       caf::SRSlice& slice,
                       bool allowEmpty = false);

  void FillSliceCRUMBS(const sbn::CRUMBSResult *crumbs,
                       caf::SRSlice& slice,
                       bool allowEmpty = false);

  void FillSliceBarycenter(const std::vector<art::Ptr<recob::Hit>> &inputHits,
                           const std::vector<art::Ptr<recob::SpacePoint>> &inputPoints,
                           caf::SRSlice &slice);

  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic);

  void FillTrackVars(const recob::Track& track,
                     unsigned producer,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);

  void FillHitVars(const recob::Hit& hit,
                   unsigned producer,
                   const recob::SpacePoint& spacepoint,
                   const recob::PFParticle& particle,
                   caf::SRHit& srhit,
                   bool allowEmpty = false);

  void FillPFPVars(const recob::PFParticle &particle,
                   const recob::PFParticle *primary,
                   const larpandoraobj::PFParticleMetadata *pfpMeta,
                   const art::Ptr<anab::T0> t0,
                   caf::SRPFP& srpfp,
                   bool allowEmpty= false);

  void FillTrackCRTHit(const std::vector<art::Ptr<anab::T0>> &t0match,
                       const std::vector<art::Ptr<sbn::crt::CRTHit>> &hitmatch,
                       bool use_ts0,
                       int64_t CRT_T0_reference_time, // ns, signed
                       double CRT_T1_reference_time, // us
                       caf::SRTrack &srtrack,
                       bool allowEmpty = false);

  void FillTrackCRTTrack(const std::vector<art::Ptr<anab::T0>> &t0match,
                       caf::SRTrack &srtrack,
                       bool allowEmpty = false);

  void FillTrackMCS(const recob::Track& track,
                    const std::array<std::vector<art::Ptr<recob::MCSFitResult>>, 4> &mcs_results,
                    caf::SRTrack& srtrack,
                    bool allowEmpty = false);

  void FillTrackRangeP(const recob::Track& track,
                     const std::array<std::vector<art::Ptr<sbn::RangeP>>, 3> &range_results,
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void FillPlaneChi2PID(const anab::ParticleID &particle_id, caf::SRTrkChi2PID &srpid);
  void FillTrackChi2PID(const std::vector<art::Ptr<anab::ParticleID>> particleIDs,
                        const geo::GeometryCore *geom,
                        caf::SRTrack& srtrack,
                        bool allowEmpty = false);

  void FillTrackPlaneCalo(const anab::Calorimetry &calo, 
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend, 
                     const detinfo::DetectorPropertiesData &dprop,
                     caf::SRTrackCalo &srcalo);

  void FillTrackScatterClosestApproach(const art::Ptr<sbn::ScatterClosestApproach> closestapproach,
                           caf::SRTrack& srtrack,
                           bool allowEmpty = false);

  void FillTrackStoppingChi2Fit(const art::Ptr<sbn::StoppingChi2Fit> stoppingChi2,
                                caf::SRTrack& srtrack,
                                bool allowEmpty = false);

  void FillTrackDazzle(const art::Ptr<sbn::MVAPID> dazzle,
                        caf::SRTrack& srtrack,
                        bool allowEmpty = false);

  void FillTrackCalo(const std::vector<art::Ptr<anab::Calorimetry>> &calos,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend,
                     const geo::GeometryCore *geom, const detinfo::DetectorPropertiesData &dprop,
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void SetNuMuCCPrimary(std::vector<caf::StandardRecord> &recs,
                        std::vector<caf::SRTrueInteraction> &srneutrinos);
  void ApplyNumuCCMatching(std::vector<caf::StandardRecord> &recs,
                           const std::vector<caf::SRTrueInteraction> &srneutrinos,
                           unsigned truth_ind);

  void FillCRTHit(const sbn::crt::CRTHit &hit,
                  bool use_ts0,
                  int64_t CRT_T0_reference_time, // ns, signed
                  double CRT_T1_reference_time, // us
                  caf::SRCRTHit &srhit,
                  bool allowEmpty = false);
  void FillCRTTrack(const sbn::crt::CRTTrack &track,
                  bool use_ts0,
                  caf::SRCRTTrack &srtrack,
                  bool allowEmpty = false);

  void FillOpFlash(const recob::OpFlash &flash,
                  std::vector<recob::OpHit const*> const& hits,
                  int cryo,
                  caf::SROpFlash &srflash,
                  bool allowEmpty = false);
  void FillCRTPMTMatch(const sbn::crt::CRTPMTMatching &match,
		  caf::SRCRTPMTMatch &srmatch,
		  bool allowEmpty = false);

  void FillTPCPMTBarycenterMatch(const sbn::TPCPMTBarycenterMatch *matchInfo,
                           caf::SRSlice& slice);

  template<class T, class U>
  void CopyPropertyIfSet( const std::map<std::string, T>& props, const std::string& search, U& value );
}

#endif
