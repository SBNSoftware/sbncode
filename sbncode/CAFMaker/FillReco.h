#ifndef CAF_FILLRECO_H
#define CAF_FILLRECO_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesData.h"

// LArSoft includes
#include "larcorealg/Geometry/fwd.h"

#include "larevt/SpaceCharge/SpaceCharge.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/MVAOutput.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "larrecodnn/CVN/func/Result.h"
#include "sbnobj/Common/Reco/Stub.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"
#include "sbnobj/Common/Reco/MVAPID.h"
#include "sbnobj/Common/Reco/CNNScore.h"
#include "sbnobj/Common/Reco/ScatterClosestApproach.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/TPCPMTBarycenterMatch.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "sbnobj/SBND/CRT/CRTSpacePoint.hh"
#include "sbnobj/SBND/CRT/CRTTrack.hh"
#include "sbnobj/Common/CRT/CRTPMTMatching.hh"
#include "sbnobj/Common/CRT/CRTHitT0TaggingInfo.hh"
#include "sbnobj/SBND/Timing/TimingInfo.hh"
#include "sbnobj/SBND/Timing/FrameShiftInfo.hh"

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
                      const geo::WireReadoutGeom& wireReadout,
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

  void FillSliceOpT0Finder(const std::vector<art::Ptr<sbn::OpT0Finder>> &opt0_v,
                           caf::SRSlice &slice);

  void FillSliceBarycenter(const std::vector<art::Ptr<recob::Hit>> &inputHits,
                           const std::vector<art::Ptr<recob::SpacePoint>> &inputPoints,
                           caf::SRSlice &slice);

  /**
   * @brief Fills the results from NuGraph at slice level
   * @param inputHits (pointers to) the hits associated to the slice
   * @param sliceHitsMap maps position of hits in collection input to NuGraph (slice only) to the one input to Pandora (all gaus hits)
   * @param ngFilterResult NuGraph filter result, for each hit
   * @param ngSemanticResult NuGraph semnatic result, for each hit (MIP track, HIP, shower, Michel electron, diffuse activity)
   * @param[out] slice the destination slice object
   *
   * Hits with filter value (`ngFilterResult`) lower than `ng_filter_cut` are counted as background.
   */
  void FillSliceNuGraph(const std::vector<art::Ptr<recob::Hit>> &inputHits,
			const std::vector<unsigned int> &sliceHitsMap,
			const std::vector<art::Ptr<anab::FeatureVector<1>>> &ngFilterResult,
			const std::vector<art::Ptr<anab::FeatureVector<5>>> &ngSemanticResult,
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

  void FillCNNScores(const recob::PFParticle &particle,
                     const sbn::PFPCNNScore *cnnscore,
                     caf::SRPFP& srpfp,
                     bool allowEmpty = false);

  /**
   * @brief Fills the results from NuGraph at slice level
   * @param sliceHitsMap maps position of hits in collection input to NuGraph (slice only) to the one input to Pandora (all gaus hits)
   * @param ngFilterResult NuGraph filter result, for each hit
   * @param ngSemanticResult NuGraph semnatic result, for each hit (MIP track, HIP, shower, Michel electron, diffuse activity)
   * @param pfpHits Vector of hits associated to the PFParticle
   * @param[out] srpfp the destination PFParticle object
   *
   * Hits with filter value (`ngFilterResult`) lower than `ng_filter_cut` are counted as background.
   */
  void FillPFPNuGraph(const std::vector<unsigned int> &sliceHitsMap,
		      const std::vector<art::Ptr<anab::FeatureVector<1>>> &ngFilterResult,
		      const std::vector<art::Ptr<anab::FeatureVector<5>>> &ngSemanticResult,
		      const std::vector<art::Ptr<recob::Hit>> &pfpHits,
		      caf::SRPFP& srpfp,
		      bool allowEmpty= false);

  void FillTrackCRTHit(const std::vector<art::Ptr<anab::T0>> &t0match,
                       const std::vector<art::Ptr<sbn::crt::CRTHit>> &hitmatch,
                       const std::vector<art::Ptr<sbn::crt::CRTHitT0TaggingInfo>> &hitmatchinfo,
                       bool use_ts0,
                       int64_t CRT_T0_reference_time, // ns, signed
                       double CRT_T1_reference_time, // us
                       const std::map<std::pair<int, int>, sim::AuxDetSimChannel> &crtsimchanmap,
                       caf::SRTrack &srtrack,
                       bool allowEmpty = false);

  void FillTrackCRTTrack(const std::vector<art::Ptr<anab::T0>> &t0match,
                       caf::SRTrack &srtrack,
                       bool allowEmpty = false);

  void FillTrackCRTSpacePoint(const anab::T0 &t0match,
                              const art::Ptr<sbnd::crt::CRTSpacePoint> &spacepointmatch,
                              caf::SRTrack &srtrack,
                              bool allowEmpty = false);

  void FillTrackSBNDCRTTrack(const anab::T0 &t0match,
                             const art::Ptr<sbnd::crt::CRTTrack> &trackmatch,
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
                        caf::SRTrack& srtrack,
                        bool allowEmpty = false);

  void FillTrackPlaneCalo(const anab::Calorimetry &calo, 
                     const recob::Track& track,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     const std::vector<const recob::TrackHitMeta*>& thms,
                     bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend, 
                     const detinfo::DetectorPropertiesData &dprop,
                     spacecharge::SpaceCharge const& sce,
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
                     const recob::Track& track,
                     const std::vector<art::Ptr<recob::Hit>> &hits,
                     const std::vector<const recob::TrackHitMeta*>& thms,
                     bool fill_calo_points, float fillhit_rrstart, float fillhit_rrend,
                     const detinfo::DetectorPropertiesData &dprop,
                     spacecharge::SpaceCharge const& sce,
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
                  const std::map<std::pair<int, int>, sim::AuxDetSimChannel> &crtsimchanmap,
                  caf::SRCRTHit &srhit,
                  bool allowEmpty = false);
  void FillCRTTrack(const sbn::crt::CRTTrack &track,
                  bool use_ts0,
                  caf::SRCRTTrack &srtrack,
                  bool allowEmpty = false);

  void FillCRTSpacePoint(const sbnd::crt::CRTSpacePoint &spacepoint,
                         caf::SRCRTSpacePoint &srspacepoint,
                         bool allowEmpty = false);

  void FillSBNDCRTTrack(const sbnd::crt::CRTTrack &track,
                        caf::SRSBNDCRTTrack &srsbndcrttrack,
                        bool allowEmpty = false);

  void FillICARUSOpFlash(const recob::OpFlash &flash,
                  std::vector<recob::OpHit const*> const& hits,
                  int cryo,
                  caf::SROpFlash &srflash,
                  bool allowEmpty = false);

  void FillSBNDOpFlash(const recob::OpFlash &flash,
                  std::vector<recob::OpHit const*> const& hits,
                  int tpc,
                  caf::SROpFlash &srflash,
                  bool allowEmpty = false);
                  
  void FillCRTPMTMatch(const sbn::crt::CRTPMTMatching &match,
		       caf::SRCRTPMTMatch &srmatch,
		       bool allowEmpty = false);

  void FillTPCPMTBarycenterMatch(const sbn::TPCPMTBarycenterMatch *matchInfo,
                           caf::SRSlice& slice);

  void FillCVNScores(const lcvn::Result *cvnResult,
                     caf::SRSlice& slice);

  void FillPFPRazzled(const art::Ptr<sbn::MVAPID> razzled,
                      caf::SRPFP& srpfp,
                      bool allowEmpty = false);

  void FillSBNDFrameShiftInfo(const sbnd::timing::FrameShiftInfo &frame,
                        caf::SRSBNDFrameShiftInfo &srsbndframe,
                        bool allowEmpty = false);

  void FillSBNDTimingInfo(const sbnd::timing::TimingInfo &timing,
                        caf::SRSBNDTimingInfo &srsbndtiming,
                        bool allowEmpty = false);

  template<class T, class U>
  void CopyPropertyIfSet( const std::map<std::string, T>& props, const std::string& search, U& value );
}

#endif
