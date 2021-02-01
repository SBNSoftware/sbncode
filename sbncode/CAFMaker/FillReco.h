
#ifndef CAF_FILLRECO_H
#define CAF_FILLRECO_H

#include <array>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "sbnobj/Common/Reco/RangeP.h"
#include "sbnobj/Common/Reco/ShowerSelectionVars.h"
#include "sbnobj/Common/CRT/CRTHit.hh"
#include "sbnobj/Common/CRT/CRTTrack.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbnanaobj/StandardRecord/SRSlice.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"


namespace caf
{

  void FillShowerVars(const recob::Shower& shower,
                      const recob::PFParticle &particle,
                      const recob::Vertex* vertex,
                      const recob::PFParticle *primary,
                      const std::vector<art::Ptr<recob::Hit>> &hits,
                      const geo::GeometryCore *geom,
                      unsigned producer,
                      caf::SRShower& srshower,
                      bool allowEmpty = false);

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

  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic);

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle& particle,
                     const recob::PFParticle *primary,
                     unsigned producer,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);

  void FillTrackCRTHit(const std::vector<art::Ptr<anab::T0>> &t0match,
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
                     const std::array<std::vector<art::Ptr<sbn::RangeP>>, 2> &range_results,
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void FillTrackChi2PID(const std::vector<art::Ptr<anab::ParticleID>> particleIDs,
                        const geo::GeometryCore *geom,
                        caf::SRTrack& srtrack,
                        bool allowEmpty = false);

  void FillTrackCalo(const std::vector<art::Ptr<anab::Calorimetry>> &calos,
                     const geo::GeometryCore *geom,
                     const std::array<float, 3> &calo_constants,
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void SetNuMuCCPrimary(std::vector<caf::StandardRecord> &recs,
                        std::vector<caf::SRTrueInteraction> &srneutrinos);
  void ApplyNumuCCMatching(std::vector<caf::StandardRecord> &recs,
                           const std::vector<caf::SRTrueInteraction> &srneutrinos,
                           unsigned truth_ind);

  void FillCRTHit(const sbn::crt::CRTHit &hit,
                  bool use_ts0,
                  caf::SRCRTHit &srhit,
                  bool allowEmpty = false);
  void FillCRTTrack(const sbn::crt::CRTTrack &track,
                  bool use_ts0,
                  caf::SRCRTTrack &srtrack,
                  bool allowEmpty = false);
}

#endif
