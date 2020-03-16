
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
#include "sbncode/LArRecoProducer/Products/RangeP.h"
#include "sbncode/LArRecoProducer/Products/CRTHit.hh"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "sbncode/StandardRecord/SRSlice.h"
#include "sbncode/StandardRecord/StandardRecord.h"


namespace caf
{

  void FillShowerVars(const recob::Shower& shower,
                      caf::SRShower& srshower,
                      bool allowEmpty = false);

  void FillSliceVars(const recob::Slice& slice,
                     const recob::PFParticle *primary,
                     caf::SRSlice& srslice,
                     bool allowEmpty = false);

  void FillSliceFlashMatch(const anab::T0 *fmatch,
                           caf::SRSlice &srslice,
                           bool allowEmpty = false);

  void FillSliceMetadata(const larpandoraobj::PFParticleMetadata *primary_meta,
                        caf::SRSlice &srslice,
                        bool allowEmpty = false);

  void FillSliceTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                      const std::vector<caf::SRTrueInteraction> &srneutrinos,
                      const cheat::ParticleInventoryService &inventory_service,
                      caf::SRSlice &srslice,
                      bool allowEmpty = false);

  void FillSliceVertex(const recob::Vertex *vertex,
                       caf::SRSlice& slice,
                       bool allowEmpty = false);

  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic);

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle& particle,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);

  void FillTrackCRTHit(const std::vector<art::Ptr<sbn::crt::CRTHit>> &hitmatch, 
                       const std::vector<const anab::T0*> &t0match, 
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
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void FillTrackTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                     caf::SRTrack& srtrack,
                     bool allowEmpty = false);

  void SetNuMuCCPrimary(std::vector<caf::StandardRecord> &recs,
                        std::vector<caf::SRTrueInteraction> &srneutrinos);
}

#endif
