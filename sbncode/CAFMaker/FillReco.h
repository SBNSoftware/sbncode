
#ifndef CAF_FILLRECO_H
#define CAF_FILLRECO_H

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "LArReco/TrajectoryMCSFitter.h"
#include "LArReco/TrackMomentumCalculator.h"

#include "sbncode/StandardRecord/SRSlice.h"
#include "sbncode/StandardRecord/StandardRecord.h"

namespace caf
{

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

  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic);

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle& particle,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);

  void FillTrackMCS(const recob::Track& track,
                    const trkf::TrajectoryMCSFitter *mcs_calculator,
                    caf::SRTrack& srtrack,
                    bool allowEmpty = false);

  void FillTrackRangeP(const recob::Track& track,
                     const trkf::TrackMomentumCalculator *range_calculator,
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
}

#endif
