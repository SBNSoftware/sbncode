
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

  struct SliceData {
    unsigned particle_id_offset;
    unsigned slice_id_offset;
    const recob::PFParticle *primary;
    const larpandoraobj::PFParticleMetadata *primary_meta;
    const anab::T0 *fmatch;
  };

  void FillSliceVars(const recob::Slice& slice,
                     const SliceData sdata,
                     caf::SRSlice& srslice,
                     bool allowEmpty = false);

  bool SelectSlice(const caf::SRSlice &slice, bool cut_clear_cosmic);

  struct TrackData {
    unsigned particle_index_offset;

    // data products
    std::vector<art::Ptr<anab::ParticleID>> particleIDs;
    std::vector<art::Ptr<anab::Calorimetry>> calos;
    std::vector<art::Ptr<recob::Hit>> hits;

    // algorithms
    const trkf::TrajectoryMCSFitter *mcs_calculator;
    const trkf::TrackMomentumCalculator *range_calculator;
    const geo::GeometryCore *geom;
  };

  void FillTrackVars(const recob::Track& track,
                     const recob::PFParticle& particle,
                     const TrackData tdata,
                     caf::SRTrack& srtrk,
                     bool allowEmpty = false);
}

#endif
