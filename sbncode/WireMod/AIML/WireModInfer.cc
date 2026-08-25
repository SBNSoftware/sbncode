/**
 * WireModInfer.cc - WaireMLod interface for inference at SBN
 *
 * Loads weights for the Julia WaireMLod network from binary file
 * and runs model inference on hits to modify them to be more ``data-like``
 */

#include <math.h>

#include "sbncode/WireMod/AIML/WireModInfer.hh"

std::vector<recob::Hit> sys::WaireMLod::produceNew(const std::vector<art::Ptr<recob::Hit>> old_hits, 
                                                   const cheat::BackTrackerService* back_tracker,
                                                   const cheat::ParticleInventoryService* particles,
                                                   const detinfo::DetectorClocksData* det_clock,
                                                   const geo::WireReadoutGeom* wire_geom) const
{
  std::vector<recob::Hit> new_hits;

  // Per-hit features (truth-level: backtracker position, MC-trajectory direction,
  // sim-electron dQ/dx) come from featuresForHit — the SAME method the shard
  // producer uses, so what the model sees here is identical to what's binned in
  // the closure shard. A non-back-trackable hit (data/overlay) or one without a
  // usable MC trajectory has valid==false and is passed through unmodified.
  for (const auto& old_hit : old_hits)
  {
    HitFeatures hf = featuresForHit(old_hit, back_tracker, particles, det_clock, wire_geom);
    if (not hf.valid)
    {
      new_hits.push_back(*old_hit);
      continue;
    }

    // Infer and construct a new hit; guard against NaN/nonsense outputs.
    auto [integral, width] = infer(hf.nraw);
    if (std::isnan(integral) || std::isnan(width))
    {
      new_hits.push_back(*old_hit);
      continue;
    }
    WMHit new_hit(*old_hit, integral / old_hit->Integral(), width / old_hit->RMS());
    new_hits.push_back(static_cast<recob::Hit>(new_hit));
  }

  return new_hits;
}
