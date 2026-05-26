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
                                                   const detinfo::DetectorPropertiesData* det_prop,
                                                   const geo::WireReadoutGeom* wire_geom) const
{
  std::vector<recob::Hit> new_hits;

  // loop over the hits and get the X,Y,Z;
  // direction in X,Y,Z;
  // direction in Y,Z relative to the plane;
  // Channel; ThetaXW; & dQ/dx
  //
  // Caller-facing array is N_RAW = 11 primitives. The 15-feature network
  // input (sin/cos θXW, sin/cos 2θ_plane, block_is_odd) is built inside
  // WaireMLod::infer via derive_features.
  std::array<float, sys::N_RAW> hit_features;
  auto& [hitX, hitY, hitZ, hitDirX, hitDirY, hitDirZ, hitDirYRel, hitDirZRel, hitChan, hitTheta, hitdQdx]
    = hit_features; // alias for convenience
  for (const auto& old_hit : old_hits)
  {
    // Try the backtracker to get the position
    // If this fails it's a data or overlay hit and so shouldn't be modified
    float nElec = 0;
    try
    {
      std::vector<double> hitXYZ = back_tracker->HitToXYZ(*det_clock, old_hit);
      hitX = hitXYZ.at(0);
      hitY = hitXYZ.at(1);
      hitZ = hitXYZ.at(2);
      // Get simIDEs for hit
      // Loop over the simIDEs and add up the charge 
      std::vector<const sim::IDE*> simIDEs = back_tracker->HitToSimIDEs_Ps(*det_clock, old_hit);
      for (auto const& ide : simIDEs)
        nElec += ide->numElectrons;
    }
    catch (...)
    {
      // use the old hit unmodified
      new_hits.push_back(*old_hit);
      continue;
    }

    // Get MCParticle for hit
    TruthMatchUtils::G4ID particleID = TruthMatchUtils::TrueParticleID(*det_clock, old_hit, true);
    if (not TruthMatchUtils::Valid(particleID))
    {
      new_hits.push_back(*old_hit);
      continue;
    }
    const simb::MCTrajectory trajectory = particles->TrackIdToParticle_P(particleID)->Trajectory();
    float hitTrajMinDistSq = std::numeric_limits<float>::max();
    hitDirX = std::numeric_limits<float>::min();
    hitDirY = std::numeric_limits<float>::min();
    hitDirZ = std::numeric_limits<float>::min();
    for (auto const& tp : trajectory)
    {
      float hitTrajDistSq = std::pow(tp.first.X() - hitX, 2)
                          + std::pow(tp.first.Y() - hitY, 2)
                          + std::pow(tp.first.Z() - hitZ, 2);
      if (hitTrajDistSq < hitTrajMinDistSq)
      {
        hitTrajMinDistSq = hitTrajDistSq;
        hitDirX = tp.second.Px() / tp.second.Vect().Mag();
        hitDirY = tp.second.Py() / tp.second.Vect().Mag();
        hitDirZ = tp.second.Pz() / tp.second.Vect().Mag();
      }
    }
    if (hitDirX == std::numeric_limits<float>::min() ||
        hitDirY == std::numeric_limits<float>::min() ||
        hitDirZ == std::numeric_limits<float>::min() )
      {
        new_hits.push_back(*old_hit);
        continue;
      }

    // Now that the hit is known to be back trackable, we can start getting things
    hitChan = old_hit->Channel();

    // Now get the plane angle and rotate the dirs for the plane rel versions
    // and get the ThetaXW
    geo::View_t hitView = old_hit->View();
    geo::TPCID hitTPC = old_hit->WireID(); // WireID inherits from TPCID
    float planeTh = wire_geom->WireAngleToVertical(hitView, hitTPC); // Gets plane.ThetaZ() in [0, pi]
    planeTh = (planeTh < M_PI/2) ? planeTh : planeTh - M_PI;         // Want [-pi/2, pi/2]
    if (std::min({std::abs(planeTh - PLANE_THETA_1),
                  std::abs(planeTh - PLANE_THETA_2),
                  std::abs(planeTh - PLANE_THETA_3)}) > 1e-5f)
      throw std::runtime_error("Plane Angle Mismatch! Got "+std::to_string(planeTh)+" for view "+std::to_string(hitView));
    float sinTh = std::sin(planeTh);
    float cosTh = std::cos(planeTh);
    hitDirYRel = hitDirY * cosTh - hitDirZ * sinTh;
    hitDirZRel = hitDirY * sinTh + hitDirZ * cosTh;
    float vertTh = wire_geom->WireAngleToVertical(hitView, hitTPC) - 0.5*M_PI;
    float cosG = std::abs(hitDirY * std::sin(vertTh) + hitDirZ * std::cos(vertTh));
    hitTheta = std::atan(hitDirX / std::max(cosG, 1e-5f)) / DEG2RAD; // model expects degrees /not/ radians
    float pitch = wire_geom->Plane(hitTPC, hitView).WirePitch() / std::max(cosG, 1e-5f);
    float gain = (old_hit->Channel() % CHANNELS_PER_BLOCK < 2304) ? GAIN[0]
               : (old_hit->Channel() % CHANNELS_PER_BLOCK < 8064) ? GAIN[1]
               :                                                    GAIN[2];
    hitdQdx = (nElec * gain) / pitch;

    // Now that the features are set, infer and construct a new hit
    // Check that there aren't NaNs or other nonsense values
    auto [integral, width] = infer(hit_features);
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
