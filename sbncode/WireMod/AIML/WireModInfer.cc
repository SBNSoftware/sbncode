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

  // loop over the hits and get the X,Y,Z;
  // direction in X,Y,Z;
  // direction in Y,Z relative to the plane;
  // Channel; ThetaXW; & dQ/dx
  std::array<float, N_FEATURES> hit_features;
  auto& [hitX, hitY, hitZ, hitDirX, hitDirY, hitDirZ, hitDirYRel, hitDirZRel, hitChan, hitTheta, hitdQdx]
    = hit_features; // alias for convenience
  for (const auto& old_hit : old_hits)
  {
    // Try the backtracker to get the position
    // If this fails it's a data or overlay hit and so shouldn't be modified
    try
    {
      std::vector<double> hitXYZ = back_tracker->HitToXYZ(*det_clock, old_hit);
      hitX = hitXYZ.at(0);
      hitY = hitXYZ.at(1);
      hitZ = hitXYZ.at(2);
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

    // Now that the hit is known to be back trackable, we can start getting things
    hitChan = old_hit->Channel();

    // Now get the plane angle and rotate the dirs for the plane rel versions
    // and get the ThetaXW
    geo::View_t hitView = old_hit->View();
    geo::TPCID hitTPC = old_hit->WireID(); // WireID inherits from TPCID
    float planeTh = wire_geom->WireAngleToVertical(hitView, hitTPC) - 0.5*M_PI;
    float sinTh = std::sin(planeTh);
    float cosTh = std::cos(planeTh);
    hitDirYRel = hitDirY * cosTh - hitDirZ * sinTh;
    hitDirZRel = hitDirY * sinTh + hitDirZ * cosTh;
    hitTheta = std::atan(hitDirX / hitDirZ);
    hitTheta *= 180.0/M_PI; // model expects degrees /not/ radians
    float pitch = wire_geom->Plane(hitTPC, hitView).WirePitch() / std::abs(hitDirZRel);
    hitdQdx = old_hit->Integral() / pitch;

    // Now that the features are set, infer and construct a new hit
    // Check that there aren't NaNs or other nonsense values
    auto [integral, width] = infer(hit_features);
    if (std::isnan(integral) || std::isnan(width))
    {
      new_hits.push_back(*old_hit);
      continue;
    }
    WMHit new_hit(*old_hit, integral / old_hit->Integral(), width / old_hit->RMS()); 
    //float amplitude = integral / (width * SQRT2PI);
    //float sigma_amplitude = old_hit->SigmaPeakAmplitude() * amplitude / old_hit->PeakAmplitude();
    //float ROISumADC = old_hit->ROISummedADC() * old_hit->Integral() / integral;
    //float HitSumADC = old_hit->HitSummedADC() * old_hit->Integral() / integral;
    //float sigma_integral = old_hit->SigmaIntegral() * integral / old_hit->Integral();
    //recob::Hit new_hit(hitChan,                     // Hit channel unchanged
    //                   old_hit->StartTick(),        // Keep timing info the same
    //                   old_hit->EndTick(),          // Keep timing info the same
    //                   old_hit->PeakTime(),         // Keep timing info the same
    //                   old_hit->SigmaPeakTime(),    // Keep timing info the same
    //                   width,                       // *** The New Width ***
    //                   amplitude,                   // *** The New Amplitude ***
    //                   sigma_amplitude,             // *** Scaled Amplitude Sigma ***
    //                   ROISumADC,                   // *** Scaled ROI ADC Sum ***
    //                   HitSumADC,                   // *** Scaled Hit ADC Sum ***
    //                   integral,                    // *** The New Integral ***
    //                   sigma_integral,              // *** Scaled Integral Sigma ***
    //                   old_hit->Multiplicity(),     // Keep multiplicity
    //                   old_hit->LocalIndex(),       // Keep local index
    //                   old_hit->GoodnessOfFit(),    // Keep goodness of fit
    //                   old_hit->DegreesOfFreedom(), // Keep dof
    //                   hitView,                     // Keep view
    //                   old_hit->SignalType(),       // Keep signal type
    //                   old_hit->WireID());          // Keep wire ID
    new_hits.push_back(static_cast<recob::Hit>(new_hit));
  }

  return new_hits;
}
