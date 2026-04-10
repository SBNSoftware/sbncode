/**
 * WireModInfer.cc - WaireMLod interface for inference at SBN
 *
 * Loads weights for the Julia WaireMLod network from binary file
 * and runs model inference on hits to modify them to be more ``data-like``
 */

#include "sbncode/WireMod/AIML/WireModInfer.hh"

std::vector<recob::Hit> sys::WaireMLod::produceNew(const std::vector<recob::Hit> old_hits, 
                                                   cheat::BackTrackerService* back_tracker,
                                                   detinfo::DetectorClocksData* det_clock,
                                                   geo::WireReadoutGeom const* wire_geom) const
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
      new_hits.push_back(old_hit);
      continue;
    }

    // Now that the hit is known to be back trackable, get the sim channel
    hitChan = old_hit.Channel();
    art::Ptr<sim::SimChannel> simChan = back_tracker->FindSimChannel(hitChan);

    // Sum up the charge for the hit
    // Take the sume of the charge on the sim channel btwn -1 and +1 sigma on the peak time
    // Clamp to zero to avoid out-of-bounds issues
    // while here also get the charge avged position at each tick
    float hitQ = 0.0;
    int start_tdc = det_clock->TPCTick2TDC(old_hit.PeakTimeMinusRMS(1));
    int end_tdc   = det_clock->TPCTick2TDC(old_hit.PeakTimePlusRMS(1));
    if (start_tdc < 0) start_tdc = 0;
    if (end_tdc   < 0) end_tdc   = 0;
    std::vector<float> ideXs;
    std::vector<float> ideYs;
    std::vector<float> ideZs;
    for (int tdc = start_tdc; tdc < end_tdc - 1; ++tdc)
    {
      float tdcQ = simChan->Charge(tdc);
      // The IDEsBetween method isn't available in v10_02_00 of lardataobj...
      // for (auto const& [_, ides] : simChan->IDEsBetween(tdc, tdc + 1))
      // {
      //   if (ides.size() == 0)
      //     continue;
      //   float ideX = 0;
      //   float ideY = 0;
      //   float ideZ = 0;
      //   for (auto const& ide : ides)
      //   {
      //     ideX += ide.x * ide.numElectrons;
      //     ideY += ide.y * ide.numElectrons;
      //     ideZ += ide.z * ide.numElectrons;
      //   }
      //   ideX /= tdcQ;
      //   ideY /= tdcQ;
      //   ideZ /= tdcQ;
      //   ideXs.push_back(ideX);
      //   ideYs.push_back(ideY);
      //   ideZs.push_back(ideZ);
      // }
      float ideX = 0;
      float ideY = 0;
      float ideZ = 0;
      // I am pretty sure this is what we want in absence of IDEsBetween
      for (auto const& ide : simChan->TrackIDsAndEnergies(tdc, tdc + 1))
      {
        ideX += ide.x * ide.numElectrons;
        ideY += ide.y * ide.numElectrons;
        ideZ += ide.z * ide.numElectrons;
      }
      ideX /= tdcQ;
      ideY /= tdcQ;
      ideZ /= tdcQ;
      ideXs.push_back(ideX);
      ideYs.push_back(ideY);
      ideZs.push_back(ideZ);
      hitQ += tdcQ;
    }

    // if there are no ides somehow, don't modify the hit and use the return the old one
    // if there's only ides at one tick we can't get a direction so don't modify it
    if (ideXs.size() < 2)
    {
      new_hits.push_back(old_hit);
      continue;
    }

    // Take the direction to be going from the start X,Y,Z to the last X,Y,Z
    // This should be a reasonable approximation if there's minimal scatter on the hit time scale
    float deltaX = ideXs.back() - ideXs.front();
    float deltaY = ideXs.back() - ideXs.front();
    float deltaZ = ideXs.back() - ideXs.front();
    float displacement = std::sqrt(deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ);
    hitDirX = deltaX / displacement;
    hitDirY = deltaY / displacement;
    hitDirZ = deltaZ / displacement;
    hitdQdx =   hitQ / displacement;

    // Now get the plane angle and rotate the dirs for the plane rel versions
    // and get the ThetaXW
    geo::View_t hitView = old_hit.View();
    geo::TPCID hitTPC = old_hit.WireID(); // WireID inherits from TPCID
    float planeTh = wire_geom->WireAngleToVertical(hitView, hitTPC);
    float sinTh = std::sin(planeTh);
    float cosTh = std::sin(planeTh);
    hitDirYRel = hitDirY * cosTh - hitDirZ * sinTh;
    hitDirZRel = hitDirY * sinTh + hitDirZ * cosTh;
    hitTheta = std::atan(hitDirX / hitDirZ);

    // Now that the features are set, infer and construct a new hit
    auto [integral, width] = infer(hit_features);
    float amplitude = integral / (width * SQRT2PI);
    float sigma_amplitude = old_hit.SigmaPeakAmplitude() * amplitude / old_hit.PeakAmplitude();
    float ROISumADC = old_hit.ROISummedADC() * old_hit.Integral() / integral;
    float HitSumADC = old_hit.HitSummedADC() * old_hit.Integral() / integral;
    float sigma_integral = old_hit.SigmaIntegral() * integral / old_hit.Integral();
    recob::Hit new_hit(hitChan,                    // Hit channel unchanged
                       start_tdc,                  // Keep timing info the same
                       end_tdc,                    // Keep timing info the same
                       old_hit.PeakTime(),         // Keep timing info the same
                       old_hit.SigmaPeakTime(),    // Keep timing info the same
                       width,                      // *** The New Width ***
                       amplitude,                  // *** The New Amplitude ***
                       sigma_amplitude,            // *** Scaled Amplitude Sigma ***
                       ROISumADC,                  // *** Scaled ROI ADC Sum ***
                       HitSumADC,                  // *** Scaled Hit ADC Sum ***
                       integral,                   // *** The New Integral ***
                       sigma_integral,             // *** Scaled Integral Sigma ***
                       old_hit.Multiplicity(),     // Keep multiplicity
                       old_hit.LocalIndex(),       // Keep local index
                       old_hit.GoodnessOfFit(),    // Keep goodness of fit
                       old_hit.DegreesOfFreedom(), // Keep dof
                       hitView,                    // Keep view
                       old_hit.SignalType(),       // Keep signal type
                       old_hit.WireID());          // Keep wire ID
    new_hits.push_back(new_hit);
  }

  return new_hits;
}
