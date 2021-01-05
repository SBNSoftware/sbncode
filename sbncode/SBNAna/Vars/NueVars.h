#pragma once
// Specific variables for the Nue event selection.

#include "CAFAna/Core/Var.h"

namespace ana
{
  // Reconstructed variables
  extern const Var kLargestRecoShowerIdx;
  extern const Var kRecoShower_BestEnergy;
  extern const Var kRecoShower_BestdEdx;
  extern const Var kRecoShower_ConversionGap;
  extern const Var kRecoShower_Density;
  extern const Var kRecoShower_Energy;
  extern const Var kRecoShower_Length;
  extern const Var kRecoShower_OpenAngle;
  extern const Var kRecoShower_StartX;
  extern const Var kRecoShower_StartY;
  extern const Var kRecoShower_StartZ;
  extern const Var kRecoShower_EndX;
  extern const Var kRecoShower_EndY;
  extern const Var kRecoShower_EndZ;

  extern const Var kRecoShower_TruePdg;
  extern const Var kLongestTrackTruePdg;
  extern const Var kRecoShowers_EnergyCut;

  extern const Var kLongestTrackIdx;
  extern const Var kLongestTrackLength;
} // namespace

