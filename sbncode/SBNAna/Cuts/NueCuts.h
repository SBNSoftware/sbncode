#pragma once

#include "CAFAna/Core/Cut.h"
#include "SBNAna/Vars/NueVars.h"

namespace ana
{
  extern const Cut kRecoShower;
  extern const Cut kNueBasicCut;

  extern const Cut kNueHasTrackCut;
  extern const Cut kNueTrackLenCut;

  extern const Cut kNueNumShowersCut;

  extern const Cut kShowerEnergyCut;
  extern const Cut kShowerdEdxCut;
  extern const Cut kShowerConvGapCut;
  extern const Cut kShowerDensityCut;

  extern const Cut kNueContainedND;
  extern const Cut kNueContainedFD;

  const Cut kNueCut = kRecoShower && kNueBasicCut && (kRecoShower_ConversionGap < 2.1);
}
