#pragma once

#include "CAFAna/Core/Cut.h"
#include "SBNAna/Cuts/Cuts.h"
#include "SBNAna/Vars/NueVars.h"

namespace ana {
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

const Cut kPreNueSelND = kFiducialVolumeND && kSlcNuScoreCut && kSlcFlashMatchCut;
const Cut kRecoNueSel = kRecoShower && kShowerEnergyCut;
const Cut kFullNueSel = kNueTrackLenCut && kShowerConvGapCut && kShowerdEdxCut && kShowerDensityCut;

// Specific FD cuts although they're currenly basically the same as for the ND
const Cut kNueFlashScoreFDCut   = kSlcFlashMatchCut;
const Cut kNuePandoraScoreFDCut = kSlcNuScoreCut;
const Cut kRecoShowerFD = kRecoShower && kNueNumShowersCut && kShowerdEdxCut && kShowerConvGapCut && kNueTrackLenCut && kShowerDensityCut && kShowerEnergyCut;
const Cut kNueFD        = kContainedFD && kNueFlashScoreFDCut && kNuePandoraScoreFDCut && kRecoShowerFD;

}
