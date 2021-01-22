#pragma once

#include "CAFAna/Core/MultiVar.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"

namespace ana
{
  // Return the run number
  extern const SpillVar kRun;

  // Return the event number
  extern const SpillVar kEvt;

  extern const SpillMultiVar kCRTHitX;	
  extern const SpillMultiVar kCRTHitY;  
  extern const SpillMultiVar kCRTHitZ;  
  extern const SpillMultiVar kCRTHitPE;  
  extern const SpillMultiVar kCRTHitTime;  
  extern const SpillMultiVar kCRTHitTimeFD;  

  // // Return the slice number
  // extern const SliceVar kSlc;

  // Return spill count
  extern const Var kCounting;
  extern const SpillVar kSpillCounting;

  // Slice verteces
  extern const Var kSlcVtxX;
  extern const Var kSlcVtxY;
  extern const Var kSlcVtxZ;
  extern const Var kSlcNuScore;
  extern const Var kSlcHasFlash;
  extern const Var kSlcFlashScore;

}
