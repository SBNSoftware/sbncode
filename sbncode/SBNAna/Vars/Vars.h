#pragma once

#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"

namespace ana
{
  // Return the run number
  extern const SpillVar kRun;

  // Return the event number
  extern const SpillVar kEvt;

  // extern const SpillVar kCRTHitX;	
  // extern const SpillVar kCRTHitY;  
  // extern const SpillVar kCRTHitZ;  
  // extern const SpillVar kCRTHitPE;  
  // extern const SpillVar kCRTHitTime;  

  // Return the slice number
  //  extern const SliceVar kSlc;

  // Return event count
  extern const SpillVar kCountingSpill;
  extern const Var kCounting;

  // Slice verteces
  extern const Var kSlcVtxX;
  extern const Var kSlcVtxY;
  extern const Var kSlcVtxZ;

}
