#pragma once

#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"

namespace ana
{
  // Return the run number
  extern const SpillVar kRun;

  // Return the event number
  extern const SpillVar kEvt;
	
  // Return the slice number
  //  extern const SliceVar kSlc;

  // Return event count
  extern const Var kCounting;

  // Select beam mode
  extern const SpillCut kIsRHC;

}
