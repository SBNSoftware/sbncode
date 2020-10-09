#pragma once

#include "CAFAna/Core/Var.h"
#include "CAFAna/Core/Cut.h"

namespace ana
{
	// Return the run number
	extern const Var kRun;

	// Return the event number
	extern const Var kEvt;
	
	// Return the slice number
	extern const Var kSlc;

	// Return event count
	extern const Var kCounting;

	// Select beam mode
	extern const Cut kIsRHC;

}