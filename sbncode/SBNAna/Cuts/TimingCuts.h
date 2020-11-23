#pragma once

#include "CAFAna/Core/Cut.h"

namespace ana
{
	/// Does the event fall inside the beam spill window
	extern const Cut kInBeamSpill;

	/// Is the event far from the start and end of the spill window
	extern const Cut kInTimingSideband;

}