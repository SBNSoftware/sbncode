#pragma once

#include "CAFAna/Core/Cut.h"
#include "SBNAna/Vars/NueVars.h"

namespace ana
{
	extern const Cut kRecoShower;
	extern const Cut kNueBasicCut;
	extern const Cut kNueContainedND;
	extern const Cut kNueContainedFD;
	
	const Cut kNueFDCut = kRecoShower && kNueBasicCut && (kRecoShower_ConversionGap < 2.1);	
}
