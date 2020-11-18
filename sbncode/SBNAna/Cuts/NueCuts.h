#pragma once

#include "CAFAna/Core/Cut.h"
#include "SBNAna/Vars/NueVars.h"

namespace ana
{
	extern const Cut kRecoShower;
	extern const Cut kNueBasicCut;
	extern const Cut kNueContainedFD;
	
	const Cut kNueCut = kRecoShower && kNueBasicCut && (kRecoShower_ConversionGap < 2.1);	
}