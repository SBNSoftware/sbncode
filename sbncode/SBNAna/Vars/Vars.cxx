#include "SBNAna/Vars/Vars.h"

#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <vector>

namespace ana
{

	// const Var kRun = SIMPLEVAR(hdr.run);
	// const Var kEvt = SIMPLEVAR(hdr.evt);
	// const Var kSlc = SIMPLEVAR(hdr.subevt);

	const Var kCounting([](const caf::SRSliceProxy *slc)
	{
		return 1;
	});

	// // For when we have spill beam mode info
	// const Cut kIsRHC([](const caf::SRSliceProxy* sr) {return sr->spill.isRHC;});

}
