#include "SBNAna/Vars/Vars.h"

#include "CAFAna/Core/Utilities.h"
#include "StandardRecord/Proxy/SRProxy.h"

#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <vector>

namespace ana
{

	const SpillVar kRun = SIMPLESPILLVAR(hdr.run);
	const SpillVar kEvt = SIMPLESPILLVAR(hdr.evt);
  //	const Var kSlc = SIMPLEVAR(hdr.subevt);

        const Var kCounting = kUnweighted;

	// // For when we have spill beam mode info
	// const SpillCut kIsRHC([](const caf::SRProxy* sr) {return sr->spill.isRHC;});

}
