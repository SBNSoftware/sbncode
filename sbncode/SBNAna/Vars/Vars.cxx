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

	// // TO DO: Fix error 'const class caf::Proxy<std::vector<caf::SRCRTHit> > position'
	// const SpillVar kCRTHitX = SIMPLESPILLVAR(crt_hits.position.x);
	// const SpillVar kCRTHitY = SIMPLESPILLVAR(crt_hits.position.y);
	// const SpillVar kCRTHitZ = SIMPLESPILLVAR(crt_hits.position.z);
	// const SpillVar kCRHitPE = SIMPLESPILLVAR(crt_hits.pe);
	// const SpillVar kCRHitTime = SIMPLESPILLVAR(crt_hits.time);

	const Var kCounting = kUnweighted;

	// // For when we have spill beam mode info
	// const SpillCut kIsRHC([](const caf::SRProxy* sr) {return sr->spill.isRHC;});

 	const Var kSlcVtxX([](const caf::SRSliceProxy *slc)
		   {
		     return slc->vertex.x;
		   });

	const Var kSlcVtxY([](const caf::SRSliceProxy *slc)
		   {
		     return slc->vertex.y;
		   });

	const Var kSlcVtxZ([](const caf::SRSliceProxy *slc)
		   {
		     return slc->vertex.z;
		   });

}
