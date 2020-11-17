#include "SBNAna/Cuts/NueCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	// Basic reconstruction 
	const Cut kRecoShower(
		[](const caf::SRSliceProxy* slc)
		{
		  return (slc->reco.nshw > 0 && // need a shower
				  slc->reco.shw[0].bestplane_energy > 0 && // nothing is terribly wrong
				  slc->reco.shw[0].len > 0 );
	}
	);

	// Basic reconstruction 
	const Cut kNueBasicCut(
		[](const caf::SRSliceProxy* slc)
		{
		  return (slc->reco.shw[0].bestplane_energy < 250. &&
				  slc->reco.shw[0].bestplane_dEdx < 2.7 &&
				  slc->reco.shw[0].len < 42.);
	}
	);

}
