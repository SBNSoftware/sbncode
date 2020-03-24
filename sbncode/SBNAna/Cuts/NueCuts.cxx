#include "SBNAna/Cuts/NueCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	// Basic reconstruction 
	const Cut kNueBasic(
		[](const caf::SRProxy* sr)
		{
		  return (//sr->trk.cosmic.nshw > 0 && // need a shower
			  //sr->reco.shw[0].energy > 0 && // nothing is terribly wrong
			sr->reco.shw[0].len > 0 );
	}
	);

}
