#include "CAFAna/Cuts/NumuCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	// Basic reconstruction 
	const Cut kNueBasic(
		[](const caf::SRProxy* sr)
		{
		return (sr->trk.cosmic.nshw > 0 && // need a shower
			sr->reco.shw[0].energy > 0 && // nothing is terribly wrong
			sr->reco.shw[0].len > 0 );
	}
	);

	// True cuts
	const Cut kIsNueCC(
		[](const caf::SRProxy* sr)
		{
			return (sr->mc.nu[0].iscc &&
			abs(sr->mc.nu[0].pdg==12) );
		}
		);

	const Cut kIsNueNC(
		[](const caf::SRProxy* sr)
		{
			return (sr->mc.nu[0].isnc &&
			abs(sr->mc.nu[0].pdg==12) );
		}
		);

}