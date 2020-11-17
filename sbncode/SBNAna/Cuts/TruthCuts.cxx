#include "SBNAna/Cuts/TruthCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	const Cut kIsAntiNu([](const caf::SRSliceProxy* slc){
		if(int(slc->mc.nnu) == 0 ) return false;
		assert( int(slc->mc.nnu) == 1);
		return slc->mc.nu[0].pdg < 0;
	});

	const Cut kIsNu([](const caf::SRSliceProxy* slc){
		if(int(slc->mc.nnu) == 0) return false;
		assert(int(slc->mc.nnu) == 1);
		return slc->mc.nu[0].pdg > 0;
	});

	const Cut kHasNu([](const caf::SRSliceProxy* slc){
		if(int(slc->mc.nnu) == 0) return false;
		assert(int(slc->mc.nnu) == 1);
		return true;
	});


	const Cut kIsNue([](const caf::SRSliceProxy* slc){
		return (int(slc->mc.nnu) == 1 && abs(slc->mc.nu[0].pdg) ==12);
	});

	const Cut kIsNumu([](const caf::SRSliceProxy* slc){
		return (int(slc->mc.nnu) == 1 && abs(slc->mc.nu[0].pdg) ==14);
	});

	const Cut kIsNutau([](const caf::SRSliceProxy* slc){
		return (int(slc->mc.nnu) == 1 && abs(slc->mc.nu[0].pdg) ==16);
	});


	const Cut kIsNC([](const caf::SRSliceProxy* slc){
		if(int(slc->mc.nnu) == 0) return false;
		assert(int(slc->mc.nnu) == 1);
		return !slc->mc.nu[0].iscc;
	});


}
