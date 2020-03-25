#include "SBNAna/Cuts/TruthCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	const Cut kIsAntiNu([](const caf::SRProxy* sr){
		if(int(sr->mc.nnu) == 0 ) return false;
		assert( int(sr->mc.nnu) == 1);
		return sr->mc.nu[0].pdg < 0;
	});

	const Cut kIsNu([](const caf::SRProxy* sr){
		if(int(sr->mc.nnu) == 0) return false;
		assert(int(sr->mc.nnu) == 1);
		return sr->mc.nu[0].pdg > 0;
	});

	const Cut kHasNu([](const caf::SRProxy* sr){
		if(int(sr->mc.nnu) == 0) return false;
		assert(int(sr->mc.nnu) == 1);
		return true;
	});


	const Cut kIsNue([](const caf::SRProxy* sr){
		return (int(sr->mc.nnu) == 1 && abs(sr->mc.nu[0].pdg) ==12);
	});

	const Cut kIsNumu([](const caf::SRProxy* sr){
		return (int(sr->mc.nnu) == 1 && abs(sr->mc.nu[0].pdg) ==14);
	});

	const Cut kIsNutau([](const caf::SRProxy* sr){
		return (int(sr->mc.nnu) == 1 && abs(sr->mc.nu[0].pdg) ==16);
	});


	const Cut kIsNC([](const caf::SRProxy* sr){
		if(int(sr->mc.nnu) == 0) return false;
		assert(int(sr->mc.nnu) == 1);
		return !sr->mc.nu[0].iscc;
	});


}
