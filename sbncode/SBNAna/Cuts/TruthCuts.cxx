#include "SBNAna/Cuts/TruthCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

	const Cut kIsAntiNu([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return slc->truth.pdg < 0;
	});

	const Cut kIsNu([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return slc->truth.pdg > 0;
	});

	const Cut kHasNu([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0;
	});


	const Cut kIsNue([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 12;
	});

	const Cut kIsNumu([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 14;
	});

	const Cut kIsNutau([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 16;
	});


	const Cut kIsNC([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return !slc->truth.iscc;
	});


}
