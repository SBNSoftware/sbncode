#include "SBNAna/Vars/TruthVars.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include <cassert>

namespace ana
{
	  const Var kTruthEnergy(
      [](const caf::SRSliceProxy *slc)
      {
      	return ( slc->truth.E );
      }
      );

}