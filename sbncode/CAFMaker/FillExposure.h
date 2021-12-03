#ifndef CAF_FILLEXPOSURE_H
#define CAF_FILLEXPOSURE_H

#include "canvas/Persistency/Common/Ptr.h"
#include "sbnanaobj/StandardRecord/SRBNBInfo.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include <vector>

namespace caf
{
  void FillExposure(const std::vector<art::Ptr<sbn::BNBSpillInfo> > bnb_spill_info,
		    std::vector<caf::SRBNBInfo>& BNBInfo,
		    double& subRunPOT);
  
}

#endif
