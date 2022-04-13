#ifndef CAF_FILLEXPOSURE_H
#define CAF_FILLEXPOSURE_H

#include "sbnanaobj/StandardRecord/SRBNBInfo.h"
#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"
#include "sbnanaobj/StandardRecord/SRNuMIInfo.h"
#include "sbnobj/Common/POTAccounting/NuMISpillInfo.h"

#include <vector>

namespace caf
{
  void FillExposure(const std::vector<sbn::BNBSpillInfo>& bnb_spill_info,
		    std::vector<caf::SRBNBInfo>& BNBInfo,
		    double& subRunPOT);

  caf::SRBNBInfo makeSRBNBInfo(sbn::BNBSpillInfo const& info);

  void FillExposureNuMI(const std::vector<sbn::NuMISpillInfo>& numi_spill_info,
		    std::vector<caf::SRNuMIInfo>& NuMIInfo,
		    double& subRunPOT);

  
}

#endif
