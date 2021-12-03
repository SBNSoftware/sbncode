#ifndef CAF_FILLEXPOSURE_H
#define CAF_FILLEXPOSURE_H

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRHeader.h"
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
