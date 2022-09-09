#ifndef CAF_FILLTRIGGER_H
#define CAF_FILLTRIGGER_H

#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnanaobj/StandardRecord/SRTrigger.h"
#include "lardataobj/RawData/TriggerData.h"

#include <vector>

namespace caf
{

  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
                   const std::vector<raw::Trigger>& trig_info,
                   std::vector<caf::SRTrigger>& triggerInfo);

}

#endif
