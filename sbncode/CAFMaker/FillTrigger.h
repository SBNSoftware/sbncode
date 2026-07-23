#ifndef CAF_FILLTRIGGER_H
#define CAF_FILLTRIGGER_H

#include "sbncode/CAFMaker/TimeRefShifters.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnanaobj/StandardRecord/SRTrigger.h"
#include "sbnanaobj/StandardRecord/SRSBNDTimingInfo.h"
#include "lardataobj/RawData/TriggerData.h"

namespace caf
{

  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
                   const raw::Trigger& trig_info,
                   caf::SRTrigger& triggerInfo,
                   caf::TimeRefShifter<> const& shifter = {});

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, caf::SRTrigger& triggerInfo);
}

#endif
