#ifndef CAF_FILLTRIGGER_H
#define CAF_FILLTRIGGER_H

#include "sbncode/CAFMaker/TimeRefShifters.h"
#include "sbnobj/Common/Trigger/ExtraTriggerInfo.h"
#include "sbnanaobj/StandardRecord/SRTrigger.h"
#include "sbnanaobj/StandardRecord/SRSoftwareTrigger.h"
#include "sbnanaobj/StandardRecord/SRSBNDTimingInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "art/Framework/Principal/Handle.h"
#include "sbndaq-artdaq-core/Obj/SBND/pmtSoftwareTrigger.hh"
#include "sbncode/BeamSpillInfoRetriever/POTTools.h"

#include <vector>

namespace caf
{
  void FillTrigger(const sbn::ExtraTriggerInfo& addltrig_info,
                   const raw::Trigger& trig_info,
                   caf::SRTrigger& triggerInfo,
                   caf::TimeRefShifter<> const& shifter = {});

  void FillTriggerICARUS(const sbn::ExtraTriggerInfo& addltrig_info,
                         caf::SRTrigger& triggerInfo);

  void FillTriggerSBND(caf::SRSBNDTimingInfo& timingInfo, caf::SRTrigger& triggerInfo);

  void FillTriggerEmulation(art::Handle<std::vector<int>> const& monpulsesFlat,
                             art::Handle<std::vector<int>> const& monpulseSizes,
                             art::Handle<int> const& numPairs,
                             art::Handle<bool> const& passedTrig,
                             caf::SRTrigger& triggerInfo);
  void FillSoftwareTriggerSBND(const sbnd::trigger::pmtSoftwareTrigger& softInfo, caf::SRSoftwareTrigger& caf_softInfo);
  
  void FillPTBTriggersSBND(const std::vector<sbn::pot::PTBInfo_t>& ptb_triggers, caf::SRTrigger& triggerInfo);
}

#endif
